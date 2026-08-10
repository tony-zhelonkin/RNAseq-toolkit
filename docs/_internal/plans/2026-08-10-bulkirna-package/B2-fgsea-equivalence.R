#!/usr/bin/env Rscript
# B2 numerical-equivalence diagnostic — NOT a package test.
#
# clusterProfiler::GSEA() was already `by = "fgsea"` under the hood
# (scripts/GSEA/GSEA_processing/run_gsea.R:181), so going direct to
# fgsea::fgseaMultilevel() must be the same computation. This script runs the
# OLD run_gsea() and the NEW gs_test() on the committed real-symbol fixture,
# joins on pathway, and prints the max absolute difference in NES, pval and
# padj plus any pathway present in one and not the other.
#
# Run from the package root:
#   docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
#     -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
#     -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 \
#     Rscript docs/_internal/plans/2026-08-10-bulkirna-package/B2-fgsea-equivalence.R

suppressPackageStartupMessages({
  library(msigdbr)
})

COLLECTION <- "H"
SPECIES    <- "Mus musculus"

de <- readRDS("tests/fixtures/de_real_symbols.rds")
cat(sprintf("Fixture: %d genes\n", nrow(de)))

## ---- OLD path: scripts/ + clusterProfiler --------------------------------
suppressMessages(source("scripts/GSEA/GSEA_processing/run_gsea.R"))
old_obj <- suppressMessages(suppressWarnings(
  run_gsea(de, rank_metric = "t", species = SPECIES, collection = COLLECTION)
))
old <- as.data.frame(old_obj@result)
cat(sprintf("OLD run_gsea(): %d pathways\n", nrow(old)))

## ---- NEW path: bulkiRNA gs_ranks() + gs_test() ---------------------------
suppressMessages(devtools::load_all(".", quiet = TRUE))

m2g <- msigdbr(species = SPECIES, collection = COLLECTION)
m2g <- unique(as.data.frame(m2g)[, c("gs_name", "gene_symbol")])
sets <- split(m2g$gene_symbol, m2g$gs_name)
db <- structure(
  sets,
  pathway_names = stats::setNames(names(sets), names(sets)),
  database = "Hallmark", species = SPECIES, gene_id_type = "symbol",
  class = "gs_db"
)

ranks <- gs_ranks(de, metric = "t")
new <- suppressWarnings(gs_test(ranks, db, contrast = "fixture"))
cat(sprintf("NEW gs_test():  %d pathways\n", nrow(new)))

## ---- Compare --------------------------------------------------------------
only_old <- setdiff(old$ID, new$pathway_id)
only_new <- setdiff(new$pathway_id, old$ID)
cat(sprintf("Pathways only in OLD: %d%s\n", length(only_old),
            if (length(only_old)) paste0(" [", paste(only_old, collapse = ", "), "]") else ""))
cat(sprintf("Pathways only in NEW: %d%s\n", length(only_new),
            if (length(only_new)) paste0(" [", paste(only_new, collapse = ", "), "]") else ""))

j <- merge(
  old[, c("ID", "NES", "pvalue", "p.adjust")],
  as.data.frame(new)[, c("pathway_id", "stat", "p_value", "padj")],
  by.x = "ID", by.y = "pathway_id"
)
cat(sprintf("Joined on %d pathways\n\n", nrow(j)))

d_nes  <- abs(j$NES - j$stat)
d_p    <- abs(j$pvalue - j$p_value)
d_padj <- abs(j$p.adjust - j$padj)

cat(sprintf("max |dNES|  = %.3e   (at %s)\n", max(d_nes),  j$ID[which.max(d_nes)]))
cat(sprintf("max |dp|    = %.3e   (at %s)\n", max(d_p),    j$ID[which.max(d_p)]))
cat(sprintf("max |dpadj| = %.3e   (at %s)\n", max(d_padj), j$ID[which.max(d_padj)]))
cat(sprintf("max relative |dp| (p > 0) = %.3e\n",
            max((d_p / pmax(j$pvalue, .Machine$double.xmin))[j$pvalue > 0])))

cat("\nLargest 5 NES differences:\n")
print(utils::head(j[order(-d_nes), c("ID", "NES", "stat", "pvalue", "p_value")], 5))

cat("\nPlanted signal, both engines:\n")
planted <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_MYC_TARGETS_V1")
print(j[j$ID %in% planted, ])
