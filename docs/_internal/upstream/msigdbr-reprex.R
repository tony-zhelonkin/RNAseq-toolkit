#!/usr/bin/env Rscript
# Reproducible example for the msigdbr ortholog-cache bug.
#
# The bug is order-dependent, so each mode MUST run in its own fresh R session.
#
#   for m in reactome_only hallmark_first symmetry_rev workaround; do
#     Rscript msigdbr-reprex.R $m
#   done
#
# In a container (there is no R on the host):
#
#   docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
#     -v /path/to/.msigdb-cache:/cache -v "$PWD":/out -w /out \
#     scdock-r-dev:v0.5.11 \
#     bash -c 'for m in reactome_only hallmark_first symmetry_rev workaround; do
#                Rscript /out/msigdbr-reprex.R $m; done'
#
# Measured with msigdbr 26.1.0 / babelgene 22.9 / R 4.5.3. Only the first line
# and the workaround are correct:
#
#   reactome_only   Reactome = 10762
#   hallmark_first  H = 4393   then Reactome =  3688   <-- 7074 genes lost
#   symmetry_rev    Reactome = 10762  then H =  3688   <-- 705 genes lost
#   workaround      H = 4393   then Reactome = 10762   <-- cache dropped
#
# 3688 is the intersection of the two collections' mapped gene spaces. Whichever
# collection is queried second is reduced to that intersection, in either order.

suppressPackageStartupMessages(library(msigdbr))

n_genes <- function(...) {
  length(unique(msigdbr(db_species = "HS", species = "mouse", ...)$gene_symbol))
}

REACTOME <- list(collection = "C2", subcollection = "CP:REACTOME")
mode <- commandArgs(TRUE)[1]

if (identical(mode, "reactome_only")) {
  cat(sprintf("reactome_only   Reactome = %d\n", do.call(n_genes, REACTOME)))

} else if (identical(mode, "hallmark_first")) {
  h <- n_genes(collection = "H")
  cat(sprintf("hallmark_first  H = %d   then Reactome = %d\n",
              h, do.call(n_genes, REACTOME)))

} else if (identical(mode, "symmetry_rev")) {
  r <- do.call(n_genes, REACTOME)
  cat(sprintf("symmetry_rev    Reactome = %d  then H = %d\n",
              r, n_genes(collection = "H")))

} else if (identical(mode, "workaround")) {
  h <- n_genes(collection = "H")
  # Drop the species-keyed ortholog cache so the next query rebuilds it.
  key <- paste0("orthologs", babelgene::species("mouse")$taxon_id)
  env <- asNamespace("msigdbr")$pkg_env
  if (exists(key, envir = env, inherits = FALSE)) rm(list = key, envir = env)
  cat(sprintf("workaround      H = %d   then Reactome = %d\n",
              h, do.call(n_genes, REACTOME)))

} else {
  stop("`mode` must be one of: reactome_only, hallmark_first, symmetry_rev, workaround")
}
