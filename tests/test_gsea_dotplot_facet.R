#!/usr/bin/env Rscript
## ============================================================================
## test_gsea_dotplot_facet.R
##
## Test suite for gsea_dotplot_facet()
## Tests the shared dotplot contract:
## 1. Top-N selection is by significance WITHIN each Up/Down facet
## 2. Base dots use fill-based NES rendering (shape 21), not colour-only points
## 3. Significant pathways are indicated by black outline overlay
## ============================================================================

library(ggplot2)

source("scripts/custom_minimal_theme.R")
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R")

suppressPackageStartupMessages({
  library(methods)
  library(clusterProfiler)
})

test_passed <- function(name) cat(sprintf("[PASS] %s\n", name))
test_failed <- function(name, reason) cat(sprintf("[FAIL] %s: %s\n", name, reason))

make_mock_gsea_facet <- function() {
  genes <- paste0("Gene", seq_len(400))
  ranks <- sort(rnorm(length(genes), sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  ids <- paste0("MOCK_PATH_", seq_len(8))
  sets <- lapply(ids, function(pid) sample(genes, size = 30))
  names(sets) <- ids

  result_df <- data.frame(
    ID = ids,
    Description = c(
      "Up low padj", "Up high NES", "Up second low padj", "Up third low padj",
      "Down low padj", "Down high |NES|", "Down second low padj", "Down third low padj"
    ),
    setSize = 30L,
    NES = c(1.2, 3.8, 0.9, 0.7, -1.1, -4.1, -0.8, -0.6),
    pvalue = c(0.002, 0.40, 0.006, 0.02, 0.003, 0.35, 0.008, 0.03),
    p.adjust = c(0.004, 0.45, 0.009, 0.03, 0.005, 0.40, 0.01, 0.04),
    qvalue = c(0.004, 0.45, 0.009, 0.03, 0.005, 0.40, 0.01, 0.04),
    rank = sample(50:300, 8),
    leading_edge = rep("tags=30%, list=50%", 8),
    core_enrichment = vapply(sets, function(gs) paste(gs[1:5], collapse = "/"), character(1)),
    stringsAsFactors = FALSE,
    row.names = ids
  )

  new("gseaResult",
      result = result_df,
      geneSets = sets,
      geneList = ranks,
      readable = FALSE,
      keytype = "SYMBOL",
      organism = "Mus musculus",
      params = list(exponent = 1))
}

cat("=== Mock-Based Tests for gsea_dotplot_facet() ===\n")
mock_obj <- make_mock_gsea_facet()
output_dir <- "tests/output/gsea_dotplot_facet"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p <- tryCatch(
  gsea_dotplot_facet(
    mock_obj,
    showCategory = 2,
    padj_cutoff = 0.01,
    title = "Facet contract test"
  ),
  error = function(e) { cat(sprintf("[FAIL] Plot build: %s\n", e$message)); NULL }
)

if (!is.null(p) && inherits(p, "gg")) {
  ## Test 1: top-by-significance within each facet, not top-by-|NES|
  selected <- p$data
  up_labels <- as.character(selected$Description[selected$Direction == "Up"])
  down_labels <- as.character(selected$Description[selected$Direction == "Down"])

  expected_up <- c("Up low padj", "Up second low padj")
  expected_down <- c("Down low padj", "Down second low padj")

  if (setequal(up_labels, expected_up) && setequal(down_labels, expected_down)) {
    test_passed("Facet selection uses top p.adjust within each direction")
  } else {
    test_failed("Facet selection",
                sprintf("Up={%s}; Down={%s}",
                        paste(up_labels, collapse = ", "),
                        paste(down_labels, collapse = ", ")))
  }

  ## Test 2: base rendering uses fill aesthetic with shape 21
  layer1 <- p$layers[[1]]
  has_fill <- !is.null(layer1$mapping$fill) || !is.null(p$mapping$fill)
  if (identical(layer1$geom$objname, "point") && has_fill && identical(layer1$aes_params$shape, 21)) {
    test_passed("Base layer uses fill-mapped shape-21 points")
  } else {
    test_failed("Base layer contract", "Expected geom_point(shape=21) with fill mapping")
  }

  ## Test 3: outline layer counts padj < threshold among selected rows
  built <- ggplot_build(p)$data
  expected_outlines <- sum(selected$qvalue < 0.01)
  actual_outlines <- if (length(built) >= 2) nrow(built[[2]]) else 0L

  if (actual_outlines == expected_outlines) {
    test_passed(sprintf("Outline layer marks %d selected significant pathways", actual_outlines))
  } else {
    test_failed("Outline count",
                sprintf("Expected %d outlines, got %d", expected_outlines, actual_outlines))
  }

  ggsave(file.path(output_dir, "facet_contract.pdf"), p, width = 10, height = 9)
  cat("Saved: tests/output/gsea_dotplot_facet/facet_contract.pdf\n")
}
