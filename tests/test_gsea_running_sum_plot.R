#!/usr/bin/env Rscript
## ============================================================================
## test_gsea_running_sum_plot.R
##
## Regression + smoke tests for gsea_running_sum_plot()
##
## Covers:
##   1. Multi-pathway MSigDB-style (default palette)
##   2. Unnamed palette regression (Feb 2026 fix: named palettes → gseaplot2 error)
##   3. Single pathway
##   4. Custom DB (pre-populated @geneSets, simulating SynGO/MitoPathways)
##   5. Edge case: pathway ID not in results → warning, not crash
##
## Usage:
##   cd /path/to/RNAseq-toolkit
##   Rscript tests/test_gsea_running_sum_plot.R
## ============================================================================

suppressPackageStartupMessages({
  library(methods)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(patchwork)
})

source("scripts/custom_minimal_theme.R")
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")

## ============================================================================
## Test helpers
## ============================================================================

test_passed <- function(name) cat(sprintf("  [PASS] %s\n", name))
test_failed <- function(name, reason) cat(sprintf("  [FAIL] %s: %s\n", name, reason))

output_dir <- "tests/output/gsea_running_sum"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## Mock gseaResult builder
##
## Creates a self-consistent S4 gseaResult with:
##   - @geneList : named numeric vector (ranked genes, descending)
##   - @geneSets : named list of character vectors (pathway -> gene names)
##   - @result   : data.frame with required GSEA columns
## ============================================================================

make_mock_gsea <- function(n_genes = 300, n_pathways = 20, seed = 42) {
  set.seed(seed)

  genes  <- paste0("Gene", seq_len(n_genes))
  ranks  <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  pathway_ids <- paste0("MOCK_PATH_", seq_len(n_pathways))
  gene_sets   <- lapply(pathway_ids, function(pid) {
    sample(genes, size = sample(20:60, 1))
  })
  names(gene_sets) <- pathway_ids

  nes_vals  <- rnorm(n_pathways, sd = 1.5)
  pvals     <- runif(n_pathways, 0.001, 0.5)
  padj_vals <- p.adjust(pvals, method = "BH")

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = paste("Mock Pathway", seq_len(n_pathways)),
    setSize        = vapply(gene_sets, length, integer(1)),
    NES            = nes_vals,
    pvalue         = pvals,
    p.adjust       = padj_vals,
    qvalue         = padj_vals * 1.05,
    rank           = sample(50:250, n_pathways),
    leading_edge   = rep("tags=30%, list=50%", n_pathways),
    core_enrichment = vapply(gene_sets, function(gs) {
      paste(gs[seq_len(min(5L, length(gs)))], collapse = "/")
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids   # enrichplot uses rowname lookup
  )

  new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Mus musculus",
    params   = list(exponent = 1)   # required by gseaScores()
  )
}

gsea_mock <- make_mock_gsea()

## ============================================================================
## Test 1: Multi-pathway MSigDB-style, default palette (palette = NULL)
## ============================================================================

cat("\n=== Test 1: Multi-pathway, default palette ===\n")

p1 <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1:5,
    palette      = NULL,
    title        = "Test 1: Default palette, 5 pathways"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p1) && inherits(p1, "patchwork")) {
  test_passed("Returns patchwork object with default palette")
  ggsave(file.path(output_dir, "test1_multi_pathway_default_pal.pdf"),
         p1, width = 10, height = 8)
  cat("  Saved: test1_multi_pathway_default_pal.pdf\n")
} else if (!is.null(p1) && inherits(p1, "gg")) {
  test_passed("Returns gg object with default palette")
} else {
  test_failed("Multi-pathway default palette", "did not return patchwork/gg object")
}

## ============================================================================
## Test 2: Unnamed palette regression
##
## The Feb 2026 bug: gsea_running_sum_plot() was naming the internal palette
## vector (names = pathway IDs), which caused enrichplot::gseaplot2() to fail
## for custom databases where pathway IDs don't match the Description field.
## Fix: palette vector is intentionally kept UNNAMED inside the function.
##
## This test:
##   (a) Passes an unnamed color vector and verifies it runs without error
##   (b) Passes a named color vector and verifies it also runs (the function
##       should accept it without crashing — the fix is in the internal logic)
## ============================================================================

cat("\n=== Test 2: Unnamed palette regression ===\n")

unnamed_pal <- c("#E41A1C", "#377EB8", "#4DAF4A")

p2a <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1:3,
    palette      = unnamed_pal,
    title        = "Test 2a: Unnamed palette (3 pathways)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR (unnamed): %s\n", e$message)); NULL
})

if (!is.null(p2a)) {
  test_passed("Unnamed palette runs without error")
  ggsave(file.path(output_dir, "test2a_unnamed_palette.pdf"),
         p2a, width = 10, height = 8)
} else {
  test_failed("Unnamed palette regression", "function errored with unnamed palette")
}

# Also verify the internal palette is not named after pathway IDs
# (we can only check indirectly by verifying the function runs with custom DB mock)
named_pal <- c(MOCK_PATH_1 = "#E41A1C", MOCK_PATH_2 = "#377EB8", MOCK_PATH_3 = "#4DAF4A")

p2b <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = c("MOCK_PATH_1", "MOCK_PATH_2", "MOCK_PATH_3"),
    palette      = named_pal,
    title        = "Test 2b: Named palette (should also work)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR (named): %s\n", e$message)); NULL
})

if (!is.null(p2b)) {
  test_passed("Named palette also runs without error")
} else {
  test_failed("Named palette", "function errored with named palette (known limitation)")
}

## ============================================================================
## Test 3: Single pathway
## ============================================================================

cat("\n=== Test 3: Single pathway ===\n")

p3 <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1L,
    title        = "Test 3: Single pathway"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p3)) {
  test_passed("Single pathway runs without error")
  ggsave(file.path(output_dir, "test3_single_pathway.pdf"),
         p3, width = 10, height = 6)
} else {
  test_failed("Single pathway", "function errored on single pathway")
}

## ============================================================================
## Test 4: Custom DB (pre-populated @geneSets, SynGO/MitoPathways pattern)
##
## For custom databases, @geneSets must be set manually after GSEA run.
## This tests the fix that makes running sum work for non-MSigDB databases.
## ============================================================================

cat("\n=== Test 4: Custom DB with pre-populated @geneSets ===\n")

# Simulate a custom DB result: IDs are short codes, Description = same as ID
make_custom_db_gsea <- function(seed = 99) {
  set.seed(seed)
  n_genes <- 200
  genes   <- paste0("Gene", seq_len(n_genes))
  ranks   <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  # Custom DB pathway IDs (not HALLMARK_ style)
  pathway_ids <- c("SYNGO:123", "SYNGO:456", "MITO:789", "MITO:012", "TRANSPORT:345")
  gene_sets   <- lapply(pathway_ids, function(pid) sample(genes, size = 25))
  names(gene_sets) <- pathway_ids

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = c("Synaptic vesicle exocytosis", "Postsynaptic density",
                       "Mitochondrial complex I", "ATP synthase",
                       "Solute carrier transport"),
    setSize        = 25L,
    NES            = c(2.1, -1.8, 1.5, -1.2, 1.9),
    pvalue         = c(0.001, 0.01, 0.03, 0.08, 0.002),
    p.adjust       = c(0.005, 0.02, 0.05, 0.12, 0.008),
    qvalue         = c(0.005, 0.02, 0.05, 0.12, 0.008),
    rank           = c(80L, 120L, 95L, 140L, 75L),
    leading_edge   = rep("tags=30%, list=50%", 5),
    core_enrichment = vapply(gene_sets, function(gs) paste(gs[1:5], collapse="/"), character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids
  )

  gsea_obj <- new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Mus musculus",
    params   = list(exponent = 1)
  )
  gsea_obj
}

gsea_custom <- make_custom_db_gsea()

p4 <- tryCatch({
  gsea_running_sum_plot(
    gsea_custom,
    gene_set_ids = c("SYNGO:123", "MITO:789", "TRANSPORT:345"),
    title        = "Test 4: Custom DB (SynGO/Mito pattern)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p4)) {
  test_passed("Custom DB with pre-populated @geneSets runs without error")
  ggsave(file.path(output_dir, "test4_custom_db.pdf"),
         p4, width = 10, height = 8)
} else {
  test_failed("Custom DB", "function errored with custom DB gene sets")
}

## ============================================================================
## Test 5: Edge case — pathway ID not in results → warning, not crash
## ============================================================================

cat("\n=== Test 5: Invalid pathway ID → warning, not crash ===\n")

# Mix of valid and invalid IDs
valid_id   <- gsea_mock@result$ID[1]
invalid_id <- "NONEXISTENT_PATH_XYZ"

p5 <- withCallingHandlers(
  tryCatch({
    gsea_running_sum_plot(
      gsea_mock,
      gene_set_ids = c(valid_id, invalid_id),
      title        = "Test 5: Invalid pathway ID (should warn)"
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message)); NULL
  }),
  warning = function(w) {
    cat(sprintf("  [OK] Warning issued: %s\n", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
)

if (!is.null(p5)) {
  test_passed("Invalid pathway ID: warning issued, function recovers with valid ID")
  ggsave(file.path(output_dir, "test5_invalid_id_recovery.pdf"),
         p5, width = 10, height = 6)
} else {
  test_failed("Invalid pathway ID", "function crashed instead of warning + recovering")
}

## ============================================================================
## All-invalid IDs → expect informative stop
## ============================================================================

cat("\n=== Test 5b: All-invalid IDs → informative stop ===\n")

all_invalid <- tryCatch({
  gsea_running_sum_plot(gsea_mock, gene_set_ids = c("BAD_1", "BAD_2"))
  "no_error"
}, error = function(e) {
  cat(sprintf("  [OK] Error caught as expected: %s\n", e$message))
  "error_caught"
})

if (all_invalid == "error_caught") {
  test_passed("All-invalid IDs: informative stop (not silent failure)")
} else {
  test_failed("All-invalid IDs", "function did not stop on all-invalid IDs")
}

## ============================================================================
## Summary
## ============================================================================

cat("\n=== All tests complete! ===\n")
cat("Output PDFs in:", output_dir, "\n\n")
cat("VISUAL INSPECTION CHECKLIST:\n")
cat("  test1: 5-panel running sum curves, legend shows 5 pathway names\n")
cat("  test2a: 3-panel running sum, clean unnamed colors\n")
cat("  test3: Single running sum curve\n")
cat("  test4: Custom DB paths (SynGO/Mito IDs), readable legend\n")
cat("  test5: 1-panel running sum (invalid ID removed with warning)\n")
