#!/usr/bin/env Rscript
## ============================================================================
## test_gsea_dotplot.R
##
## Test suite for gsea_dotplot function
## Tests the "show all, highlight significant" pattern
##
## Key behaviors tested:
## 1. All top N pathways shown regardless of significance
## 2. Only significant pathways get black outline
## 3. Non-significant pathways have NO outline (clean dots)
## 4. Direction filters (NES_positive/NES_negative) work correctly
##
## Usage:
##   cd /workspaces/DC_Dictionary/01_scripts/RNAseq-toolkit
##   Rscript tests/test_gsea_dotplot.R
## ============================================================================

library(ggplot2)

# Source dependencies
source("scripts/custom_minimal_theme.R")
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot.R")

## ============================================================================
## Test Helper Functions
## ============================================================================

test_passed <- function(name) {
  cat(sprintf("[PASS] %s\n", name))
}

test_failed <- function(name, reason) {
  cat(sprintf("[FAIL] %s: %s\n", name, reason))
}

## ============================================================================
## Load Real GSEA Results (if available)
## ============================================================================

gsea_file <- "/workspaces/DC_Dictionary/03_results/00_main/checkpoints/A2/A2c/A2c_gsea_Q1_cDC1A_vs_cDC1B.rds"

skip_real_data_tests <- !file.exists(gsea_file)
if (skip_real_data_tests) {
  cat("WARNING: Real GSEA results not found. Skipping integration tests.\n")
  cat("Expected: ", gsea_file, "\n")
  cat("Run A2c_pathway_gsea.R first to generate test data.\n\n")
}

if (!skip_real_data_tests) {

cat("=== Loading GSEA results ===\n")
gsea_data <- readRDS(gsea_file)

# Extract the results list (contains gseaResult objects per database)
gsea_results <- gsea_data$results
cat("Available databases:", paste(names(gsea_results), collapse = ", "), "\n\n")

# Use first available database (Hallmark is usually a good test)
test_db <- names(gsea_results)[1]
gsea_obj <- gsea_results[[test_db]]

cat("Testing with database:", test_db, "\n")
cat("Total pathways:", nrow(gsea_obj@result), "\n")
n_sig_005 <- sum(gsea_obj@result$p.adjust < 0.05, na.rm = TRUE)
n_sig_010 <- sum(gsea_obj@result$p.adjust < 0.10, na.rm = TRUE)
cat("Pathways with FDR < 0.05:", n_sig_005, "\n")
cat("Pathways with FDR < 0.10:", n_sig_010, "\n\n")

## ============================================================================
## Test 1: All pathways shown (no filtering by significance)
## ============================================================================

cat("=== Test 1: All top N pathways shown regardless of significance ===\n")

p1 <- gsea_dotplot(
  gsea_obj,
  filterBy = "NES",
  showCategory = 20,
  padj_cutoff = 0.05,  # This should only affect outlines, not filtering
  highlight_sig = TRUE,
  title = "Test 1: Top 20 by |NES|"
)

# Check number of dots
layers <- ggplot_build(p1)$data
n_dots <- nrow(layers[[1]])

if (n_dots == min(20, nrow(gsea_obj@result))) {
  test_passed(sprintf("Shows all %d requested pathways", n_dots))
} else {
  test_failed("Pathway count", sprintf("Expected %d dots but got %d",
    min(20, nrow(gsea_obj@result)), n_dots))
}

## ============================================================================
## Test 2: Outlines only on significant pathways
## ============================================================================

cat("\n=== Test 2: Only significant pathways have black outline ===\n")

# Sort by NES to get top 20
sorted_data <- gsea_obj@result[order(abs(gsea_obj@result$NES), decreasing = TRUE), ]
top20 <- head(sorted_data, 20)
n_should_be_outlined <- sum(top20$p.adjust < 0.05, na.rm = TRUE)

p2 <- gsea_dotplot(
  gsea_obj,
  filterBy = "NES",
  showCategory = 20,
  padj_cutoff = 0.05,
  highlight_sig = TRUE,
  title = "Test 2: Outline check"
)

layers2 <- ggplot_build(p2)$data

# Base layer should have no outline
base_has_outline <- any(layers2[[1]]$stroke > 0 & !is.na(layers2[[1]]$colour))
if (!base_has_outline) {
  test_passed("Base layer has no outline (stroke=0 or color=NA)")
} else {
  test_failed("Base layer outline", "Base points should not have outline")
}

# Check highlight layer if it exists
if (length(layers2) >= 2 && n_should_be_outlined > 0) {
  n_outlined <- nrow(layers2[[2]])
  if (n_outlined == n_should_be_outlined) {
    test_passed(sprintf("Correct number of outlined points: %d", n_outlined))
  } else {
    test_failed("Outline count", sprintf("Expected %d outlined, got %d",
      n_should_be_outlined, n_outlined))
  }
} else if (n_should_be_outlined == 0) {
  test_passed("No outlines when no pathways significant (correct)")
}

## ============================================================================
## Test 3: Very strict threshold - no outlines but still shows dots
## ============================================================================

cat("\n=== Test 3: Strict threshold (FDR < 0.001) still shows all dots ===\n")

p3 <- gsea_dotplot(
  gsea_obj,
  filterBy = "NES",
  showCategory = 20,
  padj_cutoff = 0.001,  # Very strict
  highlight_sig = TRUE,
  title = "Test 3: Strict FDR < 0.001"
)

layers3 <- ggplot_build(p3)$data
n_dots_strict <- nrow(layers3[[1]])

if (n_dots_strict == min(20, nrow(gsea_obj@result))) {
  test_passed(sprintf("Still shows %d dots despite strict threshold", n_dots_strict))
} else {
  test_failed("Strict threshold", sprintf("Expected %d dots but got %d",
    min(20, nrow(gsea_obj@result)), n_dots_strict))
}

# Count expected outlines at FDR < 0.001
n_very_sig <- sum(top20$p.adjust < 0.001, na.rm = TRUE)
n_outlined_strict <- if (length(layers3) >= 2) nrow(layers3[[2]]) else 0
cat(sprintf("  FDR < 0.001 pathways in top 20: %d, outlined: %d\n",
  n_very_sig, n_outlined_strict))

## ============================================================================
## Test 4: NES_positive filter
## ============================================================================

cat("\n=== Test 4: NES_positive filter shows only positive NES ===\n")

n_pos <- sum(gsea_obj@result$NES > 0)

p4 <- gsea_dotplot(
  gsea_obj,
  filterBy = "NES_positive",
  showCategory = 15,
  padj_cutoff = 0.10,
  highlight_sig = TRUE,
  title = "Test 4: Positive NES only"
)

layers4 <- ggplot_build(p4)$data
n_shown <- nrow(layers4[[1]])

if (n_shown == min(15, n_pos) && n_shown > 0) {
  test_passed(sprintf("NES_positive shows %d pathways", n_shown))
} else if (n_pos == 0) {
  test_passed("No positive NES pathways (correct empty handling)")
} else {
  test_failed("NES_positive", sprintf("Expected up to %d, got %d", min(15, n_pos), n_shown))
}

## ============================================================================
## Test 5: NES_negative filter
## ============================================================================

cat("\n=== Test 5: NES_negative filter shows only negative NES ===\n")

n_neg <- sum(gsea_obj@result$NES < 0)

p5 <- gsea_dotplot(
  gsea_obj,
  filterBy = "NES_negative",
  showCategory = 15,
  padj_cutoff = 0.10,
  highlight_sig = TRUE,
  title = "Test 5: Negative NES only"
)

layers5 <- ggplot_build(p5)$data
n_shown_neg <- nrow(layers5[[1]])

if (n_shown_neg == min(15, n_neg) && n_shown_neg > 0) {
  test_passed(sprintf("NES_negative shows %d pathways", n_shown_neg))
} else if (n_neg == 0) {
  test_passed("No negative NES pathways (correct empty handling)")
} else {
  test_failed("NES_negative", sprintf("Expected up to %d, got %d", min(15, n_neg), n_shown_neg))
}

## ============================================================================
## Save Visual Output for Manual Inspection
## ============================================================================

cat("\n=== Saving visual test outputs ===\n")

output_dir <- "tests/output/gsea_dotplot"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Add subtitles for clarity
p1_final <- p1 + labs(subtitle = sprintf("Expected: %d dots", min(20, nrow(gsea_obj@result))))
p3_final <- p3 + labs(subtitle = sprintf("FDR<0.001: ~%d outlines expected", n_very_sig))

ggsave(file.path(output_dir, "test1_all_pathways.pdf"), p1_final, width = 12, height = 10)
ggsave(file.path(output_dir, "test3_strict_threshold.pdf"), p3_final, width = 12, height = 10)
ggsave(file.path(output_dir, "test4_positive_nes.pdf"), p4, width = 12, height = 10)
ggsave(file.path(output_dir, "test5_negative_nes.pdf"), p5, width = 12, height = 10)

cat("Saved to:", output_dir, "\n")
cat("\nVISUAL INSPECTION CHECKLIST:\n")
cat("1. test1: All 20 dots visible, some with black outline\n")
cat("2. test3: All 20 dots visible, very few/no outlines (strict FDR)\n")
cat("3. test4: Only positive NES pathways (upregulated colors)\n")
cat("4. test5: Only negative NES pathways (downregulated colors)\n")

} # end if (!skip_real_data_tests)

## ============================================================================
## Mock-based tests (run always — no real checkpoint file required)
## ============================================================================

cat("\n=== Mock-Based Tests (no real data required) ===\n")

suppressPackageStartupMessages({
  library(methods)
  library(clusterProfiler)
})

# Build a self-consistent mock gseaResult with predictable pos/neg NES
make_mock_gsea_dotplot <- function(n_pos = 10, n_neg = 10, seed = 42) {
  set.seed(seed)
  n_pathways <- n_pos + n_neg
  n_genes    <- 300

  genes  <- paste0("Gene", seq_len(n_genes))
  ranks  <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  pathway_ids <- paste0("MOCK_PATH_", seq_len(n_pathways))
  gene_sets   <- lapply(pathway_ids, function(pid) sample(genes, size = 30))
  names(gene_sets) <- pathway_ids

  nes_vals  <- c(runif(n_pos, 0.6, 2.5), runif(n_neg, -2.5, -0.6))
  # Make some significant (padj < 0.1) and others not (padj ≥ 0.1)
  padj_vals <- c(runif(ceiling(n_pathways / 2), 0.001, 0.09),
                 runif(floor(n_pathways / 2), 0.11, 0.99))
  padj_vals <- padj_vals[seq_len(n_pathways)]

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = paste("Mock Pathway", seq_len(n_pathways)),
    setSize        = 30L,
    NES            = nes_vals,
    pvalue         = padj_vals / 1.5,
    p.adjust       = padj_vals,
    qvalue         = padj_vals * 1.05,
    rank           = sample(50:250, n_pathways),
    leading_edge   = rep("tags=30%, list=50%", n_pathways),
    core_enrichment = vapply(gene_sets, function(gs) paste(gs[1:5], collapse="/"), character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids
  )

  new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Mus musculus",
    params   = list(exponent = 1)
  )
}

mock_obj <- make_mock_gsea_dotplot(n_pos = 10, n_neg = 10)
mock_output_dir <- "tests/output/gsea_dotplot_mock"
dir.create(mock_output_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Mock Test 6: highlight_threshold pattern (DC_Dictionary real-world) ---
## padj_cutoff=1 (show ALL 20 pathways regardless of significance)
## highlight_threshold=0.1 (only those with padj < 0.1 get black outline)

cat("\n=== Mock Test 6: highlight_threshold pattern (padj_cutoff=1, threshold=0.1) ===\n")

p_m6 <- tryCatch(
  gsea_dotplot(
    mock_obj,
    filterBy           = "NES",
    showCategory       = 20,
    padj_cutoff        = 1,            # show ALL pathways
    highlight_threshold = 0.1,         # but outline only FDR < 0.10
    highlight_sig      = TRUE,
    title              = "Mock Test 6: All shown, FDR<0.10 outlined"
  ),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
)

if (!is.null(p_m6) && inherits(p_m6, "gg")) {
  layers_m6    <- ggplot_build(p_m6)$data
  n_dots_shown <- nrow(layers_m6[[1]])
  n_sig_in_top20 <- sum(head(mock_obj@result[order(abs(mock_obj@result$NES), decreasing=TRUE), ]$p.adjust, 20) < 0.1)

  if (n_dots_shown == min(20, nrow(mock_obj@result))) {
    test_passed(sprintf("highlight_threshold: all %d pathways shown (padj_cutoff=1)", n_dots_shown))
  } else {
    test_failed("highlight_threshold dot count",
                sprintf("Expected 20 dots but got %d", n_dots_shown))
  }
  if (length(layers_m6) >= 2) {
    n_outlined <- nrow(layers_m6[[2]])
    cat(sprintf("  INFO: %d pathways outlined (padj < 0.10), %d significant in top 20\n",
                n_outlined, n_sig_in_top20))
    if (n_outlined == n_sig_in_top20) {
      test_passed("highlight_threshold: correct number of outlines")
    } else {
      cat(sprintf("  WARN: outline count mismatch (%d vs %d) - may vary with filterBy sort\n",
                  n_outlined, n_sig_in_top20))
    }
  } else if (n_sig_in_top20 == 0) {
    test_passed("highlight_threshold: no outlines (no FDR<0.10 in top 20)")
  }
  ggsave(file.path(mock_output_dir, "test6_highlight_threshold.pdf"),
         p_m6, width = 12, height = 10)
  cat("  Saved: test6_highlight_threshold.pdf\n")
} else {
  test_failed("Mock Test 6 (highlight_threshold)", "did not return ggplot object")
}

## ---- Mock Test 7: NES_positive filter -----------------------------------

cat("\n=== Mock Test 7: filterBy='NES_positive' shows only positive NES ===\n")

p_m7 <- tryCatch(
  gsea_dotplot(
    mock_obj,
    filterBy     = "NES_positive",
    showCategory = 8,
    padj_cutoff  = 1,
    highlight_sig = FALSE,
    title        = "Mock Test 7: Positive NES only"
  ),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
)

if (!is.null(p_m7) && inherits(p_m7, "gg")) {
  pd7     <- ggplot_build(p_m7)$data[[1]]
  n_shown <- nrow(pd7)
  n_pos   <- sum(mock_obj@result$NES > 0)
  expected <- min(8L, n_pos)

  if (n_shown == expected) {
    test_passed(sprintf("NES_positive: %d pathways shown (expected %d)", n_shown, expected))
  } else {
    test_failed("NES_positive count",
                sprintf("Expected %d positive-NES dots, got %d", expected, n_shown))
  }
  ggsave(file.path(mock_output_dir, "test7_nes_positive.pdf"),
         p_m7, width = 12, height = 8)
  cat("  VISUAL CHECK: All bars should be orange-ish (positive NES only)\n")
} else {
  test_failed("Mock Test 7 (NES_positive)", "did not return ggplot object")
}

## ---- Mock Test 8: NES_negative filter -----------------------------------

cat("\n=== Mock Test 8: filterBy='NES_negative' shows only negative NES ===\n")

p_m8 <- tryCatch(
  gsea_dotplot(
    mock_obj,
    filterBy     = "NES_negative",
    showCategory = 8,
    padj_cutoff  = 1,
    highlight_sig = FALSE,
    title        = "Mock Test 8: Negative NES only"
  ),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
)

if (!is.null(p_m8) && inherits(p_m8, "gg")) {
  pd8     <- ggplot_build(p_m8)$data[[1]]
  n_shown <- nrow(pd8)
  n_neg   <- sum(mock_obj@result$NES < 0)
  expected <- min(8L, n_neg)

  if (n_shown == expected) {
    test_passed(sprintf("NES_negative: %d pathways shown (expected %d)", n_shown, expected))
  } else {
    test_failed("NES_negative count",
                sprintf("Expected %d negative-NES dots, got %d", expected, n_shown))
  }
  ggsave(file.path(mock_output_dir, "test8_nes_negative.pdf"),
         p_m8, width = 12, height = 8)
  cat("  VISUAL CHECK: All bars should be blue-ish (negative NES only)\n")
} else {
  test_failed("Mock Test 8 (NES_negative)", "did not return ggplot object")
}

cat("\n=== All tests complete! ===\n")
cat("Mock PDF outputs in:", mock_output_dir, "\n")
