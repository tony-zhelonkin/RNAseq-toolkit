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

if (!file.exists(gsea_file)) {
  cat("WARNING: Real GSEA results not found. Skipping integration tests.\n")
  cat("Expected: ", gsea_file, "\n")
  cat("Run A2c_pathway_gsea.R first to generate test data.\n\n")
  quit(status = 0)
}

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

cat("\n=== All tests complete! ===\n")
