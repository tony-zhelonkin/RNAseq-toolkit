#!/usr/bin/env Rscript
## ============================================================================
## test_volcano_plots.R
##
## Comprehensive tests for volcano plot functions to ensure correct behavior
## across various edge cases, especially dashed line alignment with color boundaries
##
## Run this script to verify volcano plot functions work correctly
## ============================================================================

library(testthat)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Source the volcano plot functions
source("scripts/DE/plot_standard_volcano.R")
source("scripts/DE/volcano_helpers.R")

## ============================================================================
## Test Data Generators
## ============================================================================

#' Generate synthetic DE results for testing
#'
#' @param n_genes Total number of genes
#' @param n_sig Number of significant genes (by FDR)
#' @param fdr_cutoff FDR cutoff for significance
#' @param seed Random seed for reproducibility
#' @return Data frame with logFC, P.Value, adj.P.Val
generate_test_de_results <- function(n_genes = 1000,
                                     n_sig = 100,
                                     fdr_cutoff = 0.05,
                                     seed = 123) {
  set.seed(seed)

  # Generate fold changes
  logFC <- rnorm(n_genes, mean = 0, sd = 2)

  # Generate p-values (mixture: some very small, most uniform)
  p_vals <- c(
    runif(n_sig, min = 1e-10, max = 0.001),      # Significant genes
    runif(n_genes - n_sig, min = 0.001, max = 1) # Non-significant
  )
  p_vals <- sample(p_vals)  # Shuffle

  # Calculate FDR (simplified BH procedure)
  adj_p_vals <- p.adjust(p_vals, method = "BH")

  # Create data frame with gene names as rownames
  df <- data.frame(
    logFC = logFC,
    P.Value = p_vals,
    adj.P.Val = adj_p_vals,
    row.names = paste0("Gene_", seq_len(n_genes))
  )

  return(df)
}

#' Generate edge case: NO significant genes
generate_no_sig_genes <- function(n_genes = 500) {
  set.seed(456)
  df <- data.frame(
    logFC = rnorm(n_genes, 0, 1),
    P.Value = runif(n_genes, 0.1, 1),
    adj.P.Val = runif(n_genes, 0.1, 1),
    row.names = paste0("Gene_", seq_len(n_genes))
  )
  return(df)
}

#' Generate edge case: ALL genes significant
generate_all_sig_genes <- function(n_genes = 100) {
  set.seed(789)
  df <- data.frame(
    logFC = rnorm(n_genes, 0, 3),
    P.Value = runif(n_genes, 1e-10, 0.001),
    adj.P.Val = runif(n_genes, 1e-10, 0.01),
    row.names = paste0("Gene_", seq_len(n_genes))
  )
  return(df)
}

#' Generate edge case: Very few significant genes (boundary case)
generate_few_sig_genes <- function(n_genes = 1000, n_sig = 3) {
  set.seed(321)

  # Most genes non-significant
  p_vals <- runif(n_genes, 0.1, 1)
  # A few with very small p-values
  p_vals[1:n_sig] <- c(1e-8, 1e-6, 1e-4)

  adj_p_vals <- p.adjust(p_vals, method = "BH")

  df <- data.frame(
    logFC = rnorm(n_genes, 0, 2),
    P.Value = p_vals,
    adj.P.Val = adj_p_vals,
    row.names = paste0("Gene_", seq_len(n_genes))
  )
  return(df)
}

#' Generate edge case: Genes exactly at FDR threshold
generate_boundary_genes <- function(n_genes = 500, fdr_cutoff = 0.05) {
  set.seed(654)

  # Create a sharp boundary at FDR = 0.05
  n_below <- 50
  n_at <- 10
  n_above <- n_genes - n_below - n_at

  p_vals <- c(
    runif(n_below, 1e-10, 1e-3),    # Clearly significant
    runif(n_at, 0.01, 0.02),        # Near boundary
    runif(n_above, 0.1, 1)          # Clearly non-significant
  )

  adj_p_vals <- p.adjust(p_vals, method = "BH")

  df <- data.frame(
    logFC = rnorm(n_genes, 0, 2),
    P.Value = p_vals,
    adj.P.Val = adj_p_vals,
    row.names = paste0("Gene_", seq_len(n_genes))
  )
  return(df)
}

## ============================================================================
## Helper Functions for Testing
## ============================================================================

#' Extract the horizontal line position from a ggplot volcano
#'
#' @param plot ggplot object
#' @return y-intercept of horizontal dashed line
get_horizontal_line <- function(plot) {
  # Extract geom_hline layer
  hline_layer <- NULL
  for (layer in plot$layers) {
    if ("GeomHline" %in% class(layer$geom)) {
      hline_layer <- layer
      break
    }
  }

  if (is.null(hline_layer)) {
    stop("No horizontal line found in plot")
  }

  return(hline_layer$data$yintercept)
}

#' Check if dashed line aligns with color boundary
#'
#' @param de_results Data frame with DE results
#' @param p_cutoff FDR cutoff
#' @param horiz_line Position of horizontal line
#' @param tolerance Allowed difference (in -log10 units)
#' @return List with test result and diagnostics
check_line_alignment <- function(de_results, p_cutoff, horiz_line, tolerance = 0.01) {
  # Find genes passing FDR threshold
  sig_genes <- de_results$adj.P.Val <= p_cutoff

  if (sum(sig_genes) == 0) {
    # No significant genes - line should be at p_cutoff
    expected_line <- -log10(p_cutoff)
  } else {
    # Line should be at max raw p-value among significant genes
    expected_line <- -log10(max(de_results$P.Value[sig_genes], na.rm = TRUE))
  }

  diff <- abs(horiz_line - expected_line)
  aligned <- diff < tolerance

  return(list(
    aligned = aligned,
    expected = expected_line,
    actual = horiz_line,
    difference = diff,
    n_sig = sum(sig_genes),
    max_sig_p = if (sum(sig_genes) > 0) max(de_results$P.Value[sig_genes]) else NA
  ))
}

## ============================================================================
## Test Suite
## ============================================================================

test_that("Standard volcano with typical data", {
  de_res <- generate_test_de_results(n_genes = 1000, n_sig = 100)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: Typical data"
  )

  expect_s3_class(plot, "ggplot")

  # Check line alignment
  hline <- get_horizontal_line(plot)
  alignment <- check_line_alignment(de_res, 0.05, hline)

  expect_true(alignment$aligned,
              info = sprintf("Line misaligned: expected %.3f, got %.3f (diff: %.5f)",
                           alignment$expected, alignment$actual, alignment$difference))
})

test_that("Standard volcano with NO significant genes", {
  de_res <- generate_no_sig_genes(n_genes = 500)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: No significant genes"
  )

  expect_s3_class(plot, "ggplot")

  hline <- get_horizontal_line(plot)
  alignment <- check_line_alignment(de_res, 0.05, hline)

  # When no significant genes, line should be at p_cutoff
  expect_equal(alignment$n_sig, 0)
  expect_true(alignment$aligned,
              info = sprintf("With 0 sig genes: expected %.3f, got %.3f",
                           alignment$expected, alignment$actual))
})

test_that("Standard volcano with ALL genes significant", {
  de_res <- generate_all_sig_genes(n_genes = 100)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: All significant"
  )

  expect_s3_class(plot, "ggplot")

  hline <- get_horizontal_line(plot)
  alignment <- check_line_alignment(de_res, 0.05, hline)

  expect_true(alignment$n_sig > 0)
  expect_true(alignment$aligned,
              info = sprintf("With all sig genes: expected %.3f, got %.3f",
                           alignment$expected, alignment$actual))
})

test_that("Standard volcano with very FEW significant genes", {
  de_res <- generate_few_sig_genes(n_genes = 1000, n_sig = 3)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: Very few significant"
  )

  expect_s3_class(plot, "ggplot")

  hline <- get_horizontal_line(plot)
  alignment <- check_line_alignment(de_res, 0.05, hline)

  expect_true(alignment$n_sig <= 5)
  expect_true(alignment$aligned,
              info = sprintf("With %d sig genes: expected %.3f, got %.3f (max_sig_p = %.2e)",
                           alignment$n_sig, alignment$expected, alignment$actual,
                           alignment$max_sig_p))
})

test_that("Standard volcano with genes at FDR boundary", {
  de_res <- generate_boundary_genes(n_genes = 500, fdr_cutoff = 0.05)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: Boundary genes"
  )

  expect_s3_class(plot, "ggplot")

  hline <- get_horizontal_line(plot)
  alignment <- check_line_alignment(de_res, 0.05, hline)

  expect_true(alignment$aligned,
              info = sprintf("Boundary case: expected %.3f, got %.3f",
                           alignment$expected, alignment$actual))
})

test_that("Standard volcano with raw p-value decision", {
  de_res <- generate_test_de_results(n_genes = 1000, n_sig = 100)

  plot <- create_standard_volcano(
    de_results = de_res,
    decision_by = "p",
    p_cutoff = 0.01,
    fc_cutoff = 1,
    title = "Test: Raw p-value"
  )

  expect_s3_class(plot, "ggplot")

  hline <- get_horizontal_line(plot)

  # With raw p-value, line should be exactly at -log10(p_cutoff)
  expected_line <- -log10(0.01)
  expect_equal(hline, expected_line, tolerance = 1e-6)
})

test_that("Vertical volcano alignment", {
  de_res <- generate_test_de_results(n_genes = 1000, n_sig = 100)

  plot <- create_vertical_volcano(
    de_results = de_res,
    decision_by = "fdr",
    p_cutoff = 0.05,
    fc_cutoff = 1,
    title = "Test: Vertical volcano"
  )

  expect_s3_class(plot, "ggplot")

  # For vertical volcano, check vertical line instead
  vline_layer <- NULL
  for (layer in plot$layers) {
    if ("GeomVline" %in% class(layer$geom)) {
      vline_layer <- layer
      break
    }
  }

  expect_false(is.null(vline_layer), info = "Vertical line should exist")

  vline <- vline_layer$data$xintercept
  alignment <- check_line_alignment(de_res, 0.05, vline)

  expect_true(alignment$aligned,
              info = sprintf("Vertical volcano: expected %.3f, got %.3f",
                           alignment$expected, alignment$actual))
})

## ============================================================================
## Visual Inspection Tests (generate plots for manual review)
## ============================================================================

cat("\n=== Running visual inspection tests ===\n")
cat("Generating plots for manual review in tests/output/...\n\n")

# Create output directory
dir.create("tests/output", showWarnings = FALSE, recursive = TRUE)

# Test case 1: Typical data
cat("1. Typical data (100 sig genes)...\n")
de_typical <- generate_test_de_results(1000, 100)
p1 <- create_standard_volcano(de_typical, title = "Typical (100/1000 sig)")
ggsave("tests/output/volcano_typical.pdf", p1, width = 8, height = 6)

# Test case 2: No significant genes
cat("2. No significant genes...\n")
de_none <- generate_no_sig_genes(500)
p2 <- create_standard_volcano(de_none, title = "No significant genes (0/500)")
ggsave("tests/output/volcano_no_sig.pdf", p2, width = 8, height = 6)

# Test case 3: Very few significant
cat("3. Very few significant (3 genes)...\n")
de_few <- generate_few_sig_genes(1000, 3)
p3 <- create_standard_volcano(de_few, title = "Very few significant (3/1000)")
ggsave("tests/output/volcano_few_sig.pdf", p3, width = 8, height = 6)

# Test case 4: All significant
cat("4. All genes significant...\n")
de_all <- generate_all_sig_genes(100)
p4 <- create_standard_volcano(de_all, title = "All significant (100/100)")
ggsave("tests/output/volcano_all_sig.pdf", p4, width = 8, height = 6)

# Test case 5: Boundary genes
cat("5. Genes at FDR boundary...\n")
de_boundary <- generate_boundary_genes(500)
p5 <- create_standard_volcano(de_boundary, title = "Boundary genes")
ggsave("tests/output/volcano_boundary.pdf", p5, width = 8, height = 6)

# Test case 6: Vertical volcano
cat("6. Vertical volcano...\n")
p6 <- create_vertical_volcano(de_typical, title = "Vertical (100/1000 sig)")
ggsave("tests/output/volcano_vertical.pdf", p6, width = 6, height = 8)

cat("\n=== Visual inspection plots saved to tests/output/ ===\n")
cat("Please review these plots to ensure dashed lines align with color boundaries.\n\n")

## ============================================================================
## Run all tests
## ============================================================================

cat("\n=== Running automated tests ===\n\n")
test_results <- test_dir(".", reporter = "summary")

cat("\n=== Test Summary ===\n")
print(test_results)
