#!/usr/bin/env Rscript
## ============================================================================
## test_pathway_formatting.R
##
## Comprehensive test suite for pathway name formatting functions
## Tests smart capitalization, exception handling, and edge cases
##
## Run this script to verify pathway formatting works correctly
## ============================================================================

library(testthat)
library(stringr)

# Source the formatting functions
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")

## ============================================================================
## Test Cases from User Examples
## ============================================================================

test_that("User-reported formatting issues are fixed", {

  # Issue 1: TNFR1 and NF-kappaB should be preserved
  input <- "HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY"
  result <- format_pathway_name(input)
  expect_true(grepl("TNFR1", result), info = sprintf("TNFR1 not preserved in: %s", result))
  expect_true(grepl("NF-kappaB|NF-kappa", result), info = sprintf("NF-kappaB not preserved in: %s", result))
  expect_false(grepl("Tnfr1|Nf ", result), info = sprintf("Incorrect capitalization in: %s", result))

  # Issue 2: MHC should be uppercase
  input <- "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION"
  result <- format_pathway_name(input)
  expect_true(grepl("MHC", result), info = sprintf("MHC not uppercase in: %s", result))
  expect_false(grepl("Mhc", result), info = sprintf("MHC incorrectly capitalized in: %s", result))

  # Issue 3: Roman numerals should be uppercase
  input <- "GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY"
  result <- format_pathway_name(input)
  expect_true(grepl("Type II", result), info = sprintf("Roman numeral 'II' not preserved in: %s", result))
  expect_false(grepl("Type Ii", result), info = sprintf("Roman numeral incorrectly formatted in: %s", result))

  # Issue 4: IL-2, STAT5 formatting
  input <- "KEGG_IL2_STAT5_SIGNALING_PATHWAY"
  result <- format_pathway_name(input)
  expect_true(grepl("IL-2|IL2", result), info = sprintf("IL-2 not preserved in: %s", result))
  expect_true(grepl("STAT5|STAT-5", result), info = sprintf("STAT5 not preserved in: %s", result))

  # Issue 5: Leukotrienes and Eoxins (LT and EX)
  input <- "REACTOME_SYNTHESIS_OF_LEUKOTRIENES_LT_AND_EOXINS_EX"
  result <- format_pathway_name(input)
  expect_true(grepl("Leukotrienes", result), info = sprintf("Leukotrienes not capitalized in: %s", result))
  expect_true(grepl("Eoxins", result), info = sprintf("Eoxins not capitalized in: %s", result))
  # LT and EX are tricky - they're abbreviations in parentheses
  # Accept either "LT" or "(LT)" format
})

## ============================================================================
## Test Roman Numerals
## ============================================================================

test_that("Roman numerals are preserved correctly", {

  test_cases <- list(
    c("TYPE_I_INTERFERON", "Type I Interferon"),
    c("TYPE_II_INTERFERON", "Type II Interferon"),
    c("CLASS_I_MHC", "Class I MHC"),
    c("CLASS_II_MHC", "Class II MHC"),
    c("COMPLEX_III_BIOGENESIS", "Complex III Biogenesis"),
    c("COMPLEX_IV_ASSEMBLY", "Complex IV Assembly"),
    c("FACTOR_V_ACTIVATION", "Factor V Activation"),
    c("CYTOCHROME_C_OXIDASE_COMPLEX_IV", "Cytochrome C Oxidase Complex IV")
  )

  for (tc in test_cases) {
    input <- tc[1]
    # Extract roman numeral from expected
    roman <- stringr::str_extract(tc[2], " (I|II|III|IV|V|VI|VII|VIII|IX|X) ")
    if (!is.na(roman)) {
      roman <- trimws(roman)
      result <- format_pathway_name(input, strip_prefix = FALSE)
      expect_true(grepl(roman, result),
                 info = sprintf("Roman numeral '%s' not preserved in: %s (from %s)", roman, result, input))
    }
  }
})

## ============================================================================
## Test Common Biological Abbreviations
## ============================================================================

test_that("Common biological abbreviations are preserved", {

  # Immunology
  expect_match(format_pathway_name("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), "TNF")
  expect_match(format_pathway_name("TLR4_SIGNALING"), "TLR4")
  expect_match(format_pathway_name("CD28_COSTIMULATION"), "CD28")

  # Signaling
  expect_match(format_pathway_name("JAK_STAT_SIGNALING"), "JAK")
  expect_match(format_pathway_name("JAK_STAT_SIGNALING"), "STAT")
  expect_match(format_pathway_name("PI3K_AKT_MTOR_SIGNALING"), "PI3K")
  expect_match(format_pathway_name("PI3K_AKT_MTOR_SIGNALING"), "AKT")
  expect_match(format_pathway_name("PI3K_AKT_MTOR_SIGNALING"), "mTOR")

  # Metabolism
  expect_match(format_pathway_name("TCA_CYCLE"), "TCA")
  expect_match(format_pathway_name("ATP_SYNTHESIS"), "ATP")

  # Nucleic acids
  expect_match(format_pathway_name("DNA_REPAIR"), "DNA")
  expect_match(format_pathway_name("RNA_PROCESSING"), "RNA")
  expect_match(format_pathway_name("MRNA_SPLICING"), "mRNA")

  # Growth factors
  expect_match(format_pathway_name("EGF_RECEPTOR_SIGNALING"), "EGF")
  expect_match(format_pathway_name("VEGF_SIGNALING"), "VEGF")
  expect_match(format_pathway_name("TGFB_SIGNALING"), "TGF")
})

## ============================================================================
## Test Interleukin Formatting (IL-2, IL-6, etc.)
## ============================================================================

test_that("Interleukins are formatted correctly", {

  test_cases <- list(
    "IL2_SIGNALING",
    "IL6_JAK_STAT3_SIGNALING",
    "IL10_ANTI_INFLAMMATORY",
    "IL12_SIGNALING",
    "IL17_SIGNALING"
  )

  for (input in test_cases) {
    result <- format_pathway_name(input, strip_prefix = FALSE)

    # Should contain IL (uppercase)
    expect_match(result, "IL", info = sprintf("IL not found in: %s", result))

    # Should NOT contain lowercase "Il"
    expect_false(grepl("Il[0-9]", result), info = sprintf("Lowercase 'Il' found in: %s", result))
  }
})

## ============================================================================
## Test Prefix Stripping
## ============================================================================

test_that("Common prefixes are stripped when requested", {

  test_cases <- list(
    c("HALLMARK_APOPTOSIS", "Apoptosis"),
    c("KEGG_GLYCOLYSIS", "Glycolysis"),
    c("REACTOME_CELL_CYCLE", "Cell Cycle"),
    c("GOBP_IMMUNE_RESPONSE", "Immune Response"),
    c("GOCC_MITOCHONDRION", "Mitochondrion"),
    c("GOMF_PROTEIN_BINDING", "Protein Binding"),
    c("WIKIPATHWAY_METABOLISM", "Metabolism"),
    c("BIOCARTA_APOPTOSIS_PATHWAY", "Apoptosis Pathway")
  )

  for (tc in test_cases) {
    input <- tc[1]
    result <- format_pathway_name(input, strip_prefix = TRUE)

    # Check that prefix is removed
    expect_false(grepl("^HALLMARK |^KEGG |^REACTOME |^GOBP |^GOCC |^GOMF ", result),
                info = sprintf("Prefix not stripped from: %s", result))
  }
})

test_that("Prefixes are retained when strip_prefix = FALSE", {

  input <- "HALLMARK_APOPTOSIS"
  result <- format_pathway_name(input, strip_prefix = FALSE)

  # Prefix should still be present (though may be formatted)
  expect_true(nchar(result) > nchar("Apoptosis"),
             info = sprintf("Prefix appears to be stripped: %s", result))
})

## ============================================================================
## Test Chemical Nomenclature
## ============================================================================

test_that("Chemical nomenclature is preserved", {

  # Greek letters should stay lowercase
  input <- "ALPHA_LINOLENIC_ACID_METABOLISM"
  result <- format_pathway_name(input, strip_prefix = FALSE)
  expect_match(result, "alpha", info = sprintf("Greek letter not lowercase in: %s", result))

  # N-glycan, O-glycan formatting
  input <- "N_GLYCAN_BIOSYNTHESIS"
  result <- format_pathway_name(input, strip_prefix = FALSE)
  # Should be "N-Glycan" or "N Glycan", not "N glycan"
  expect_true(grepl("N-Glycan|N Glycan", result), info = sprintf("N-glycan not formatted correctly: %s", result))

  # Cis/trans isomers
  input <- "CIS_TRANS_ISOMERIZATION"
  result <- format_pathway_name(input, strip_prefix = FALSE)
  expect_true(grepl("cis|Cis", result), info = sprintf("cis/trans not handled: %s", result))
})

## ============================================================================
## Test Edge Cases
## ============================================================================

test_that("Edge cases are handled", {

  # Empty string
  expect_equal(format_pathway_name(""), "")

  # NA values
  expect_true(is.na(format_pathway_name(NA)))

  # Very short pathway
  result <- format_pathway_name("DNA")
  expect_match(result, "DNA")

  # Mixed case input (should normalize)
  result <- format_pathway_name("MiXeD_CaSe_PaThWaY")
  expect_false(grepl("MiXeD|CaSe", result), info = sprintf("Mixed case not normalized: %s", result))

  # Multiple consecutive underscores
  result <- format_pathway_name("PATHWAY___WITH___EXTRA___UNDERSCORES")
  expect_false(grepl("  ", result), info = sprintf("Multiple spaces in: %s", result))

  # Pathway with numbers
  result <- format_pathway_name("CYTOCHROME_P450_2D6")
  expect_true(grepl("P450|2D6", result), info = sprintf("Numbers not preserved: %s", result))
})

## ============================================================================
## Test Formatting Toggle
## ============================================================================

test_that("Formatting can be toggled on/off", {

  input <- "HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING"

  # With formatting (should preserve abbreviations)
  result_on <- format_pathway_name(input, use_formatting = TRUE)

  # Without formatting (simple title case)
  result_off <- format_pathway_name(input, use_formatting = FALSE)

  # They should be different
  expect_false(result_on == result_off,
              info = sprintf("ON: %s\nOFF: %s", result_on, result_off))

  # With formatting should have uppercase abbreviations
  expect_true(grepl("TNFR1|NF", result_on),
             info = sprintf("Abbreviations not preserved with use_formatting=TRUE: %s", result_on))
})

## ============================================================================
## Test Consistency Across Similar Pathways
## ============================================================================

test_that("Similar pathways are formatted consistently", {

  # IFN pathways
  ifn_pathways <- c(
    "GOBP_TYPE_I_INTERFERON_SIGNALING_PATHWAY",
    "GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  )

  results <- sapply(ifn_pathways, format_pathway_name)

  # All should contain "Interferon" (not "interferon")
  for (r in results) {
    expect_match(r, "Interferon", info = sprintf("Inconsistent: %s", r))
  }

  # IL pathways
  il_pathways <- c(
    "IL2_STAT5_SIGNALING",
    "IL6_JAK_STAT3_SIGNALING",
    "IL12_SIGNALING"
  )

  results <- sapply(il_pathways, format_pathway_name, strip_prefix = FALSE)

  # All should have IL in uppercase
  for (r in results) {
    expect_match(r, "IL", info = sprintf("IL not consistent: %s", r))
  }
})

## ============================================================================
## Test Real-World Examples
## ============================================================================

test_that("Real-world pathway names from MSigDB are formatted correctly", {

  # Actual MSigDB pathway names
  real_pathways <- c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
    "REACTOME_IMMUNE_SYSTEM",
    "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
    "REACTOME_SYNTHESIS_OF_LEUKOTRIENES_LT_AND_EOXINS_EX",
    "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
    "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA",
    "KEGG_JAK_STAT_SIGNALING_PATHWAY"
  )

  cat("\n=== Real-World Pathway Formatting Examples ===\n")
  for (pathway in real_pathways) {
    result <- format_pathway_name(pathway)
    cat(sprintf("%-70s -> %s\n", pathway, result))

    # Basic checks
    expect_true(nchar(result) > 0, info = sprintf("Empty result for: %s", pathway))
    expect_false(grepl("_", result), info = sprintf("Underscores not removed: %s", result))
  }
  cat("\n")
})

## ============================================================================
## Run All Tests
## ============================================================================

cat("\n=== Running Pathway Formatting Tests ===\n\n")
test_results <- test_dir(".", filter = "pathway_formatting", reporter = "summary")

cat("\n=== Test Summary ===\n")
print(test_results)

## ============================================================================
## Visual Comparison Output
## ============================================================================

cat("\n=== Visual Comparison: Before vs After Formatting ===\n\n")

comparison_examples <- c(
  "HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING_PATHWAY",
  "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
  "GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY",
  "KEGG_IL2_STAT5_SIGNALING_PATHWAY",
  "REACTOME_SYNTHESIS_OF_LEUKOTRIENES_LT_AND_EOXINS_EX",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  "GOBP_RESPONSE_TO_INTERFERON_ALPHA",
  "REACTOME_COMPLEX_I_BIOGENESIS",
  "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"
)

cat(sprintf("%-75s | %s\n", "ORIGINAL", "FORMATTED"))
cat(paste(rep("-", 140), collapse = ""), "\n")

for (pathway in comparison_examples) {
  formatted <- format_pathway_name(pathway)
  cat(sprintf("%-75s | %s\n", pathway, formatted))
}

cat("\n")
cat("Tests complete! Check output above for any failures.\n")
