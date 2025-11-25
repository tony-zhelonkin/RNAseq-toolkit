#!/usr/bin/env Rscript
###############################################################################
## Integrated test of all volcano plot improvements
## Tests: subtitle, calcium gene layer separation, updated gene list
###############################################################################

library(here)
library(limma)

# Source required helpers
source(here::here("01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R"))

# Load checkpoints
fit <- readRDS(here::here("03_Results/02_Analysis/checkpoints/fit_object.rds"))

# Updated calcium genes list (12 genes - corrected)
calcium_genes <- c("NNAT","CACNG3","CACNA1S","ATP2A1",
                   "RYR1","MYLK3","VDR","STIM1","STIM2",
                   "ORAI1_1","CALB1","CALR","PNPO")

message("Testing integrated volcano plot improvements...")
message(sprintf("  Using %d calcium genes (corrected list)\n", length(calcium_genes)))

# Test 1: Disease vs Control (should have significant genes)
message("Test 1: G32A_vs_Ctrl_D35 (disease vs control)")
test1 <- limma::topTable(fit, coef = "G32A_vs_Ctrl_D35", number = Inf, sort.by = "t")

p1_fdr <- create_standard_volcano(
  test1,
  decision_by = "fdr",
  p_cutoff = 0.1,
  fc_cutoff = 2,
  label_method = "top",
  highlight_gene = calcium_genes,
  title = "G32A_vs_Ctrl_D35 (Threshold: FDR ≤ 0.1)",
  subtitle = "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
)

pdf(here::here("01_Scripts/RNAseq-toolkit/tests/test1_G32A_fdr.pdf"), 10, 8)
print(p1_fdr)
dev.off()
message("  ✓ Saved: test1_G32A_fdr.pdf")

# Test 2: Maturation contrast (few/no significant genes - test label_method='top')
message("\nTest 2: Maturation_G32A_specific (maturation effect)")
test2 <- limma::topTable(fit, coef = "Maturation_G32A_specific", number = Inf, sort.by = "t")

p2_fdr <- create_standard_volcano(
  test2,
  decision_by = "fdr",
  p_cutoff = 0.1,
  fc_cutoff = 2,
  label_method = "top",  # Should show top genes even if not significant
  highlight_gene = calcium_genes,
  title = "Maturation_G32A_specific (Threshold: FDR ≤ 0.1)",
  subtitle = "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
)

pdf(here::here("01_Scripts/RNAseq-toolkit/tests/test2_Maturation_G32A_fdr.pdf"), 10, 8)
print(p2_fdr)
dev.off()
message("  ✓ Saved: test2_Maturation_G32A_fdr.pdf")

# Test 3: P-value mode for comparison
message("\nTest 3: G32A_vs_Ctrl_D35 (p-value mode)")
p3_p <- create_standard_volcano(
  test1,
  decision_by = "p",
  p_cutoff = 0.05,
  fc_cutoff = 2,
  label_method = "top",
  highlight_gene = calcium_genes,
  title = "G32A_vs_Ctrl_D35 (Threshold: p ≤ 0.05)",
  subtitle = "Highlighting by unadjusted p-value"
)

pdf(here::here("01_Scripts/RNAseq-toolkit/tests/test3_G32A_p.pdf"), 10, 8)
print(p3_p)
dev.off()
message("  ✓ Saved: test3_G32A_p.pdf")

# Diagnostic: Count calcium genes present
message("\n=== Calcium Gene Presence Diagnostic ===")
present_test1 <- sum(calcium_genes %in% rownames(test1))
present_test2 <- sum(calcium_genes %in% rownames(test2))

message(sprintf("  G32A_vs_Ctrl_D35: %d/%d calcium genes present",
                present_test1, length(calcium_genes)))
message(sprintf("  Maturation_G32A_specific: %d/%d calcium genes present",
                present_test2, length(calcium_genes)))

message("\n✅ Integration test complete!")
message("\nExpected results:")
message("  • All plots should show subtitles distinguishing p vs fdr approach")
message("  • Calcium genes should appear as BOLD BLACK labels")
message("  • Maturation plot should show top genes (not empty)")
message("  • All 12 corrected calcium genes should be visible where present")
message("\nOutput files:")
message("  01_Scripts/RNAseq-toolkit/tests/test1_G32A_fdr.pdf")
message("  01_Scripts/RNAseq-toolkit/tests/test2_Maturation_G32A_fdr.pdf")
message("  01_Scripts/RNAseq-toolkit/tests/test3_G32A_p.pdf")
