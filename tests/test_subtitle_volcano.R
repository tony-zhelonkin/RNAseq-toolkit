#!/usr/bin/env Rscript
# Quick test of subtitle functionality in plot_standard_volcano.R

library(here)
library(limma)

# Source required helpers (relative to project root)
source(here::here("01_Scripts/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R"))

# Load checkpoints from main project
fit <- readRDS(here::here("03_Results/02_Analysis/checkpoints/fit_object.rds"))

# Get one contrast
test_tbl <- limma::topTable(fit, coef = "G32A_vs_Ctrl_D35", number = Inf, sort.by = "t")

calcium_genes <- c("NNAT","CACNG3","CACNA1C","CACNA1S","ATP2A1",
                   "RYR1","MYLK3","CASR","VDR","STIM1","STIM2",
                   "ORAI1","CALB1","CALR","PNPO")

# Test p-value mode with subtitle
message("Testing p-value mode with subtitle...")
p1 <- create_standard_volcano(
  test_tbl,
  decision_by = "p",
  p_cutoff = 0.05,
  fc_cutoff = 2,
  highlight_gene = calcium_genes,
  title = "G32A_vs_Ctrl_D35 (Threshold: p ≤ 0.05)",
  subtitle = "Highlighting by unadjusted p-value"
)

pdf(here::here("01_Scripts/RNAseq-toolkit/tests/test_p_volcano.pdf"), 8, 7)
print(p1)
dev.off()
message("✓ Saved: test_p_volcano.pdf")

# Test FDR mode with subtitle
message("Testing FDR mode with subtitle...")
p2 <- create_standard_volcano(
  test_tbl,
  decision_by = "fdr",
  p_cutoff = 0.1,
  fc_cutoff = 2,
  highlight_gene = calcium_genes,
  title = "G32A_vs_Ctrl_D35 (Threshold: FDR ≤ 0.1)",
  subtitle = "Highlighting by FDR (adj.P.Val), displayed on p-value scale for resolution"
)

pdf(here::here("01_Scripts/RNAseq-toolkit/tests/test_fdr_volcano.pdf"), 8, 7)
print(p2)
dev.off()
message("✓ Saved: test_fdr_volcano.pdf")

message("\n✅ Subtitle test complete!")
message("Check: 01_Scripts/RNAseq-toolkit/tests/test_p_volcano.pdf")
message("Check: 01_Scripts/RNAseq-toolkit/tests/test_fdr_volcano.pdf")
