# Volcano Plot Test Suite

This directory contains comprehensive tests for the volcano plot functions to ensure correct behavior across various edge cases.

## Purpose

The main issue being tested is **dashed line alignment with color boundaries** in FDR mode. The horizontal/vertical dashed line should align exactly with the boundary between colored (significant) and gray (non-significant) points.

## Test Cases

### Automated Tests

The test suite includes these edge cases:

1. **Typical data**: ~100 significant genes out of 1000
2. **No significant genes**: All genes fail FDR threshold
3. **All genes significant**: All genes pass FDR threshold
4. **Very few significant**: Only 2-3 genes significant
5. **Boundary genes**: Genes with adj.P.Val very close to cutoff
6. **Raw p-value mode**: Using `decision_by = "p"` instead of FDR
7. **Vertical volcano**: Same tests for vertical orientation

### Visual Inspection Tests

The script also generates plots for manual review in `tests/output/`:

- `volcano_typical.pdf` - Standard case
- `volcano_no_sig.pdf` - No significant genes
- `volcano_few_sig.pdf` - Very few significant genes
- `volcano_all_sig.pdf` - All genes significant
- `volcano_boundary.pdf` - Genes near FDR boundary
- `volcano_vertical.pdf` - Vertical orientation

## Running the Tests

From the GSEA submodule root directory:

```bash
cd /workspaces/13036-DM_DMlab_summer_2025/01_Scripts/GSEA
Rscript tests/test_volcano_plots.R
```

This will:
1. Run automated tests with `testthat`
2. Generate PDF plots for visual inspection
3. Print a summary of test results

## Expected Behavior

**Correct alignment** means:

- In FDR mode with significant genes: dashed line at `-log10(max(P.Value[adj.P.Val <= cutoff]))`
- In FDR mode with NO significant genes: dashed line at `-log10(cutoff)`
- In raw p mode: dashed line at exactly `-log10(cutoff)`

The line should **perfectly align** with the boundary where point colors change from significant (colored) to non-significant (gray).

## Key Fix Applied

The original code had:
```r
p_thresh <- max(de_results$P.Value[sig_logic], p_cutoff, na.rm = TRUE)
```

This incorrectly returned `p_cutoff` when ANY significant gene existed, causing misalignment.

The fix:
```r
sig_pvals <- de_results$P.Value[sig_logic]
if (length(sig_pvals) > 0) {
  p_thresh <- max(sig_pvals, na.rm = TRUE)  # Actual boundary
} else {
  p_thresh <- p_cutoff  # Fallback when no sig genes
}
```

## Dependencies

- `testthat` for automated testing
- `ggplot2` for plots
- `ggrepel` for labels
- The volcano plot functions from `scripts/DE/`

## Interpreting Results

All tests should **PASS**. If any test fails, it indicates:

1. **Line misalignment**: The dashed line doesn't match the color boundary
2. **Edge case failure**: An uncommon scenario (no sig genes, all sig, etc.) breaks the function
3. **Logic error**: The significance decision and line position use different thresholds

Visual inspection plots should show:
- Dashed line exactly at the transition from colored to gray dots
- No gap between the line and the color boundary
- Consistent behavior across all edge cases
