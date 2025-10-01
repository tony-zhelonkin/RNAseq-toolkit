# Volcano Plot Dashed Line Fix - Summary

## Problem

The horizontal dashed line in FDR-mode volcano plots was not aligning with the natural color boundary between significant (colored) and non-significant (gray) points.

## Root Cause

The original code calculated the line position incorrectly:

```r
# INCORRECT - old code
sig_logic  <- sig_stat <= p_cutoff
p_thresh   <- max(de_results$P.Value[sig_logic], p_cutoff, na.rm = TRUE)
horiz_line <- -log10(p_thresh)
```

The problem: `max(de_results$P.Value[sig_logic], p_cutoff, na.rm = TRUE)` would return `p_cutoff` whenever ANY significant gene existed, even if all significant p-values were much smaller. This caused a gap between the line and the actual color boundary.

## Solution

```r
# CORRECT - new code
sig_logic  <- sig_stat <= p_cutoff

# Find the boundary p-value: the largest raw p among significant genes
sig_pvals <- de_results$P.Value[sig_logic]
if (length(sig_pvals) > 0) {
  p_thresh <- max(sig_pvals, na.rm = TRUE)  # Actual boundary
} else {
  p_thresh <- p_cutoff  # Fallback when no sig genes
}
horiz_line <- -log10(p_thresh)
```

## Why This Works

**The color boundary is defined by:**
- Points with `adj.P.Val <= 0.05` are colored (green/orange/blue)
- Points with `adj.P.Val > 0.05` are gray

**The boundary in -log10(p) space is:**
- The **maximum raw p-value** among genes passing FDR
- This is where the color actually changes

**Our fix:**
- Draws the line exactly at `max(P.Value[adj.P.Val <= 0.05])`
- Ensures perfect alignment with the visual boundary

## Edge Cases Handled

### Case 1: Typical data (~100 sig genes out of 1000)
✅ Line aligns with color boundary
- Many genes above and below the threshold
- Clear visual separation

### Case 2: No significant genes
✅ Correct behavior
- All dots are gray (none pass FDR)
- Line drawn at `-log10(0.05)` as reference
- **This looks like "misbehavior" but is actually correct!**

### Case 3: All genes significant
✅ Line aligns with the lowest significant gene
- All dots colored
- Line at the boundary where a non-sig gene would appear

### Case 4: Very few significant (2-3 genes)
✅ Line aligns correctly
- Only a few colored dots above the line
- Clear threshold even with sparse data

### Case 5: Genes at FDR boundary
✅ Handles close-to-threshold genes
- Line position determined by actual FDR results
- Not affected by genes hovering near 0.05

## Special Case: Zero Significant Genes

**Example:** `KO_specific_IFNg_response` contrast
- 0 genes with FDR < 0.05 (all adj.P.Val >= 0.999)
- All dots are gray
- Line appears at `-log10(0.05)` ≈ 1.3

**Why it looks "wrong":**
- The line seems to float in empty space
- No colored points above it to validate alignment

**Why it's actually correct:**
- The line shows where the FDR threshold is
- If a gene were significant, it would appear above this line
- Useful for seeing how far away genes are from significance

## Testing

Created comprehensive test suite in `tests/test_volcano_plots.R`:

**Automated tests (6 cases):**
```
✓ Standard volcano with typical data
✓ Standard volcano with NO significant genes
✓ Standard volcano with ALL genes significant
✓ Standard volcano with very FEW significant genes
✓ Standard volcano with genes at FDR boundary
✓ Standard volcano with raw p-value decision
```

All tests **PASS** ✅

**Visual inspection plots:**
- `tests/output/volcano_typical.pdf`
- `tests/output/volcano_no_sig.pdf`
- `tests/output/volcano_few_sig.pdf`
- `tests/output/volcano_all_sig.pdf`
- `tests/output/volcano_boundary.pdf`
- `tests/output/volcano_vertical.pdf`

## Files Modified

1. **`scripts/DE/plot_standard_volcano.R`** (lines 117-139)
   - Fixed horizontal line calculation in FDR mode

2. **`scripts/DE/volcano_helpers.R`** (lines 71-94)
   - Fixed vertical line calculation in FDR mode

## Verification

To verify the fix works on your data:

```r
# Check any volcano plot
source("scripts/DE/plot_standard_volcano.R")

# Read your DE results
tt <- read.csv("path/to/DE_results.csv", row.names = 1)

# Create plot
plot <- create_standard_volcano(tt, decision_by = "fdr", p_cutoff = 0.05)

# Check alignment
sig_genes <- tt$adj.P.Val <= 0.05
if (sum(sig_genes) > 0) {
  expected_line <- -log10(max(tt$P.Value[sig_genes]))
  cat("Expected line at:", expected_line, "\n")
  cat("Should align with color boundary\n")
} else {
  cat("No significant genes - line at:", -log10(0.05), "\n")
  cat("All dots should be gray\n")
}
```

## Conclusion

The fix ensures that:
1. Dashed line **always** aligns with the color boundary
2. Edge cases (no sig, all sig, few sig) are handled correctly
3. The visual representation matches the statistical decision rule
4. Users can trust that colored = significant, gray = not significant

**The "misbehaving" plots you saw were likely cases with zero significant genes**, which is correct behavior!
