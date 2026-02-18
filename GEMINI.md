# GEMINI.md - RNAseq-toolkit Context

**Project:** RNAseq-toolkit (Submodule)
**Type:** R Script Library (Not a formal R package)
**Language:** R
**Version:** 2.1.0 (Inferred)

---

## Project Overview

This is a modular R toolkit designed to be used as a submodule within larger RNA-seq analysis projects (e.g., `AdaW_eWAT_WL_2025`). It provides standardized, publication-ready visualization functions for Differential Expression (DE) and Gene Set Enrichment Analysis (GSEA).

**Key Philosophy:**
- **Modular:** Source only what you need.
- **Standardized:** Enforces consistent aesthetics (themes, colors) across the parent project.
- **Wrapper-based:** Simplifies complex calls to `clusterProfiler`, `limma`, and `ggplot2`.

---

## Directory Structure

```
scripts/
├── DE/                         # Differential Expression Visualization
│   ├── plot_standard_volcano.R # (PRIMARY) Enhanced volcano plots
│   ├── plotPCA.R               # 2D PCA from DGEList
│   └── plotPCA3d.R             # Interactive 3D PCA
├── GSEA/
│   ├── GSEA_processing/        # Analysis Logic
│   │   ├── run_gsea.R          # (PRIMARY) Single-contrast GSEA wrapper
│   │   ├── run_gsea_analysis.R # Multi-database pipeline
│   │   └── run_pooled_gsea.R   # Cross-contrast aggregation
│   └── GSEA_plotting/          # Visualization
│       ├── gsea_dotplot.R      # (PRIMARY) Standard dotplot
│       ├── gsea_heatmap.R      # Pathway expression heatmaps
│       └── gsea_barplot.R      # NES barplots
└── General/                    # Utilities
    ├── annotate_genes.R        # Symbol <-> Entrez conversion
    └── custom_minimal_theme.R  # Project-wide ggplot2 theme
```

---

## Golden Path: Usage in Parent Project

Since this is a script library, functions must be `source()`'d relative to the project root.

### 1. Differential Expression Visualization

```r
# Load toolkit
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")
source("01_modules/RNAseq-toolkit/scripts/custom_minimal_theme.R")

# Generate Plot
volcano_plot <- plot_standard_volcano(
  de_results = topTable_results, # Must have: logFC, adj.P.Val, symbol
  p_cutoff = 0.05,
  fc_cutoff = 1.0, # log2 scale
  title = "Condition A vs B"
)

# Save
ggsave("03_results/plots/volcano_AvsB.pdf", volcano_plot, width = 8, height = 6)
```

### 2. GSEA Analysis & Plotting

```r
# Load toolkit
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")

# Run Analysis
gsea_res <- run_gsea(
  DE_results = de_table,    # Rownames = Gene Symbols, col = rank_metric
  rank_metric = "t",        # t-statistic recommended
  species = "Mus musculus",
  category = "H"            # Hallmark
)

# Visualize
dotplot <- gsea_dotplot(
  gsea_res,
  showCategory = 20,
  padj_cutoff = 0.05
)
```

---

## Key Development Conventions

1.  **Dependencies:**
    -   Do **not** assume libraries (`dplyr`, `ggplot2`) are loaded. Explicitly load them in scripts or ensure the parent environment has them.
    -   Key libs: `limma`, `edgeR`, `clusterProfiler`, `msigdbr`, `enrichplot`.

2.  **Sourcing:**
    -   Helper scripts often source each other. `run_gsea_analysis.R` is the "master" script that sources necessary sub-scripts.
    -   When modifying, ensure paths remain valid relative to the toolkit root.

3.  **Data Formats:**
    -   **DE Tables:** Expect columns `logFC`, `P.Value`, `adj.P.Val`, and gene symbols (either as a column `symbol` or rownames).
    -   **GSEA Objects:** Operates on `gseaResult` S4 objects from `clusterProfiler`.

4.  **Testing:**
    -   No formal `testthat` suite.
    -   Use **Visual Regression Testing**: Run `tests/test_volcano_plots.R`.
    -   Check outputs in `tests/output/` (PDFs) to verify plotting logic changes.

5.  **Git Workflow:**
    -   **`dev`**: Main integration branch.
    -   **`dev-{project}`**: Project-specific branches (e.g., `dev-AdaW`).
    -   **Submodule:** Parent projects track specific branches. Do not push breaking changes to `dev` without testing.

---

## Troubleshooting & Common Issues

-   **"Function not found":** You likely missed sourcing the specific script file. There is no `library(RNAseqToolkit)`.
-   **"Pathway name mismatch":** Ensure `species` is correct ("Mus musculus" vs "Homo sapiens") in `run_gsea`.
-   **Theme errors:** Always source `custom_minimal_theme.R` before plotting if using toolkit defaults.
