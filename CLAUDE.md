# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

RNAseq-toolkit is an R library for bulk RNA-seq differential expression (DE) analysis and Gene Set Enrichment Analysis (GSEA). It wraps `limma`/`edgeR` for DE and `clusterProfiler`/`msigdbr`/`fgsea` for GSEA, providing consistent visualization functions.

## Architecture

### Directory Structure

```
scripts/
├── General/                    # Core utilities
│   ├── annotate_genes.R        # Ensembl → Symbol/ENTREZID annotation
│   ├── dge_helpers.R           # DGEList construction (build_dge)
│   └── io_helpers.R            # File I/O utilities
├── DE/                         # Differential expression
│   ├── plot_standard_volcano.R # Main volcano plot function
│   ├── volcano_helpers.R       # Shared volcano utilities
│   ├── plotPCA.R               # 2D PCA from DGEList
│   └── plotPCA3d.R             # Interactive 3D PCA (plotly)
├── GSEA/
│   ├── GSEA_processing/        # Core GSEA analysis
│   │   ├── run_gsea.R          # Single-database GSEA (clusterProfiler wrapper)
│   │   ├── run_gsea_analysis.R # Multi-database pipeline with auto-plotting
│   │   ├── run_pooled_gsea.R   # Cross-contrast GSEA aggregation
│   │   ├── normalize_gsea.R    # Convert gseaResult → standard tibble
│   │   ├── get_pathway_genes.R # Extract leading edge genes
│   │   └── calculate_pathway_scores.R
│   └── GSEA_plotting/          # Visualization functions
│       ├── gsea_dotplot.R      # Standard dotplot
│       ├── gsea_dotplot_facet.R # Up/Down faceted dotplot
│       ├── gsea_barplot.R      # NES barplot
│       ├── gsea_running_sum_plot.R # Running sum enrichment plot
│       ├── gsea_heatmap.R      # Pathway × sample heatmaps
│       └── format_pathway_names.R # Clean MSigDB pathway names
└── custom_minimal_theme.R      # Shared ggplot2 theme
```

### Key Patterns

**1. GSEA Normalize-Then-Visualize Pattern**
```r
# Convert gseaResult objects to standardized tibbles
source("scripts/GSEA/GSEA_processing/normalize_gsea.R")
df <- normalize_gsea_results(gsea_obj, database = "Hallmark", contrast = "A_vs_B")
# Returns tibble with: pathway_id, pathway_name, NES, padj, direction, etc.
```

**2. Pipeline Functions Auto-Source Dependencies**
`run_gsea_analysis()` and `run_pooled_gsea()` source their helper scripts automatically. Pass `helper_root` if sourcing from a non-standard location.

**3. Decision-by-FDR Volcano Plots**
The volcano plot uses `-log10(P.Value)` on the y-axis but decides significance by `adj.P.Val`. The dashed line is placed at the raw p-value corresponding to the FDR boundary.

**4. GSEA Dotplot "Show All, Highlight Significant" Pattern**
`gsea_dotplot()` separates pathway **selection** (what to show) from **highlighting** (significance):
- `filterBy` + `showCategory`: Control which pathways are displayed (top N by NES, p.adjust, etc.)
- `padj_cutoff` + `highlight_threshold`: Control which pathways get black outline (FDR < threshold)

This ensures top pathways are always visible, with significant ones highlighted by a black border.

```r
gsea_dotplot(
  gsea_obj,
  filterBy = "NES",           # Sort by |NES| magnitude
  showCategory = 20,          # Show top 20 pathways
  padj_cutoff = 0.10,         # Black outline for FDR < 0.10
  highlight_sig = TRUE,       # Enable significance highlighting
  use_gradient = TRUE         # Color by NES gradient
)
```

## Running Tests

```bash
cd /path/to/RNAseq-toolkit
Rscript tests/test_volcano_plots.R
Rscript tests/test_gsea_dotplot.R
Rscript tests/test_pathway_formatting.R
```

Tests generate visual outputs in `tests/output/` for manual inspection.

## Known Issues / Compatibility

- **ggplot2 4.0+**: Use `color = "transparent"` instead of `color = NA` for shape 21 points. `NA` causes points to be removed as "missing values".

## Quick Reference: gsea_dotplot()

See `docs/API_REFERENCE.md` for full documentation.

**TL;DR — Selection vs Highlighting:**
- `filterBy` + `showCategory` → Which pathways to SHOW
- `padj_cutoff` + `highlight_threshold` → Which get BLACK OUTLINE

**Common patterns:**
```r
# Top 20 by |NES|, outline FDR < 0.10
gsea_dotplot(gsea_obj, filterBy = "NES", showCategory = 20, padj_cutoff = 0.10)

# Upregulated only
gsea_dotplot(gsea_obj, filterBy = "NES_positive", showCategory = 15)

# Downregulated only
gsea_dotplot(gsea_obj, filterBy = "NES_negative", showCategory = 15)

# Custom outline threshold
gsea_dotplot(gsea_obj, filterBy = "NES", showCategory = 30,
             padj_cutoff = 0.10, highlight_threshold = 0.01)

# No outlines at all
gsea_dotplot(gsea_obj, filterBy = "NES", showCategory = 20, highlight_sig = FALSE)
```

## Key Functions

| Function | Purpose | Input → Output |
|----------|---------|----------------|
| `run_gsea()` | Single GSEA analysis | DE table → gseaResult |
| `run_gsea_analysis()` | Multi-database GSEA + plots | DE table → list of gseaResult |
| `normalize_gsea_results()` | Standardize GSEA output | gseaResult → tibble |
| `create_standard_volcano()` | DE volcano plot | DE table → ggplot |
| `gsea_dotplot()` | GSEA dotplot | gseaResult → ggplot |
| `build_dge()` | Construct DGEList | count matrix + metadata → DGEList |

## Required R Packages

**Core:** `limma`, `edgeR`, `dplyr`, `tibble`, `ggplot2`
**GSEA:** `clusterProfiler`, `msigdbr`, `enrichplot`, `fgsea`
**Annotation:** `org.Mm.eg.db` (mouse), `org.Hs.eg.db` (human), `biomaRt`
**Visualization:** `ggrepel`, `pheatmap`, `plotly`, `scales`

## Common Usage

```r
# Source from toolkit root
toolkit_dir <- "path/to/RNAseq-toolkit"
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/run_gsea.R"))
source(file.path(toolkit_dir, "scripts/DE/plot_standard_volcano.R"))

# Run GSEA (requires DE results with gene symbols as rownames and 't' column)
gsea_res <- run_gsea(
  DE_results = de_table,
  rank_metric = "t",
  species = "Mus musculus",
  category = "H"  # Hallmark
)

# Create volcano plot
volcano <- create_standard_volcano(
  de_results = de_table,
  decision_by = "fdr",
  p_cutoff = 0.05,
  fc_cutoff = 1
)
```

## Git Branching

- `main`: Stable releases
- `dev`: Integration branch
- `dev-{project}`: Project-specific development branches

Projects using this as a submodule should track their own `dev-{project}` branch and periodically merge from `dev` to get shared improvements.
