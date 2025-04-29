# R GSEA Visualizations

A comprehensive toolkit for RNA-seq data analysis, with a focus on Gene Set Enrichment Analysis (GSEA) visualization. These tools are designed to work with differential expression results from edgeR/limma and GSEA results from the ClusterProfiler R package.

## Overview

This package provides a modular set of functions for:
1. Differential Expression (DE) analysis visualization
2. Gene Set Enrichment Analysis (GSEA) processing and visualization
3. Custom visualization themes

## Script Organization

```
scripts/
├── custom_minimal_theme.R       # Custom ggplot2 theme for consistent visualization
├── DE/                          # Differential Expression visualization
│   ├── analyzePathVolcanoViz.R  # Multi-style pathway volcano plots
│   ├── plot_standard_volcano.R  # Standard volcano plots for DE results
│   ├── plotPCA.R                # 2D PCA plots
│   └── plotPCA3d.R              # 3D PCA plots using plotly
└── GSEA/                        # GSEA analysis and visualization
    ├── GSEA_processing/         # Core GSEA analysis functions
    │   ├── calculate_pathway_scores.R    # Calculate pathway scores from expression data
    │   ├── get_pathway_genes.R           # Extract genes from significant pathways
    │   ├── get_pathway_genes_all.R       # Extract pathway genes across contrasts
    │   ├── get_significant_pathways.R    # Identify significant pathways
    │   ├── run_gsea_analysis.R           # Run comprehensive GSEA analysis with visualizations
    │   ├── runGSEA.R                     # Run GSEA on DE results
    │   └── runGSEA_pool.R                # Run pooled GSEA across multiple contrasts
    └── GSEA_plotting/           # GSEA visualization functions
        ├── gsea_barplot.R                # Create barplots of NES values
        ├── gsea_dotplot.R                # Create customizable GSEA dotplots
        ├── gsea_dotplot_compare.R        # Compare GSEA results between two conditions
        ├── gsea_dotplot_facet.R          # Create faceted dotplots with up/down regulation
        ├── gsea_nes_comparison.R         # Compare NES values between two datasets
        ├── gsea_pathway_heatmap.R        # Create heatmaps for specific pathways
        └── gsea_running_sum_plot.R       # Create running sum enrichment plots
```

## Key Functions

### GSEA Processing

| Function | Description |
|----------|-------------|
| `runGSEA()` | Runs GSEA using a ranked gene list (e.g., t-statistics) and MSigDB categories |
| `run_gsea_analysis()` | Performs comprehensive GSEA analysis with multiple visualizations |
| `run_pooled_gsea()` | Performs GSEA across multiple contrasts with robust error handling |
| `get_pathway_genes()` | Retrieves core enrichment genes from significant pathways |
| `get_significant_pathways()` | Returns pathways passing a specified FDR cutoff |
| `calculate_pathway_scores()` | Computes average expression of pathway genes per sample |

### GSEA Visualization

| Function | Description |
|----------|-------------|
| `gsea_dotplot()` | Creates customizable dotplots for GSEA results |
| `gsea_dotplot_facet()` | Creates faceted dotplots separating up/down regulated pathways |
| `gsea_dotplot_compare()` | Creates side-by-side comparison of GSEA results from two datasets |
| `gsea_barplot()` | Creates barplots of NES values for top pathways |
| `gsea_running_sum_plot()` | Creates running sum enrichment plots for specified pathways |
| `gsea_pathway_heatmap()` | Creates heatmaps showing expression of genes in a specific pathway |
| `gsea_scores_heatmap()` | Creates heatmaps of pathway scores across samples |
| `gsea_visualize_all_databases()` | Creates heatmaps for multiple GSEA database results |
| `gsea_nes_comparison()` | Creates scatter plots comparing NES values between two datasets |

### DE Visualization

| Function | Description |
|----------|-------------|
| `create_volcano_plot()` | Creates customizable volcano plots for DE results |
| `analyze_pathway_volcano()` | Creates volcano plots highlighting genes from a specific pathway |
| `create_pca_plot()` | Creates 2D PCA plots from a DGEList object |
| `create_3d_pca_plot()` | Creates interactive 3D PCA plots using plotly |

## Usage Examples

### Basic GSEA Analysis

```R
# Load the package
source("scripts/GSEA/GSEA_processing/runGSEA.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot.R")

# Run GSEA on differential expression results
gsea_result <- runGSEA(
  DE_results = de_table,
  rank_metric = "t",
  species = "Mus musculus",
  category = "H",  # Hallmark gene sets
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

# Create a dotplot of the results
gsea_dotplot(
  gsea_result,
  showCategory = 15,
  filterBy = "NES",
  sortBy = "GeneRatio",
  title = "Hallmark Pathways"
)
```

### Comprehensive GSEA Analysis

```R
# Load the package
source("scripts/GSEA/GSEA_processing/run_gsea_analysis.R")

# Run comprehensive analysis with default settings
gsea_results <- run_gsea_analysis(
  de_table = de_results,
  analysis_name = "Treatment_vs_Control",
  output_dir = "results/gsea/"
)
```

### Comparing GSEA Results Between Two Conditions

```R
# Load the package
source("scripts/GSEA/GSEA_plotting/gsea_nes_comparison.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot_compare.R")

# Compare NES values between two datasets
comparison <- gsea_nes_comparison(
  gsea_obj_x = gsea_result1,
  gsea_obj_y = gsea_result2,
  x_label = "Treatment",
  y_label = "Control"
)

# Create a comparative dotplot for common pathways
dotplot_compare <- gsea_dotplot_compare(
  gsea_obj_x = gsea_result1,
  gsea_obj_y = gsea_result2,
  pathway_ids = comparison$common_pathways$common_mix,
  sample_x_name = "Treatment",
  sample_y_name = "Control"
)
```

### Pathway Heatmaps

```R
# Load the package
source("scripts/GSEA/GSEA_plotting/gsea_pathway_heatmap.R")

# Create a heatmap for a specific pathway
gsea_pathway_heatmap(
  gsea_obj = gsea_result,
  pathway_name = "HALLMARK_APOPTOSIS",
  expression_data = norm_counts,
  sample_order = sample_order,
  annotation_col = annotation_df
)
```

## Dependencies

This package requires the following R packages:

- Core packages: dplyr, ggplot2, stringr, tidyr
- GSEA analysis: clusterProfiler, msigdbr, org.Mm.eg.db (or other organism packages)
- Visualization: ggrepel, scales, pheatmap, enrichplot, plotly (for 3D PCA)

## Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.

## License

MIT License

Copyright (c) 2024 Anton Zhelonkin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
