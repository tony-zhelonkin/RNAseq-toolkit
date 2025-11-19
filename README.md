# RNAseq toolkit 

An R toolkit for RNA-seq data analysis with personal convenience scripts, focusing on Differential Expression (DE) analysis visualization and Gene Set Enrichment Analysis (GSEA) processing and visualization. These tools are designed primarily for use with `limma`/`edgeR` DE results and `clusterProfiler`/`msigdbr` GSEA results.

## Overview

This collection provides a modular set of R functions for:

1.  **Differential Expression (DE) Visualization:** Creating volcano plots (standard and pathway-focused) and PCA plots (2D and 3D).
2.  **Gene Set Enrichment Analysis (GSEA):** Running GSEA using `clusterProfiler` with gene sets from `msigdbr`, processing results across multiple contrasts, calculating pathway scores, and generating various GSEA plots (dotplots, barplots, heatmaps, running sum plots).
3.  **Custom Visualization Theme:** A minimal ggplot2 theme (`custom_minimal_theme_with_grid`) for consistent plot aesthetics.

## Script Organization

```
scripts/
├── custom_minimal_theme.R       # Custom ggplot2 theme function
├── DE/                          # Differential Expression visualization functions
│   ├── analyzePathVolcanoViz.R  # Volcano plots highlighting specific GSEA pathways
│   ├── plot_standard_volcano.R  # Standard volcano plots for DE results
│   ├── plotPCA.R                # 2D PCA plots from DGEList objects
│   └── plotPCA3d.R              # Interactive 3D PCA plots from DGEList objects
└── GSEA/                        # GSEA analysis and visualization functions
    ├── GSEA_processing/         # Core GSEA analysis functions
    │   ├── calculate_pathway_scores.R    # Calculate average pathway scores per sample
    │   ├── get_pathway_genes_all.R       # Extract core genes for top pathways across contrasts
    │   ├── get_pathway_genes.R           # Extract core genes for significant pathways (single result)
    │   ├── get_significant_pathways.R    # Get unique significant pathway IDs from a list of results
    │   ├── run_gsea_analysis.R           # Pipeline to run GSEA & plots for multiple databases
    │   ├── run_gsea.R                    # Wrapper to run clusterProfiler::GSEA with msigdbr
    │   └── run_pooled_gsea.R             # Pipeline to run GSEA across contrasts & calculate scores
    └── GSEA_plotting/           # GSEA visualization functions
        ├── gsea_barplot.R                # Barplot of NES values for top pathways
        ├── gsea_dotplot_compare.R        # Side-by-side dotplot comparing two GSEA results
        ├── gsea_dotplot_facet.R          # Faceted dotplot separating Up/Down regulated pathways
        ├── gsea_dotplot.R                # Standard customizable GSEA dotplot
        ├── gsea_nes_comparison.R         # Scatter plot comparing NES/ES between two results
        ├── gsea_pathway_heatmap.R        # Heatmap of core enrichment gene expression for one pathway
        ├── gsea_running_sum_plot.R       # GSEA running sum enrichment plots (enrichplot::gseaplot2)
        └── gsea_scores_heatmap.R         # Heatmap of pathway scores across samples (+ wrapper)
        # Note: pathway_heatmap.R and plot_pathway_heatmap.R are likely redundant
```

## Key Functions


### GSEA Processing (`scripts/GSEA/GSEA_processing/`)

| Function                     | Description                                                                 | Key Inputs                                     | Output                                      |
| :--------------------------- | :-------------------------------------------------------------------------- | :--------------------------------------------- | :------------------------------------------ |
| `run_gsea()`                 | Runs `clusterProfiler::GSEA` using gene sets from `msigdbr`.                | DE results table, rank metric, species, MSigDB collection/subcollection | `gseaResult` object                         |
| `run_gsea_analysis()`        | Runs `run_gsea` & plotting functions for multiple databases.                | DE results table, analysis name, species, database list | List of `gseaResult` objects per database |
| `run_pooled_gsea()`          | Runs `run_gsea` across multiple contrasts, aggregates results, calculates scores. | `limma` fit, contrasts matrix, `DGEList` object | List (gsea_results, pools, genes, scores) |
| `get_significant_pathways()` | Extracts unique significant pathway IDs from a list of `gseaResult` objects. | List of `gseaResult`, `padj_cutoff`            | Character vector of pathway IDs             |
| `get_pathway_genes()`        | Extracts core genes for top significant pathways from one `gseaResult`.     | `gseaResult`, `padj_cutoff`, `top`             | Named list (PathwayID -> Gene Vector)       |
| `get_pathway_genes_all()`    | Extracts core genes for top pathways across multiple contrasts/results.     | Nested list of `gseaResult`, database name, `padj_cutoff`, `top` | Named list (PathwayID -> Gene Vector)       |
| `calculate_pathway_scores()` | Calculates average expression scores for pathways across samples.           | Expression matrix (genes x samples), pathway gene list | Matrix (samples x pathways)               |

### GSEA Visualization (`scripts/GSEA/GSEA_plotting/`)

| Function                  | Description                                                              | Key Inputs                               | Output          |
| :------------------------ | :----------------------------------------------------------------------- | :--------------------------------------- | :-------------- |
| `gsea_dotplot()`          | Creates standard customizable dotplots for GSEA results.                 | `gseaResult`, `showCategory`, `filterBy`, `sortBy`, `padj_cutoff` | ggplot object |
| `gsea_dotplot_facet()`    | Creates faceted dotplots separating Up/Down regulated pathways.          | `gseaResult`, `showCategory`, `padj_cutoff` | ggplot object |
| `gsea_dotplot_compare()`  | Creates side-by-side comparison dotplot for two `gseaResult` objects.    | `gseaResult` (x2), `pathway_ids`         | ggplot object |
| `gsea_barplot()`          | Creates barplots of NES values for top pathways.                         | `gseaResult`, `top_n`, `padj_cutoff`     | ggplot object |
| `gsea_running_sum_plot()` | Creates running sum enrichment plots (`enrichplot::gseaplot2`).          | `gseaResult`, `gene_set_ids`             | ggplot object(s) |
| `gsea_pathway_heatmap()`  | Creates heatmap of core enrichment gene expression for one pathway.      | `gseaResult`, `pathway_name`, expression matrix, `sample_order` | pheatmap object |
| `gsea_scores_heatmap()`   | Creates heatmap of pathway scores across samples.                        | Scores matrix (samples x pathways), `sample_order`, annotations | pheatmap object |
| `gsea_visualize_all_databases()` | Wrapper to generate & save heatmaps for multiple score matrices. | Output from `run_pooled_gsea`, annotations, `output_dir` | List of pheatmap objects |
| `gsea_nes_comparison()`   | Creates scatter plot comparing NES/ES between two `gseaResult` objects.  | `gseaResult` (x2), labels                | List (data, plot, common_pathways, gmt) |

### DE Visualization (`scripts/DE/`)

| Function                    | Description                                                              | Key Inputs                               | Output          |
| :-------------------------- | :----------------------------------------------------------------------- | :--------------------------------------- | :-------------- |
| `create_volcano_plot()`     | Creates standard customizable volcano plots for DE results.              | DE results table, `p_cutoff`, `fc_cutoff` | ggplot object |
| `analyze_pathway_volcano()` | Creates volcano plots highlighting genes from a specific GSEA pathway.   | `pathway_name`, `gseaResult`, DE results table | ggplot object |
| `create_pca_plot()`         | Creates 2D PCA plots from a `DGEList` object.                          | `DGEList` object                         | ggplot object |
| `create_3d_pca_plot()`      | Creates interactive 3D PCA plots (`plotly`) from a `DGEList` object.     | `DGEList` object                         | plotly object |

### Custom Theme (`scripts/`)

| Function                         | Description                                  | Key Inputs | Output          |
| :------------------------------- | :------------------------------------------- | :--------- | :-------------- |
| `custom_minimal_theme_with_grid()` | Custom `ggplot2` theme for consistent plots. | None       | ggplot theme object |

## Usage Examples

*(Note: Ensure required packages are installed and scripts are sourced correctly. Saving logic has been removed from individual plotting functions; use `ggsave()` for ggplot objects or `pheatmap::save_pheatmap_pdf()` for pheatmap objects.)*

### Basic GSEA Analysis and Dotplot

```R
# Source necessary scripts
source("scripts/GSEA/GSEA_processing/run_gsea.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("scripts/custom_minimal_theme.R") # Source the theme

# Assuming 'de_results_table' has rownames (genes) and a 't' column for ranking
# Run GSEA for Hallmark pathways (Mouse)
gsea_hallmark_res <- run_gsea(
  DE_results = de_results_table,
  rank_metric = "t",
  species = "Mus musculus",
  db_species = "MM", # Mouse abbreviation for msigdbr
  collection = "MH"  # Hallmark collection for Mouse
  # subcollection = "" # Default empty for Hallmark
)

# Create a dotplot of the results (Top 15, filtered by p.adjust < 0.05, sorted by GeneRatio)
if (!is.null(gsea_hallmark_res) && nrow(gsea_hallmark_res@result) > 0) {
  p_dotplot <- gsea_dotplot(
    gsea_hallmark_res,
    showCategory = 15,
    filterBy = "p.adjust", # Filter by significance first
    sortBy = "GeneRatio",  # Sort the top significant by GeneRatio
    padj_cutoff = 0.05,
    title = "Hallmark Pathways GSEA"
  )
  print(p_dotplot)
  # To save: ggsave("hallmark_dotplot.png", p_dotplot, width = 10, height = 7)
}
```

### Comprehensive Analysis Pipeline (Multiple Databases & Plots)

```R
# Source the pipeline script and dependencies (it sources others internally)
source("scripts/GSEA/GSEA_processing/run_gsea_analysis.R")

# Assuming 'de_results_table' is ready
# Run the comprehensive analysis for default databases (Mouse)
# This will run GSEA and generate/save multiple plots to 'results/gsea/'
all_gsea_results <- run_gsea_analysis(
  de_table = de_results_table,
  analysis_name = "Treatment_vs_Control", # Used for plot titles/filenames
  species = "Mus musculus",
  rank_metric = "t",
  padj_cutoff = 0.05, # Cutoff for plots
  output_dir = "results/gsea_analysis_output/", # Specify output directory
  save_plots = TRUE # Enable saving plots
)

# Access results for a specific database
# hallmark_results <- all_gsea_results$HALLMARK
```

### Pooled Analysis Across Contrasts

```R
# Source the pipeline script and dependencies
source("scripts/GSEA/GSEA_processing/run_pooled_gsea.R")
source("scripts/GSEA/GSEA_plotting/gsea_scores_heatmap.R") # For visualization wrapper

# Assuming 'limma_fit', 'contrast_matrix', 'dge_list', 'sample_annot', 'annot_colors' are prepared
pooled_results <- run_pooled_gsea(
  fit = limma_fit,
  contrasts = contrast_matrix,
  DGEobject = dge_list,
  species = "Mus musculus",
  top_n = 25, # Get genes/scores for top 25 pathways per DB
  padj_cutoff = 0.05, # Significance threshold for pooling
  log_file = "results/pooled_gsea_log.txt"
)

# Visualize the scores using the wrapper function (saves plots)
gsea_visualize_all_databases(
  gsea_pooled_results = pooled_results,
  annotation_col = sample_annot,
  annotation_colors = annot_colors,
  sample_order = rownames(sample_annot), # Ensure order matches annotation
  output_dir = "results/pooled_gsea_heatmaps/"
)
```

### Pathway Gene Expression Heatmap

```R
# Source necessary scripts
source("scripts/GSEA/GSEA_plotting/gsea_pathway_heatmap.R")
library(pheatmap) # Needed for saving

# Assuming 'gsea_hallmark_res', 'logcpm_matrix', 'sample_ordering', 'sample_annot' exist
# Plot heatmap for a specific pathway (e.g., by ID)
pathway_id_to_plot <- "HALLMARK_INFLAMMATORY_RESPONSE" # Example ID

heatmap_obj <- gsea_pathway_heatmap(
  gsea_obj = gsea_hallmark_res,
  pathway_name = pathway_id_to_plot,
  expression_data = logcpm_matrix,
  sample_order = sample_ordering,
  annotation_col = sample_annot,
  # annotation_colors = your_color_list, # Optional
  scale = "row" # Scale genes
)

# Save the heatmap object if generated
if (!is.null(heatmap_obj)) {
  # Use pheatmap's saving function
  pheatmap::save_pheatmap_pdf(heatmap_obj, filename = paste0("heatmap_", pathway_id_to_plot, ".pdf"))
}
```

## Recent Updates & Enhancements (2025-09-30)

### Critical Bug Fixes

1. **Environment Scoping in `run_pooled_gsea.R`**
   - **Issue**: Helper functions (`get_significant_pathways`, `get_pathway_genes_all`, `calculate_pathway_scores`) were not accessible after sourcing due to incorrect environment scoping in `sapply()` wrapper
   - **Fix**: Replaced `sapply()` sourcing with explicit `for` loop that sources into function environment using `func_env <- environment()`
   - **Impact**: Pooled GSEA now correctly finds and uses helper functions

2. **Pathway Name Formatter Sourcing**
   - **Issue**: Plotting functions tried to source `format_pathway_names.R` using `sys.frame(1)$ofile` which returned NULL when called from another script
   - **Fix**: Centralized sourcing in `run_gsea_analysis.R` by adding formatter to `helper_paths` list
   - **Impact**: All plotting functions now have access to `format_pathway_name()` for clean pathway labels

### New Features

1. **GSEA Result Caching**
   - **Purpose**: Avoid redundant GSEA computation between per-contrast and pooled analyses
   - **Implementation**:
     - Added `cache_dir` and `use_cache` parameters to `run_gsea_model1()` wrapper
     - Automatically saves GSEA results as `{contrast}_gsea_results.rds` in cache directory
     - Loads from cache on subsequent runs, dramatically reducing analysis time
     - Added `cached_gsea_results` parameter to `run_pooled_gsea()` to reuse per-contrast results
   - **Usage**:
     ```r
     # First run: computes and caches
     results <- run_gsea_model1(de_table, "KO_vs_WT_LPS",
                                cache_dir = "03_Results/GSEA/cache",
                                use_cache = TRUE)

     # Second run: loads from cache instantly
     results <- run_gsea_model1(de_table, "KO_vs_WT_LPS",
                                cache_dir = "03_Results/GSEA/cache",
                                use_cache = TRUE)

     # Pooled GSEA reuses per-contrast results
     pooled <- run_pooled_gsea(...,
                               cached_gsea_results = all_gsea_results)
     ```

2. **Enhanced Pathway Filtering Logging**
   - **Purpose**: Understand which pathways are being filtered and why
   - **Implementation**: Added `verbose` parameter to `get_significant_pathways()`
   - **Logging Output**:
     ```
     [get_significant_pathways] Processing 14 GSEA result(s)
     [get_significant_pathways] Using padj cutoff: 0.05
     [get_significant_pathways]   KO_vs_WT_LPS: 50 total pathways, 35 significant, 15 filtered out
     [get_significant_pathways] Total significant pathways before deduplication: 420
     [get_significant_pathways] Unique significant pathways: 380
     [get_significant_pathways] Removed 40 duplicates
     ```

3. **Per-Contrast Cross-Database Pooled Visualization**
   - **Purpose**: Show top pathways from ALL databases for a single contrast (user-requested feature)
   - **File**: `scripts/GSEA/GSEA_plotting/plot_pooled_contrast_heatmap.R`
   - **Function**: `plot_pooled_contrast_heatmap()`
   - **Output**: Bubble plot with:
     - X-axis: Databases (Hallmark, KEGG, GO, Reactome, etc.)
     - Y-axis: Pathways (top N per database)
     - Fill color: NES (blue = downregulated, red = upregulated)
     - Size: -log10(padj) (larger = more significant)
   - **Usage**:
     ```r
     plot_pooled_contrast_heatmap(
       gsea_results_list = all_gsea_results[["KO_vs_WT_LPS"]],
       contrast_name = "KO_vs_WT_LPS",
       top_n = 10,
       padj_cutoff = 0.05,
       output_file = "KO_vs_WT_LPS_cross_database.pdf"
     )
     ```
   - **Integration**: Automatically generated for each contrast in `02b_DEG_GSEA.R`

### Clarification: Pooled GSEA Logic

**Current "Pooled GSEA" (`run_pooled_gsea.R`):**
- Pools pathways **ACROSS contrasts** (not across databases)
- Identifies pathways significant in multiple biological comparisons
- Generates sample × pathway score heatmaps showing patterns across conditions

**New "Cross-Database Pooled" (`plot_pooled_contrast_heatmap.R`):**
- Pools pathways **ACROSS databases** for a single contrast
- Shows which databases contribute to enrichment for one comparison
- Helps identify consensus pathways across gene set collections

Both approaches are complementary and serve different analytical purposes.

## Dependencies

This toolkit relies on several R packages:

-   **Core:** `dplyr`, `stringr`, `tidyr`
-   **GSEA Analysis:** `clusterProfiler`, `msigdbr`, `stats`
-   **DE:** `limma`, `edgeR` (primarily for input object types `MArrayLM`, `DGEList` and functions like `topTable`, `cpm`)
-   **Visualization:** `ggplot2`, `ggrepel`, `enrichplot`, `pheatmap`, `plotly` (for 3D PCA), `scales`, `grDevices`
-   **Organism Databases:** Specific annotation packages like `org.Mm.eg.db` or `org.Hs.eg.db` are required depending on the species analyzed.

Ensure these packages are installed.

## Contributing

Contributions, bug reports, and suggestions are welcome. Please feel free to fork the repository and submit pull requests.

## Git Workflow & Branching Strategy

This repository uses a structured branching strategy to support multiple research projects while maintaining a stable codebase.

### Branch Structure

- **`main`**: Production-ready, stable code. Only receives tested merges from `dev`.
- **`dev`**: Main development integration branch. Receives features from project-specific branches.
- **`dev-{project}`**: Project-specific branches (e.g., `dev-GVDRP1`, `dev-Project2`). Used for project-specific development.

### Workflow for Project Development

**1. Using the toolkit in your project:**

When adding RNAseq-toolkit as a submodule to your project, configure it to track your project-specific branch:

```bash
# Add submodule and configure to track project branch
git submodule add git@github.com:tony-zhelonkin/RNAseq-toolkit.git path/to/toolkit
cd path/to/toolkit
git checkout -b dev-YourProject origin/dev  # Create from dev if new project
cd ../..

# Update .gitmodules to track your project branch
# Add this line under the submodule configuration:
#   branch = dev-YourProject
```

**2. Developing toolkit features:**

```bash
cd path/to/RNAseq-toolkit
git checkout dev-YourProject

# Make your changes, test them
git add .
git commit -m "Add feature X for Project Y"
git push origin dev-YourProject
```

**3. Contributing features back to main toolkit:**

When your project develops a feature that would benefit other projects:

```bash
# Ensure dev-YourProject is up to date
git checkout dev-YourProject
git push origin dev-YourProject

# Create Pull Request: dev-YourProject → dev
# After review and merge, all projects can benefit from your improvements
```

**4. Getting updates from the main development branch:**

Periodically sync your project branch with the latest improvements:

```bash
cd path/to/RNAseq-toolkit
git checkout dev-YourProject
git fetch origin
git merge origin/dev  # Pull in improvements from other projects
git push origin dev-YourProject

# In parent project, update the submodule reference
cd ../..
git add path/to/RNAseq-toolkit
git commit -m "Update RNAseq-toolkit submodule"
```

**5. When `dev` is ready for production:**

After testing and validation, `dev` can be merged to `main`:

```bash
# Create Pull Request: dev → main
# Requires all tests to pass
# Results in new stable release
```

### Branch Responsibilities

| Branch Type | Purpose | Update Frequency | Merge Target |
|------------|---------|------------------|--------------|
| `main` | Stable releases | Only after testing | N/A (end state) |
| `dev` | Integration of features | When features mature | `main` |
| `dev-{project}` | Project-specific work | During active development | `dev` |

### Notes

- **Not all projects need to contribute code**: Some projects may only consume the toolkit without developing new features. That's perfectly fine!
- **Project branches provide isolation**: Each project can have its own version/state of the toolkit without affecting others.
- **Contribution is optional**: Only create PRs to `dev` when you've developed something useful for other projects.
- **Stay in sync**: Periodically merge `dev` into your project branch to get improvements from other projects.

## License

MIT License

Copyright (c) 2025 Anton Zhelonkin

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
