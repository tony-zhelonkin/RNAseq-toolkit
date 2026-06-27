# RNAseq-toolkit

A modular R toolkit for bulk RNA-seq analysis, focusing on Differential Expression visualization and Gene Set Enrichment Analysis. Built for reproducibility, designed for bioinformatics workflows.

**Version:** 0.3.0
**License:** MIT
**Author:** Anton Zhelonkin

---

## Table of Contents

- [Philosophy](#philosophy)
- [Quick Start](#quick-start)
- [Architecture](#architecture)
- [Module Reference](#module-reference)
  - [General Utilities](#general-utilities)
  - [Differential Expression](#differential-expression-de)
  - [GSEA Processing](#gsea-processing)
  - [GSEA Visualization](#gsea-visualization)
  - [ORA](#ora-over-representation-analysis)
  - [Shared Utilities](#shared-utilities)
- [Design Patterns](#design-patterns)
- [Usage Examples](#usage-examples)
- [Dependencies](#dependencies)
- [Git Workflow](#git-workflow)
- [Contributing](#contributing)
- [License](#license)

---

## Philosophy

This toolkit implements the **"Normalize Once, Visualize Many"** principle:

### Core Tenets

1. **Separation of Concerns**
   - Data processing scripts never plot
   - Visualization scripts never compute
   - Config files never contain logic

2. **Checkpoint Everything Expensive**
   - Any computation taking >1 minute should be cached
   - First run: 45-60 minutes; subsequent runs: 5-10 minutes
   - Use `load_or_compute()` pattern for caching

3. **Single Source of Truth**
   - Configuration lives in one YAML file
   - Colors defined once, used everywhere
   - Schemas validated at every boundary

4. **Master Tables as Bridges**
   - CSV exports connect R computation to Python visualization
   - Standardized schemas enable cross-tool compatibility
   - No language lock-in

5. **Decision by FDR, Display by Raw P**
   - Volcano plots use FDR for significance decisions
   - Y-axis displays -log10(raw p-value) for maximal resolution
   - Avoids stair-step artifacts of plotting -log10(FDR)

### Guiding Questions

Before writing code, ask:
- Does a checkpoint already exist for this computation?
- Is this the right phase for this code (processing vs visualization)?
- Is there a single source of truth for this parameter/color/threshold?
- Can this be done with existing toolkit functions?

---

## Quick Start

```r
# Source the toolkit from your project
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")

# Run GSEA on DE results
gsea_result <- run_gsea(
  DE_results = de_table,
  rank_metric = "t",
  species = "Mus musculus",
  collection = "H"  # Hallmark; auto-detects msigdbr v7.5 vs v8+ API
)

# Create dotplot
plot <- gsea_dotplot(gsea_result, showCategory = 20, padj_cutoff = 0.05)

# Create volcano plot
volcano <- create_standard_volcano(
  de_results,
  decision_by = "fdr",
  p_cutoff = 0.05,
  fc_cutoff = 2
)

# Run ORA on a gene list (e.g. significant DE genes)
source("01_modules/RNAseq-toolkit/scripts/ORA/run_ora.R")
source("01_modules/RNAseq-toolkit/scripts/ORA/ora_dotplot.R")
sig_genes  <- rownames(de_table)[de_table$adj.P.Val < 0.05]
ora_result <- run_ora(sig_genes, species = "Mus musculus", ont = "BP")
ora_plot   <- ora_dotplot(ora_result, top_n = 20, padj_cutoff = 0.05)
```

---

## Architecture

```
RNAseq-toolkit/
├── scripts/
│   ├── General/                    # Core utilities
│   │   ├── annotate_genes.R        # Ensembl -> Symbol annotation
│   │   ├── dge_helpers.R           # DGEList construction
│   │   └── io_helpers.R            # File I/O utilities
│   │
│   ├── DE/                         # Differential Expression
│   │   ├── plot_standard_volcano.R # FDR-decision volcano plots
│   │   ├── volcano_helpers.R       # Vertical volcano, multi-panel
│   │   ├── analyzePathVolcanoViz.R # Pathway-highlighted volcano
│   │   ├── plotPCA.R               # 2D PCA visualization
│   │   ├── plotPCA3d.R             # Interactive 3D PCA
│   │   ├── create_fc_b_plot.R      # FC vs B-statistic
│   │   └── create_MD_plot.R        # Mean-Difference plot
│   │
│   ├── GSEA/
│   │   ├── GSEA_processing/        # Core analysis
│   │   │   ├── run_gsea.R              # Single-database GSEA
│   │   │   ├── run_gsea_analysis.R     # Multi-database pipeline
│   │   │   ├── run_pooled_gsea.R       # Cross-contrast aggregation
│   │   │   ├── normalize_gsea.R        # gseaResult -> tibble
│   │   │   ├── pathway_utils.R         # Gene set utilities
│   │   │   ├── parse_external_genesets.R # External DB parsers
│   │   │   ├── get_pathway_genes.R     # Leading edge extraction
│   │   │   ├── get_pathway_genes_all.R # Cross-contrast genes
│   │   │   ├── get_significant_pathways.R # Pathway pooling
│   │   │   └── calculate_pathway_scores.R # Sample scores
│   │   │
│   │   ├── GSEA_plotting/          # Visualization
│   │   │   ├── gsea_dotplot.R          # Main dotplot
│   │   │   ├── gsea_dotplot_facet.R    # Up/Down faceted
│   │   │   ├── gsea_dotplot_compare.R  # Side-by-side comparison
│   │   │   ├── gsea_barplot.R          # NES barplot
│   │   │   ├── gsea_running_sum_plot.R # Enrichment curve
│   │   │   ├── gsea_heatmap.R          # Pathway heatmaps
│   │   │   ├── format_pathway_names.R  # Smart capitalization
│   │   │   ├── gsea_plotting_utils.R   # Shared helpers
│   │   │   └── ...                     # Additional plots
│   │   │
│   │   └── GSEA_plotting_python/   # Python alternatives
│   │
│   ├── ORA/                        # Over-Representation Analysis
│   │   ├── run_ora.R               # Fisher-exact ORA
│   │   └── ora_dotplot.R           # ORA dotplot visualization
│   │
│   ├── custom_minimal_theme.R      # Publication ggplot2 theme
│   └── utils_plotting.R            # DRY utilities
│
├── tests/                          # Visual regression tests
├── examples/                       # Usage examples
└── docs/                           # Extended documentation
```

---

## Module Reference

### General Utilities

| Function | File | Description |
|----------|------|-------------|
| `annotate_genes_from_ensembl()` | `annotate_genes.R` | Convert Ensembl IDs to symbols using org.*.eg.db |
| `build_dge()` | `dge_helpers.R` | Construct validated DGEList with TMM normalization |
| `read_counts_matrix()` | `io_helpers.R` | Flexible featureCounts/generic count parser |
| `read_metadata()` | `io_helpers.R` | Excel metadata reader with validation |
| `align_metadata_to_counts()` | `io_helpers.R` | Sync metadata to count matrix columns |

### Differential Expression (DE)

| Function | File | Description |
|----------|------|-------------|
| `create_standard_volcano()` | `plot_standard_volcano.R` | FDR-decision volcano with raw-p y-axis |
| `create_vertical_volcano()` | `volcano_helpers.R` | 90-degree rotated volcano for grids |
| `combine_volcano_row()` | `volcano_helpers.R` | Multi-panel volcano with unified legend |
| `analyze_pathway_volcano()` | `analyzePathVolcanoViz.R` | Highlight genes from GSEA pathway |
| `create_pca_plot()` | `plotPCA.R` | Standard 2D PCA from DGEList |
| `create_3d_pca_plot()` | `plotPCA3d.R` | Interactive 3D PCA (plotly) |
| `create_fc_b_plot()` | `create_fc_b_plot.R` | logFC vs B-statistic scatter |
| `create_MD_plot()` | `create_MD_plot.R` | Mean-Difference (MA) plot |

### GSEA Processing

| Function | File | Description |
|----------|------|-------------|
| `run_gsea()` | `run_gsea.R` | Single-database GSEA; auto-detects msigdbr v7.5/v8+ API |
| `run_gsea_analysis()` | `run_gsea_analysis.R` | Multi-database pipeline with auto-plotting |
| `run_pooled_gsea()` | `run_pooled_gsea.R` | Cross-contrast aggregation with scoring |
| `normalize_gsea_results()` | `normalize_gsea.R` | Convert gseaResult to standardized tibble |
| `filter_pathways_by_size()` | `pathway_utils.R` | Filter gene sets by size |
| `list_to_term2gene()` | `pathway_utils.R` | Convert list to TERM2GENE format |
| `parse_transportdb()` | `parse_external_genesets.R` | Parse TransportDB2.0 |
| `parse_gmt()` | `parse_external_genesets.R` | Parse GMT gene set files |
| `get_pathway_genes()` | `get_pathway_genes.R` | Extract leading edge genes |
| `get_significant_pathways()` | `get_significant_pathways.R` | Pool significant pathway IDs |
| `calculate_pathway_scores()` | `calculate_pathway_scores.R` | Sample-wise pathway scores |

### GSEA Visualization

| Function | File | Description |
|----------|------|-------------|
| `gsea_dotplot()` | `gsea_dotplot.R` | Customizable GSEA dotplot |
| `gsea_dotplot_facet()` | `gsea_dotplot_facet.R` | Separate Up/Down panels |
| `gsea_dotplot_compare()` | `gsea_dotplot_compare.R` | Side-by-side comparison |
| `gsea_barplot()` | `gsea_barplot.R` | NES horizontal barplot |
| `gsea_running_sum_plot()` | `gsea_running_sum_plot.R` | Enrichment curve; works with non-MSigDB external DBs |
| `format_pathway_name()` | `format_pathway_names.R` | Smart biological capitalization |
| `gsea_heatmap_save()` | `gsea_heatmap.R` | Sample x pathway heatmap |
| `plot_pooled_contrast_dotplot()` | `plot_pooled_contrast_dotplot.R` | Cross-database dotplot |

### ORA (Over-Representation Analysis)

| Function | File | Description |
|----------|------|-------------|
| `run_ora()` | `ORA/run_ora.R` | Fisher-exact ORA against MSigDB or custom gene sets |
| `ora_dotplot()` | `ORA/ora_dotplot.R` | Dotplot for ORA results with significance highlighting |

### Shared Utilities

| Function | File | Description |
|----------|------|-------------|
| `custom_minimal_theme_with_grid()` | `custom_minimal_theme.R` | Publication-ready ggplot2 theme |
| `ensure_dir()` | `utils_plotting.R` | Create directories idempotently |
| `save_plot()` | `utils_plotting.R` | PDF saving with device management |
| `load_checkpoint()` | `utils_plotting.R` | RDS checkpoint loader |
| `log_message()` | `utils_plotting.R` | Timestamped logging |

---

## Design Patterns

### 1. Volcano Plot: FDR Decision, Raw-P Axis

The `create_standard_volcano()` function implements a best-practice approach:

```r
# Why raw p on the y-axis?
# - Raw p-values are per-gene test statistics
# - FDR is an average over the rejected set
# - Plotting raw p gives maximal resolution
# - Avoids stair-step artifacts of -log10(FDR)

volcano <- create_standard_volcano(
  de_results,
  decision_by = "fdr",     # Use FDR for color decisions
  p_cutoff = 0.05,         # FDR threshold
  fc_cutoff = 2            # log2FC threshold
)
# Horizontal dashed line: boundary p-value where FDR = cutoff
# This aligns EXACTLY with the color transition
```

### 2. Smart Biological Term Formatting

The `format_pathway_name()` function handles biological nomenclature:

```r
format_pathway_name("HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING")
# Returns: "TNFR1-Induced NF-kappaB Signaling"

format_pathway_name("GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY")
# Returns: "Type II Interferon Signaling Pathway"
```

**Features:**
- 400+ exception dictionary for biological terms
- Multi-word pattern recognition ("nf kappa b" -> "NF-kappaB")
- Greek letter preservation
- Roman numeral handling
- Chemical prefix awareness

### 3. GSEA Result Normalization

Convert clusterProfiler results to standardized tibbles:

```r
source("scripts/GSEA/GSEA_processing/normalize_gsea.R")

# Convert gseaResult to standardized tibble
normalized <- normalize_gsea_results(
  gsea_result,
  database = "Hallmark",
  contrast = "Treatment_vs_Control"
)

# Schema:
# pathway_id, pathway_name, database, contrast, NES, pvalue, padj,
# set_size, core_enrichment, direction
```

### 4. External Gene Set Integration

Support for custom databases beyond MSigDB:

```r
source("scripts/GSEA/GSEA_processing/parse_external_genesets.R")

# Parse TransportDB
transport_sets <- parse_transportdb("TransportDB2.0.csv")

# Parse GMT file
custom_sets <- parse_gmt("custom_pathways.gmt")

# Use with clusterProfiler
gsea_result <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = transport_sets$TERM2GENE,
  TERM2NAME = transport_sets$TERM2NAME
)
```

---

## Usage Examples

### Basic GSEA Analysis

```r
# Source required scripts
source("scripts/GSEA/GSEA_processing/run_gsea.R")
source("scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")

# DE results from limma (rownames = gene symbols)
de_table <- topTable(fit, coef = "Treatment_vs_Control", n = Inf)

# Run GSEA for multiple databases
databases <- list(
  H = c("H", ""),
  KEGG = c("C2", "CP:KEGG"),
  GO_BP = c("C5", "GO:BP")
)

results <- list()
for (db_name in names(databases)) {
  results[[db_name]] <- run_gsea(
    DE_results = de_table,
    rank_metric = "t",
    species = "Mus musculus",
    collection    = databases[[db_name]][1],
    subcollection = databases[[db_name]][2]
    # run_gsea() auto-detects msigdbr v7.5 (category/subcategory) vs v8+ (collection/subcollection)
  )
}

# Create dotplots
for (db_name in names(results)) {
  plot <- gsea_dotplot(
    results[[db_name]],
    showCategory = 20,
    padj_cutoff = 0.05,
    title = paste(db_name, "Pathways")
  )
  ggsave(paste0("dotplot_", db_name, ".pdf"), plot, width = 10, height = 8)
}
```

### Comprehensive Pipeline

```r
# Use run_gsea_analysis for automated multi-database analysis
source("scripts/GSEA/GSEA_processing/run_gsea_analysis.R")

all_results <- run_gsea_analysis(
  de_table = de_table,
  analysis_name = "Treatment_vs_Control",
  species = "Mus musculus",
  rank_metric = "t",
  padj_cutoff = 0.05,
  output_dir = "results/GSEA/"
)

# Automatically generates:
# - Dotplots (standard, faceted up/down)
# - Barplots
# - Running sum plots for top pathways
# - Heatmaps
```

### Pooled Analysis Across Contrasts

```r
source("scripts/GSEA/GSEA_processing/run_pooled_gsea.R")

# Run GSEA across all contrasts and aggregate
pooled <- run_pooled_gsea(
  fit = limma_fit,
  contrasts = contrast_matrix,
  DGEobject = dge_list,
  species = "Mus musculus",
  top_n = 25,
  padj_cutoff = 0.05
)

# Returns:
# - gsea_results: per-contrast results
# - pools: significant pathway IDs per database
# - genes: leading edge genes
# - scores: sample x pathway score matrix
```

### Custom Volcano Plot

```r
source("scripts/DE/plot_standard_volcano.R")

volcano <- create_standard_volcano(
  de_results,
  decision_by = "fdr",
  p_cutoff = 0.05,
  fc_cutoff = 1,
  top_n = 10,
  highlight_gene = c("Il6", "Tnf", "Ccl2"),  # Priority labels
  title = "Treatment vs Control",
  subtitle = "FDR < 0.05, |log2FC| > 1"
)

ggsave("volcano.pdf", volcano, width = 10, height = 8)
```

---

## Dependencies

### Core Packages

```r
# Data manipulation
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# Differential expression
library(limma)
library(edgeR)

# GSEA
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(enrichplot)

# Annotation
library(org.Mm.eg.db)  # Mouse
library(org.Hs.eg.db)  # Human
library(AnnotationDbi)

# Visualization
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(plotly)  # 3D PCA
library(scales)
```

### Version Requirements

- R >= 4.0
- msigdbr >= 7.5 (`run_gsea()` auto-detects v7.5 `category`/`subcategory` vs v8+ `collection`/`subcollection` API)
- clusterProfiler >= 4.0

---

## Git Workflow

### Branch Structure

| Branch | Purpose | Merge Target |
|--------|---------|--------------|
| `main` | Stable releases | N/A |
| `dev` | Integration branch | `main` |
| `dev-{project}` | Project-specific | `dev` |

### For Project Use

```bash
# Add as submodule tracking your project branch
git submodule add git@github.com:user/RNAseq-toolkit.git
cd RNAseq-toolkit
git checkout -b dev-YourProject origin/dev

# Update .gitmodules
[submodule "RNAseq-toolkit"]
    path = path/to/toolkit
    url = git@github.com:user/RNAseq-toolkit.git
    branch = dev-YourProject
```

### Contributing Features Back

```bash
# Push to your project branch
git push origin dev-YourProject

# Create PR: dev-YourProject -> dev
# After review/merge, all projects benefit
```

---

## Contributing

1. **Fork the repository**
2. **Create a feature branch** from `dev`
3. **Follow coding standards:**
   - snake_case for functions/variables
   - 2-space indentation
   - roxygen2 documentation
4. **Add visual regression tests** for plotting functions
5. **Update documentation** if adding new features
6. **Submit PR** to `dev` branch

Follow snake_case naming, 2-space indentation, and add roxygen2 docs + a visual regression test for any new plotting function.

---

## License

MIT License - Copyright (c) 2025 Anton Zhelonkin

See [LICENSE.md](LICENSE.md) for full text.
