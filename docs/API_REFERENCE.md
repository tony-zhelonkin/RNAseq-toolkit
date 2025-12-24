# API Reference

Complete function documentation for RNAseq-toolkit.

---

## Table of Contents

1. [GSEA Processing](#gsea-processing)
2. [GSEA Visualization](#gsea-visualization)
3. [Differential Expression](#differential-expression)
4. [General Utilities](#general-utilities)
5. [Shared Utilities](#shared-utilities)

---

## GSEA Processing

### run_gsea()

**File:** `scripts/GSEA/GSEA_processing/run_gsea.R`

Run Gene Set Enrichment Analysis using clusterProfiler with MSigDB gene sets.

```r
run_gsea(
  DE_results,
  rank_metric = "t",
  species = "Mus musculus",
  db_species = NULL,
  category = "H",
  subcategory = NULL,
  pvalue_cutoff = 1,
  padj_method = "fdr",
  nperm = 100000,
  seed = 123
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `DE_results` | data.frame | required | DE results with gene IDs as rownames |
| `rank_metric` | character | "t" | Column name for ranking genes |
| `species` | character | "Mus musculus" | Species for msigdbr |
| `db_species` | character | NULL | Species code (MM, HS); auto-detected if NULL |
| `category` | character | "H" | MSigDB category (H, C2, C5, etc.) |
| `subcategory` | character | NULL | MSigDB subcategory (CP:KEGG, GO:BP, etc.) |
| `pvalue_cutoff` | numeric | 1 | P-value cutoff for GSEA |
| `padj_method` | character | "fdr" | P-value adjustment method |
| `nperm` | integer | 100000 | Number of permutations |
| `seed` | integer | 123 | Random seed for reproducibility |

**Returns:** `gseaResult` object from clusterProfiler

**Example:**

```r
# Hallmark pathways
gsea_h <- run_gsea(de_results, category = "H")

# KEGG pathways
gsea_kegg <- run_gsea(de_results, category = "C2", subcategory = "CP:KEGG")

# GO Biological Process
gsea_go <- run_gsea(de_results, category = "C5", subcategory = "GO:BP")
```

---

### run_gsea_analysis()

**File:** `scripts/GSEA/GSEA_processing/run_gsea_analysis.R`

Run GSEA for multiple databases with automatic plot generation.

```r
run_gsea_analysis(
  de_table,
  analysis_name,
  species = "Mus musculus",
  rank_metric = "t",
  databases = NULL,
  padj_cutoff = 0.05,
  output_dir = "results/GSEA/",
  save_plots = TRUE,
  helper_root = NULL
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `de_table` | data.frame | required | DE results table |
| `analysis_name` | character | required | Name for output files/titles |
| `species` | character | "Mus musculus" | Species for msigdbr |
| `rank_metric` | character | "t" | Column for gene ranking |
| `databases` | list | NULL | MSigDB databases to analyze (default: H, KEGG, GO, Reactome) |
| `padj_cutoff` | numeric | 0.05 | Significance cutoff for plots |
| `output_dir` | character | "results/GSEA/" | Output directory |
| `save_plots` | logical | TRUE | Whether to save plots |
| `helper_root` | character | NULL | Root path for helper scripts |

**Returns:** Named list of `gseaResult` objects per database

**Default Databases:**

```r
list(
  H = c("H", ""),
  C2_KEGG = c("C2", "CP:KEGG"),
  C2_REACTOME = c("C2", "CP:REACTOME"),
  C2_WIKIPATHWAYS = c("C2", "CP:WIKIPATHWAYS"),
  C3_TF = c("C3", "TFT:GTRD"),
  C5_BP = c("C5", "GO:BP"),
  C5_MF = c("C5", "GO:MF"),
  C5_CC = c("C5", "GO:CC")
)
```

---

### run_pooled_gsea()

**File:** `scripts/GSEA/GSEA_processing/run_pooled_gsea.R`

Run GSEA across multiple contrasts and aggregate results.

```r
run_pooled_gsea(
  fit,
  contrasts,
  DGEobject,
  species = "Mus musculus",
  databases = NULL,
  top_n = 25,
  padj_cutoff = 0.05,
  cached_gsea_results = NULL,
  log_file = NULL
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fit` | MArrayLM | required | limma fit object |
| `contrasts` | matrix | required | Contrast matrix |
| `DGEobject` | DGEList | required | DGEList with expression data |
| `species` | character | "Mus musculus" | Species for msigdbr |
| `databases` | list | NULL | MSigDB databases |
| `top_n` | integer | 25 | Top pathways per database |
| `padj_cutoff` | numeric | 0.05 | Significance threshold |
| `cached_gsea_results` | list | NULL | Pre-computed GSEA results |
| `log_file` | character | NULL | Path for log output |

**Returns:** List with:
- `gsea_results`: Per-contrast GSEA results
- `pools`: Significant pathway IDs per database
- `genes`: Leading edge genes per pathway
- `scores`: Sample x pathway score matrix

---

### normalize_gsea_results()

**File:** `scripts/GSEA/GSEA_processing/normalize_gsea.R`

Convert gseaResult to standardized tibble format.

```r
normalize_gsea_results(
  gsea_result,
  database = "Unknown",
  contrast = "Unknown",
  format_names = TRUE
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_result` | gseaResult | required | clusterProfiler result |
| `database` | character | "Unknown" | Database name |
| `contrast` | character | "Unknown" | Contrast name |
| `format_names` | logical | TRUE | Apply format_pathway_name() |

**Returns:** tibble with standardized schema:

```r
tibble(
  pathway_id = character(),
  pathway_name = character(),
  database = character(),
  contrast = character(),
  NES = numeric(),
  pvalue = numeric(),
  padj = numeric(),
  set_size = integer(),
  core_enrichment = character(),
  direction = character()
)
```

---

### parse_transportdb()

**File:** `scripts/GSEA/GSEA_processing/parse_external_genesets.R`

Parse TransportDB2.0 CSV file into GSEA-compatible format.

```r
parse_transportdb(
  file_path,
  species = "Mus musculus"
)
```

**Returns:** List with `TERM2GENE` and `TERM2NAME` data frames

---

### parse_gmt()

**File:** `scripts/GSEA/GSEA_processing/parse_external_genesets.R`

Parse GMT gene set file.

```r
parse_gmt(file_path)
```

**Returns:** List with `TERM2GENE` and `TERM2NAME` data frames

---

### get_pathway_genes()

**File:** `scripts/GSEA/GSEA_processing/get_pathway_genes.R`

Extract leading edge genes for significant pathways.

```r
get_pathway_genes(
  gsea_result,
  padj_cutoff = 0.05,
  top = 20
)
```

**Returns:** Named list: pathway_id -> character vector of genes

---

### calculate_pathway_scores()

**File:** `scripts/GSEA/GSEA_processing/calculate_pathway_scores.R`

Calculate sample-wise pathway activity scores.

```r
calculate_pathway_scores(
  expression_matrix,
  pathway_genes,
  method = "mean"
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expression_matrix` | matrix | required | Genes x samples expression |
| `pathway_genes` | list | required | Named list of gene vectors |
| `method` | character | "mean" | Aggregation method |

**Returns:** Matrix (samples x pathways)

---

## GSEA Visualization

### gsea_dotplot()

**File:** `scripts/GSEA/GSEA_plotting/gsea_dotplot.R`

Create customizable GSEA dotplot.

```r
gsea_dotplot(
  gsea_result,
  n_top = 20,
  filterBy = "p.adjust",
  sortBy = "NES",
  padj_cutoff = 0.05,
  NES_cutoff = NULL,
  show_only = NULL,
  title = NULL,
  color_scale = "diverging",
  highlight_sig = TRUE
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_result` | gseaResult | required | GSEA result object |
| `n_top` | integer | 20 | Number of pathways to show |
| `filterBy` | character | "p.adjust" | Column for filtering |
| `sortBy` | character | "NES" | Column for sorting |
| `padj_cutoff` | numeric | 0.05 | Significance threshold |
| `NES_cutoff` | numeric | NULL | Minimum absolute NES |
| `show_only` | character | NULL | "NES_positive" or "NES_negative" |
| `title` | character | NULL | Plot title |
| `color_scale` | character | "diverging" | Color scheme |
| `highlight_sig` | logical | TRUE | Highlight significant points |

**Returns:** ggplot object

---

### gsea_dotplot_facet()

**File:** `scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R`

Create faceted dotplot separating up/down regulated pathways.

```r
gsea_dotplot_facet(
  gsea_result,
  n_top = 10,
  padj_cutoff = 0.05,
  title = NULL
)
```

**Returns:** ggplot object with Up/Down facets

---

### gsea_barplot()

**File:** `scripts/GSEA/GSEA_plotting/gsea_barplot.R`

Create horizontal barplot of NES values.

```r
gsea_barplot(
  gsea_result,
  top_n = 15,
  padj_cutoff = 0.05,
  title = NULL
)
```

**Returns:** ggplot object

---

### gsea_running_sum_plot()

**File:** `scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R`

Create classic GSEA enrichment curve.

```r
gsea_running_sum_plot(
  gsea_result,
  gene_set_id,
  title = NULL,
  description = NULL
)
```

**Note:** Wrapper for `enrichplot::gseaplot2()`. If `description` is provided (e.g., a formatted pathway name), it is used as the plot title.

**Returns:** ggplot object

---

### format_pathway_name()

**File:** `scripts/GSEA/GSEA_plotting/format_pathway_names.R`

Format pathway names with smart biological capitalization.

```r
format_pathway_name(
  text,
  use_formatting = TRUE,
  strip_prefix = TRUE
)
```

**Arguments:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `text` | character | required | Pathway name(s) |
| `use_formatting` | logical | TRUE | Apply smart formatting |
| `strip_prefix` | logical | TRUE | Remove database prefixes |

**Examples:**

```r
format_pathway_name("HALLMARK_TNFR1_INDUCED_NF_KAPPA_B_SIGNALING")
# "TNFR1-Induced NF-kappaB Signaling"

format_pathway_name("GOBP_TYPE_II_INTERFERON_SIGNALING_PATHWAY")
# "Type II Interferon Signaling Pathway"

format_pathway_name("REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING")
# "Class I MHC-Mediated Antigen Processing"
```

**Exception Dictionary Includes:**
- Interleukins: IL-1, IL-2, IL-6, etc.
- Interferons: IFN-alpha, IFN-beta, IFN-gamma
- Signaling: NF-kappaB, TGF-beta, STAT1-5
- Immunology: MHC, TCR, BCR, CD4, CD8
- Metabolism: ATP, NADH, mTOR, AMPK
- Cell Biology: DNA, RNA, mRNA, miRNA

---

## Differential Expression

### create_standard_volcano()

**File:** `scripts/DE/plot_standard_volcano.R`

Create publication-quality volcano plot with FDR-based decisions.

```r
create_standard_volcano(
  de_results,
  decision_by = c("fdr", "p"),
  p_cutoff = 0.05,
  fc_cutoff = 2,
  top_n = 5,
  highlight_gene = NULL,
  label_method = "top",
  x_breaks = 1,
  title = "Volcano plot",
  subtitle = NULL,
  caption = NULL,
  color_palette = c(
    "NS" = "#7F7F7F",
    "Log2FC" = "#0173B2",
    "p-value" = "#029E73",
    "p-value & Log2FC" = "#D55E00"
  ),
  show_grid = FALSE,
  max.overlaps = 10
)
```

**Key Design Features:**

1. **FDR for decisions, raw-p for display:**
   - `decision_by = "fdr"`: Colors based on adj.P.Val
   - Y-axis shows -log10(P.Value) for resolution
   - Horizontal line at boundary where FDR = cutoff

2. **Priority labeling:**
   - `highlight_gene` genes use `max.overlaps = Inf`
   - Never suppressed by ggrepel overlap detection
   - Regular labels respect `max.overlaps` setting

**Returns:** ggplot object

---

### create_vertical_volcano()

**File:** `scripts/DE/volcano_helpers.R`

Create 90-degree rotated volcano for grid layouts.

```r
create_vertical_volcano(
  de_results,
  decision_by = "fdr",
  p_cutoff = 0.05,
  fc_cutoff = 2,
  ...
)
```

---

### combine_volcano_row()

**File:** `scripts/DE/volcano_helpers.R`

Combine multiple volcanoes with unified legend.

```r
combine_volcano_row(
  volcano_list,
  ncol = NULL,
  shared_legend = TRUE
)
```

---

### create_pca_plot()

**File:** `scripts/DE/plotPCA.R`

Create 2D PCA plot from DGEList.

```r
create_pca_plot(
  dge,
  color_by = "Group",
  shape_by = NULL,
  title = "PCA Plot",
  n_genes = 500
)
```

---

### create_3d_pca_plot()

**File:** `scripts/DE/plotPCA3d.R`

Create interactive 3D PCA with plotly.

```r
create_3d_pca_plot(
  dge,
  color_by = "Group",
  title = "3D PCA"
)
```

**Returns:** plotly object

---

## General Utilities

### annotate_genes_from_ensembl()

**File:** `scripts/General/annotate_genes.R`

Convert Ensembl IDs to gene symbols.

```r
annotate_genes_from_ensembl(
  ensembl_ids,
  species = "Mus musculus",
  strip_version = TRUE
)
```

---

### build_dge()

**File:** `scripts/General/dge_helpers.R`

Build validated DGEList with TMM normalization.

```r
build_dge(
  counts,
  samples,
  genes = NULL,
  normalize = TRUE
)
```

---

### read_counts_matrix()

**File:** `scripts/General/io_helpers.R`

Read count matrix from various formats.

```r
read_counts_matrix(
  file_path,
  format = "auto",
  gene_id_col = 1
)
```

---

## Shared Utilities

### custom_minimal_theme_with_grid()

**File:** `scripts/custom_minimal_theme.R`

Publication-ready ggplot2 theme.

```r
custom_minimal_theme_with_grid(base_size = 12)
```

**Features:**
- White background
- No panel grid
- Visible axis lines
- Clean, minimal appearance

---

### ensure_dir()

**File:** `scripts/utils_plotting.R`

Create directory if it doesn't exist.

```r
ensure_dir(path)
```

---

### save_plot()

**File:** `scripts/utils_plotting.R`

Save plot with proper device management.

```r
save_plot(
  plot,
  filename,
  width = 10,
  height = 8
)
```

---

### load_checkpoint()

**File:** `scripts/utils_plotting.R`

Load RDS checkpoint with validation.

```r
load_checkpoint(
  checkpoint_path,
  description = "checkpoint"
)
```

---

### log_message()

**File:** `scripts/utils_plotting.R`

Timestamped logging.

```r
log_message(message, level = "INFO")
```
