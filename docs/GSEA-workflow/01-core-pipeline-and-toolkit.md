# 01 - Core Pipeline and Toolkit Reference

**Scope:** How data flows from raw counts to master CSV tables in the 12868-EH project, covering the core pipeline (`1.1.core_pipeline.R`), configuration system, RNAseq-toolkit function inventory, checkpoint caching, and master table assembly.

**Last updated:** 2026-04-14

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Script Organization and Naming Conventions](#2-script-organization-and-naming-conventions)
3. [Configuration System](#3-configuration-system)
4. [Data Loading](#4-data-loading)
5. [limma-voom Differential Expression](#5-limma-voom-differential-expression)
6. [Checkpoint Caching Pattern](#6-checkpoint-caching-pattern)
7. [GSEA Execution](#7-gsea-execution)
8. [Master Table Assembly](#8-master-table-assembly)
9. [RNAseq-toolkit Function Inventory](#9-rnaseq-toolkit-function-inventory)
10. [Color and Theme System](#10-color-and-theme-system)

---

## 1. Architecture Overview

The pipeline follows a five-phase model where each phase depends on the previous:

```
Phase 1: Core Analysis (1.x.*.R)
  Data Load -> Filter/Normalize -> DE (limma-voom) -> GSEA (multi-database)
  Output: checkpoints/*.rds

Phase 2: Master Tables (1.5.create_master_tables.R)
  Checkpoints -> Normalized CSVs
  Output: tables/master_de_table.csv, master_gsea_table.csv

Phase 3-5: Visualization
  Master CSVs -> R plots (2.x.*), Python figures (3.x.*), Interactive (HTML)
```

**Central principle:** Expensive computations are cached as RDS checkpoints. All downstream visualization reads from master CSV tables (the "single source of truth"), never from raw RDS objects directly.

### Data Flow Diagram

```
featureCounts.merged.gene_counts.txt
        |
        v
  [1.1.core_pipeline.R]
        |
        +-> checkpoints/1.1_dge_raw.rds
        +-> checkpoints/1.1_dge_normalized.rds
        +-> checkpoints/1.1_fit_object.rds
        +-> checkpoints/1.1_de_results.rds
        +-> checkpoints/1.1_gsea_H_C2.rds
        +-> checkpoints/1.1_gsea_C5.rds
        +-> checkpoints/1.1_all_gsea_results.rds
        +-> plots/QC/*.pdf
        |
        v
  [1.5.create_master_tables.R]
        |
        +-> tables/master_de_table.csv
        +-> tables/master_gsea_table.csv
        +-> tables/master_gsea_significant.csv
        +-> tables/significant_de_genes.csv
        +-> tables/de_summary.csv
        +-> tables/gsea_summary_stats.csv
```

---

## 2. Script Organization and Naming Conventions

### Naming Pattern

```
{phase}.{order}.{description}.{ext}
```

| Phase | Purpose | Examples |
|-------|---------|---------|
| `1.x` | Core analysis, heavy computation | `1.1.core_pipeline.R`, `1.5.create_master_tables.R` |
| `2.x` | R visualization scripts | `2.1.volcano_plots.R`, `2.3.gsea_heatmaps.R` |
| `3.x` | Python visualization / interactive | `3.1.publication_figures.py`, `3.5.pathway_explorer.py` |

### Toolkit Directory Layout

```
01_scripts/RNAseq-toolkit/
  scripts/
    General/
      annotate_genes.R        # Ensembl -> Symbol/Entrez annotation
      dge_helpers.R           # DGEList construction (build_dge)
      io_helpers.R            # Count matrix I/O, metadata alignment
    DE/
      plot_standard_volcano.R # Standard + vertical volcano plots
      volcano_helpers.R       # Vertical volcano, combine_volcano_row
      plotPCA.R               # 2D PCA from DGEList
      plotPCA3d.R             # Interactive 3D PCA (plotly)
    GSEA/
      GSEA_processing/
        run_gsea.R                # Single-database GSEA wrapper
        run_gsea_analysis.R       # Multi-database pipeline + auto-plotting
        run_pooled_gsea.R         # Cross-contrast GSEA aggregation
        normalize_gsea.R          # gseaResult -> standardized tibble
        get_significant_pathways.R# Extract significant pathway IDs
        get_pathway_genes.R       # Leading edge genes (single result)
        get_pathway_genes_all.R   # Leading edge genes (multi-contrast)
        calculate_pathway_scores.R# Sample-level pathway scores
        parse_external_genesets.R # TransportDB, GMT, GMX, mitoXplorer parsers
        pathway_utils.R           # Size filtering, format conversion
      GSEA_plotting/
        gsea_dotplot.R            # Standard dotplot with NES gradient
        gsea_dotplot_facet.R      # Up/Down faceted dotplot
        gsea_barplot.R            # NES barplot
        gsea_running_sum_plot.R   # Three-panel running sum enrichment plot
        gsea_heatmap.R            # Expression z-score + NES heatmaps
        gsea_meta_heatmap.R       # Cross-database NES matrix heatmap
        gsea_plotting_utils.R     # save_gsea_plot, get_db_plot_params, smart_wrap
        format_pathway_names.R    # Smart biological name capitalization
        plot_pooled_contrast_dotplot.R
    ORA/
      run_ora.R                   # GO/KEGG over-representation analysis
      ora_dotplot.R               # ORA dotplot visualization
    custom_minimal_theme.R        # Shared ggplot2 theme
```

---

## 3. Configuration System

Configuration uses a two-layer system: a shared YAML file (for cross-language settings) and an R config file (for R-specific constants and functions).

### 3.1 pipeline.yaml (Shared Source of Truth)

**Path:** `02_analysis/config/pipeline.yaml`

Defines settings consumed by both R and Python scripts:

```yaml
project:
  id: "12868-EH"
  species: "Mus musculus"
  species_db: "MM"
  genome_build: "mm10"

colors:
  diverging:
    negative: "#2166AC"   # Blue (downregulated)
    neutral: "#F7F7F7"    # White
    positive: "#B35806"   # Orange (upregulated)
  groups:
    NTC: "#1f77b4"
    IL2RA_KO: "#ff7f0e"
  databases:              # Okabe-Ito colorblind-safe palette
    Hallmark: "#E69F00"
    KEGG: "#56B4E9"
    Reactome: "#009E73"
    # ... (13 database colors total)

analysis:
  de_fdr_cutoff: 0.05
  gsea_nperm: 100000
  gsea_seed: 123

schemas:
  master_gsea_table:
    required_columns:
      - pathway_id
      - pathway_name
      - database
      - nes
      - pvalue
      - padj
      - core_enrichment
```

### 3.2 config.R (R Configuration)

**Path:** `02_analysis/config/config.R`

Sources at the top of every analysis script via:

```r
source("02_analysis/config/config.R")
```

This file provides:

**Directory paths** -- All paths derived from `PROJECT_ROOT`:

```r
DIR_DATA_PROCESSED <- file.path(PROJECT_ROOT, "00_data/processed")
DIR_CHECKPOINTS    <- file.path(DIR_RESULTS, "checkpoints")
DIR_TABLES         <- file.path(DIR_RESULTS, "tables")
DIR_PLOTS          <- file.path(DIR_RESULTS, "plots")
DIR_TOOLKIT        <- file.path(DIR_SCRIPTS, "RNAseq-toolkit/scripts")
```

**Input file paths:**

```r
FILE_COUNTS_FEATURECOUNTS <- file.path(DIR_DATA_PROCESSED, "featureCounts.merged.gene_counts.txt")
FILE_MITOCARTA            <- file.path(DIR_DATA_REFERENCES, "mitochonria/Mouse.MitoCarta3.0.xls")
FILE_TRANSPORTDB          <- file.path(DIR_DATA_REFERENCES, "TransportDB2.0.csv")
```

**Experimental design:**

```r
SAMPLE_GROUPS <- list(
  NTC      = paste0("12868-EH-", sprintf("%04d", 1:5)),
  IL2RA_KO = paste0("12868-EH-", sprintf("%04d", 6:10))
)

CONTRASTS <- list(
  IL2RAKO_vs_NTC = c("IL2RA_KO", "NTC")  # c(numerator, denominator)
)
```

**Analysis parameters:**

```r
DE_FDR_CUTOFF      <- 0.05
DE_LOGFC_CUTOFF    <- 2.0
GSEA_PVALUE_CUTOFF <- 1.0     # Store all pathways, filter later
GSEA_NPERM         <- 100000
GSEA_SEED          <- 123
RANK_METRIC        <- "t"     # t-statistic from limma
```

**MSigDB database definitions:**

```r
MSIGDB_DATABASES <- list(
  H              = c("H", ""),              # Hallmark
  C2_KEGG        = c("C2", "CP:KEGG"),      # KEGG
  C2_REACTOME    = c("C2", "CP:REACTOME"),  # Reactome
  C2_WIKIPATHWAYS= c("C2", "CP:WIKIPATHWAYS"),
  C3_TF          = c("C3", "TFT:GTRD"),     # TF targets
  C5_BP          = c("C5", "GO:BP"),         # GO Biological Process
  C5_MF          = c("C5", "GO:MF"),         # GO Molecular Function
  C5_CC          = c("C5", "GO:CC")          # GO Cellular Component
)
```

**Checkpoint file naming:**

```r
CHECKPOINT_DGE_RAW        <- "1.1_dge_raw.rds"
CHECKPOINT_DGE_NORMALIZED <- "1.1_dge_normalized.rds"
CHECKPOINT_FIT_OBJECT     <- "1.1_fit_object.rds"
CHECKPOINT_DE_RESULTS     <- "1.1_de_results.rds"
CHECKPOINT_GSEA_H_C2      <- "1.1_gsea_H_C2.rds"
CHECKPOINT_GSEA_C5        <- "1.1_gsea_C5.rds"
CHECKPOINT_ALL_GSEA       <- "1.1_all_gsea_results.rds"
```

**Helper functions defined in config.R:**

- `load_or_compute()` -- Checkpoint caching (see [Section 6](#6-checkpoint-caching-pattern))
- `source_toolkit()` -- Sources all toolkit R scripts
- `load_packages()` -- Loads all required R packages

### 3.3 color_config.R (Color Palette System)

**Path:** `02_analysis/config/color_config.R`

Loads colors from `pipeline.yaml` when available, falls back to hardcoded values. Provides:

| Constant | Purpose | Example |
|----------|---------|---------|
| `DIVERGING_COLORS` | NES/logFC heatmaps (3-point) | `$negative="#2166AC"`, `$neutral="#F7F7F7"`, `$positive="#B35806"` |
| `DIVERGING_COLORS_5PT` | Fine-grained heatmaps (5-point) | Adds `$neg_weak`, `$pos_weak` |
| `GROUP_COLORS` | Experimental groups | `NTC="#1f77b4"`, `IL2RA_KO="#ff7f0e"` |
| `DATABASE_COLORS` | Pathway database categories (Okabe-Ito) | 13 named colors for Hallmark, KEGG, etc. |
| `ROBUSTNESS_TIER_COLORS` | GATOM sensitivity analysis | `core`, `moderate`, `exploratory` |
| `RUNNING_SUM_PALETTE` | Multi-pathway running sum plots | 9 distinct colors |

**Key functions:**

| Function | Returns | Use |
|----------|---------|-----|
| `nes_color_scale(limits, n_colors)` | `circlize::colorRamp2` | ComplexHeatmap NES coloring |
| `logfc_color_scale(limits)` | `circlize::colorRamp2` | ComplexHeatmap logFC coloring |
| `nes_ggplot_scale(limits, name)` | `scale_color_gradient2` | ggplot2 NES color aesthetic |
| `nes_ggplot_fill_scale(limits, name)` | `scale_fill_gradient2` | ggplot2 NES fill aesthetic |
| `get_diverging_palette(n)` | character vector | pheatmap-compatible palette |
| `get_correlation_palette(n)` | character vector | Sample correlation heatmaps |
| `zscore_color_scale(limits, for_pheatmap)` | colorRamp2 or character | Expression z-score heatmaps |
| `print_color_summary()` | console output | Debug/verify loaded colors |

---

## 4. Data Loading

### 4.1 Count Matrix Input

The pipeline reads featureCounts output with a special two-row header:

```
        Sample1  Sample2  ...
Group   NTC      NTC      ...
Gene1   1234     5678     ...
Gene2   ...      ...      ...
```

From `1.1.core_pipeline.R`, S1.1:

```r
# Read first two lines to get metadata
header_lines <- readLines(FILE_COUNTS_FEATURECOUNTS, n = 2)
sample_ids <- strsplit(header_lines[1], "\t")[[1]][-1]
group_assignments <- strsplit(header_lines[2], "\t")[[1]][-1]

# Read the actual data (skip first 2 rows)
counts_data <- read.delim(FILE_COUNTS_FEATURECOUNTS, skip = 2,
                          header = FALSE, check.names = FALSE)
colnames(counts_data) <- c("Gene_ID", sample_ids)
count_matrix <- as.matrix(counts_data[, -1])
rownames(count_matrix) <- counts_data$Gene_ID
```

### 4.2 Gene Annotation (Ensembl to Symbol)

Ensembl IDs are mapped to gene symbols using `org.Mm.eg.db`:

```r
ensembl_ids <- gsub("\\..*$", "", rownames(dge))  # Strip version suffix

gene_symbols <- mapIds(org.Mm.eg.db,
                      keys = ensembl_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# Duplicate symbols get Ensembl ID appended
new_rownames <- ifelse(!is.na(gene_symbols), gene_symbols, rownames(dge))
dup_symbols <- duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE)
new_rownames[dup_symbols] <- paste0(new_rownames[dup_symbols], "_", ensembl_ids[dup_symbols])
```

The gene annotation is stored in `dge$genes` as a data frame with columns: `ensembl_id`, `ensembl_id_base`, `symbol`.

### 4.3 DGEList Construction

The result is an edgeR `DGEList` object with:
- `dge$counts`: raw count matrix (genes x samples)
- `dge$samples`: sample metadata with `Sample`, `Group` columns
- `dge$genes`: gene annotation data frame

### 4.4 Toolkit I/O Helpers

The toolkit provides alternative data loading in `scripts/General/io_helpers.R`:

| Function | Purpose | Input |
|----------|---------|-------|
| `read_counts_matrix(fp)` | Read featureCounts/generic count files | TSV/CSV file path |
| `read_metadata(xlsx_fp)` | Read Excel metadata, standardize columns | Excel file path |
| `align_metadata_to_counts(md, counts_cols)` | Match metadata rows to count matrix columns | metadata df + column names |
| `write_annotated_matrix(mat, md, add_cols, outfile)` | Write annotated matrix with metadata header | matrix + metadata + annotations |

And in `scripts/General/dge_helpers.R`:

```r
build_dge(count_mat, samples_df, genes_df, round_nonint = TRUE, norm_method = "TMM")
```
- Rounds non-integer counts (warns user)
- Validates sample/gene order matches
- Applies TMM normalization

---

## 5. limma-voom Differential Expression

### 5.1 Filtering and Normalization (S1.2)

```r
# Filter lowly expressed genes using edgeR's filterByExpr
filter_design <- model.matrix(~ 0 + Group, data = dge$samples)
keep <- filterByExpr(dge, design = filter_design)
dge_filt <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization
dge_filt <- calcNormFactors(dge_filt, method = "TMM")

# Store log-CPM for visualization
dge_filt$logCPM <- cpm(dge_filt, log = TRUE, prior.count = 2)
```

### 5.2 Design Matrix and Contrasts (S1.4)

The pipeline uses a cell-means parameterization (no intercept):

```r
design <- model.matrix(~ 0 + Group, data = dge_normalized$samples)
colnames(design) <- levels(dge_normalized$samples$Group)
# Result: columns "NTC" and "IL2RA_KO"
```

voom transformation with quality weights:

```r
v <- voomWithQualityWeights(dge_normalized, design, plot = FALSE)
```

Contrasts are built from the `CONTRASTS` config:

```r
# CONTRASTS = list(IL2RAKO_vs_NTC = c("IL2RA_KO", "NTC"))
# Produces: "IL2RA_KO - NTC"
contrast_strings <- sapply(CONTRASTS, function(x) paste(x[1], "-", x[2]))
contrast_matrix <- makeContrasts(contrasts = contrast_strings, levels = design)
```

### 5.3 Model Fitting

```r
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
```

The returned object is a list containing:
- `fit`: the eBayes-moderated fit object (`MArrayLM`)
- `voom`: the voom-transformed expression object
- `design`: the design matrix
- `contrast_matrix`: the contrast matrix

### 5.4 Extracting DE Results (S1.4b)

```r
results <- topTable(fit_results$fit, coef = contrast_name,
                    number = Inf, sort.by = "none")
```

Returns a data frame with columns: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`. Gene symbols are in the rownames.

---

## 6. Checkpoint Caching Pattern

Every expensive computation uses `load_or_compute()`, defined in `config.R`:

```r
load_or_compute <- function(checkpoint_file,
                            compute_fn,
                            force_recompute = FALSE,
                            description = "Result") {

  # Resolve relative paths against DIR_CHECKPOINTS
  if (!grepl("^/", checkpoint_file) && !grepl("^[A-Z]:", checkpoint_file)) {
    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)
  } else {
    checkpoint_path <- checkpoint_file
  }

  # Load from cache if available
  if (file.exists(checkpoint_path) && !force_recompute) {
    message(sprintf("[CACHE] Loading %s from: %s", description, checkpoint_file))
    return(readRDS(checkpoint_path))
  }

  # Otherwise compute, time it, and save
  message(sprintf("[COMPUTE] Computing %s...", description))
  start_time <- Sys.time()
  result <- compute_fn()
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)
  message(sprintf("[SAVE] Saving %s to: %s (took %.2f min)", description, checkpoint_file, elapsed))
  saveRDS(result, checkpoint_path)

  return(result)
}
```

**Usage pattern** (from the core pipeline):

```r
dge_normalized <- load_or_compute(
  checkpoint_file = CHECKPOINT_DGE_NORMALIZED,
  description = "DGEList (filtered & normalized)",
  compute_fn = function() {
    # ... expensive filtering/normalization code ...
    return(dge_filt)
  }
)
```

**Characteristics:**
- Checkpoint file names are relative to `DIR_CHECKPOINTS` by default
- Absolute paths are also accepted
- `force_recompute = TRUE` bypasses the cache
- Timing is logged automatically
- The `compute_fn` is a zero-argument function (closure) that captures variables from the enclosing scope

### Checkpoint Files Produced by the Core Pipeline

| Checkpoint | Content | Approximate Size |
|-----------|---------|-----------------|
| `1.1_dge_raw.rds` | Raw DGEList with gene annotations | Small |
| `1.1_dge_normalized.rds` | Filtered + TMM-normalized DGEList with logCPM | Small |
| `1.1_fit_object.rds` | List: `fit`, `voom`, `design`, `contrast_matrix` | Medium |
| `1.1_de_results.rds` | Named list of DE result data frames (one per contrast) | Medium |
| `1.1_gsea_H_C2.rds` | Named list of gseaResult objects (H, C2_KEGG, C2_REACTOME, C2_WIKIPATHWAYS, C3_TF) | Large |
| `1.1_gsea_C5.rds` | Named list of gseaResult objects (C5_BP, C5_MF, C5_CC) | Large |
| `1.1_all_gsea_results.rds` | Combined named list of all GSEA results | Large |

---

## 7. GSEA Execution

### 7.1 How GSEA Runs in the Core Pipeline

The pipeline runs GSEA in two batches:

**Batch 1 (S1.5):** Hallmark + C2 + C3

```r
gsea_h_c2 <- load_or_compute(
  checkpoint_file = CHECKPOINT_GSEA_H_C2,
  compute_fn = function() {
    db_names <- c("H", "C2_KEGG", "C2_REACTOME", "C2_WIKIPATHWAYS", "C3_TF")
    # For each database, call run_gsea() with params from MSIGDB_DATABASES
    ...
  }
)
```

**Batch 2 (S1.6):** GO terms (C5)

```r
gsea_c5 <- load_or_compute(
  checkpoint_file = CHECKPOINT_GSEA_C5,
  compute_fn = function() {
    db_names <- c("C5_BP", "C5_MF", "C5_CC")
    ...
  }
)
```

**Combine (S1.7):**

```r
all_gsea <- c(gsea_h_c2, gsea_c5)
# Result: named list with keys H, C2_KEGG, C2_REACTOME, ... C5_CC
```

### 7.2 The run_gsea() Function

**File:** `scripts/GSEA/GSEA_processing/run_gsea.R`

This is the core GSEA wrapper. It:
1. Prepares a ranked gene list from DE results (sorted by the `rank_metric` column, default `"t"`)
2. Fetches gene sets from MSigDB via `msigdbr`
3. Runs `clusterProfiler::GSEA()` with `fgsea` backend

```r
run_gsea(
  DE_results,             # data.frame with gene symbols as rownames, must have rank_metric column
  rank_metric   = "t",    # Column to rank by (t-statistic)
  species       = "Mus musculus",
  db_species    = NULL,   # "MM" or "HS" (if NULL, uses species-only mode)
  collection    = "H",    # MSigDB collection code
  subcollection = "",     # MSigDB subcollection (e.g., "CP:KEGG")
  pvalue_cutoff = 1,      # 1 = store all pathways
  padj_method   = "fdr",
  nperm         = 100000,
  seed          = 123
)
# Returns: gseaResult object from clusterProfiler
```

**msigdbr version compatibility:** The function auto-detects whether the installed msigdbr uses the v7.5 API (`category`/`subcategory`) or the v8+ API (`collection`/`subcollection`), and calls accordingly.

### 7.3 Ranked Gene List Construction

Inside `run_gsea()`:

```r
gene_vector <- DE_results[[rank_metric]]   # Extract ranking column
names(gene_vector) <- rownames(DE_results) # Assign gene names
ranked_genes <- sort(gene_vector, decreasing = TRUE)  # Sort descending
# NAs and Inf values are removed
```

The function expects **gene symbols** as rownames of `DE_results` (not Ensembl IDs).

---

## 8. Master Table Assembly

### 8.1 Script: 1.5.create_master_tables.R

This script transforms cached RDS checkpoint objects into standardized CSV tables that all downstream scripts consume.

**Inputs:** Reads from `DIR_CHECKPOINTS`:
- `1.1_dge_raw.rds` (for Ensembl ID annotations)
- `1.1_de_results.rds`
- `1.1_all_gsea_results.rds`
- `1.4_gsea_mito_unified.rds` (optional)

**Key dependency:** Sources `02_analysis/helpers/normalize_gsea.R` for the `normalize_gsea_results()` function.

### 8.2 Master DE Table

**Output:** `tables/master_de_table.csv`

**Column Schema:**

| Column | Type | Description |
|--------|------|-------------|
| `contrast` | character | Contrast formula (e.g., "IL2RA_KO - NTC") |
| `gene_symbol` | character | Gene symbol (from rownames of DE results) |
| `ensembl_id` | character | Full Ensembl ID (with version) |
| `ensembl_id_base` | character | Ensembl ID without version suffix |
| `logFC` | numeric | Log2 fold change |
| `AveExpr` | numeric | Average expression (log2-CPM) |
| `t` | numeric | Moderated t-statistic |
| `P.Value` | numeric | Raw p-value |
| `adj.P.Val` | numeric | BH-adjusted p-value (FDR) |
| `B` | numeric | Log-odds of differential expression |

**Assembly logic:**

```r
for (contrast in names(de_results)) {
  de_table <- de_results[[contrast]]
  de_table$gene_symbol <- rownames(de_table)
  # Attach Ensembl IDs from dge_raw$genes
  de_table$ensembl_id <- dge_raw$genes$ensembl_id[match(rownames(de_table), rownames(dge_raw))]
  de_table$contrast <- contrast
  master_de_list[[contrast]] <- de_table
}
master_de <- bind_rows(master_de_list)
```

**Derived tables:**
- `significant_de_genes.csv` -- Subset where `adj.P.Val < 0.05`, sorted by significance
- `de_summary.csv` -- Per-contrast counts of total, significant up, and significant down genes

### 8.3 Master GSEA Table

**Output:** `tables/master_gsea_table.csv`

**Column Schema (from `normalize_gsea_results()`):**

| Column | Type | Description |
|--------|------|-------------|
| `pathway_id` | character | Original pathway ID (e.g., "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |
| `pathway_name` | character | Cleaned/formatted pathway name |
| `database` | character | Display name: "Hallmark", "KEGG", "Reactome", "GO_BP", etc. |
| `contrast` | character | Contrast name from config |
| `nes` | numeric | Normalized Enrichment Score |
| `pvalue` | numeric | Raw p-value |
| `padj` | numeric | Adjusted p-value (FDR) |
| `set_size` | integer | Number of genes in the gene set |
| `leading_edge_size` | integer | Number of genes in the leading edge |
| `gene_ratio` | numeric | `leading_edge_size / set_size` |
| `core_enrichment` | character | Slash-separated leading edge gene symbols |
| `neg_log_padj` | numeric | `-log10(padj)`, capped at 16 |
| `direction` | character | "Up" or "Down" based on NES sign |

**Assembly logic:**

```r
# Database display name mapping
DB_DISPLAY_NAMES <- c(
  H = "Hallmark", C2_KEGG = "KEGG", C2_REACTOME = "Reactome",
  C2_WIKIPATHWAYS = "WikiPathways", C5_BP = "GO_BP",
  C5_MF = "GO_MF", C5_CC = "GO_CC", Mitochondria = "Mitochondria"
)

for (db_name in names(all_gsea)) {
  normalized <- normalize_gsea_results(
    gsea_obj    = all_gsea[[db_name]],
    database    = DB_DISPLAY_NAMES[db_name],
    contrast    = contrast_name,
    padj_cutoff = 1,            # Keep all pathways, not just significant
    strip_prefix = TRUE
  )
  all_normalized[[db_name]] <- normalized
}
master_gsea <- bind_rows(all_normalized)
```

**Derived tables:**
- `master_gsea_significant.csv` -- Subset where `padj < 0.05`, sorted by database then significance
- `gsea_summary_stats.csv` -- Per-database/contrast/direction counts and mean |NES|

### 8.4 The normalize_gsea_results() Function

There are two copies of this function:
1. **Toolkit version** at `scripts/GSEA/GSEA_processing/normalize_gsea.R` -- uses `NES` (uppercase) column name
2. **Project helper** at `02_analysis/helpers/normalize_gsea.R` -- uses `nes` (lowercase) column name

The project helper is what `1.5.create_master_tables.R` sources. Both share the same logic:

1. Accept a `gseaResult` object or a data frame
2. Detect the adjusted p-value column (`p.adjust`, `qvalue`, or `padj`)
3. Calculate `leading_edge_size` by counting `/`-separated genes in `core_enrichment`
4. Calculate `gene_ratio = leading_edge_size / setSize`
5. Clean pathway names (strip prefixes, title case, truncate)
6. Add `direction` ("Up"/"Down") and `neg_log_padj`

---

## 9. RNAseq-toolkit Function Inventory

### 9.1 GSEA Processing Functions

#### `run_gsea()`
**File:** `scripts/GSEA/GSEA_processing/run_gsea.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `DE_results` | data.frame | (required) | DE results with gene symbols as rownames |
| `rank_metric` | character | `"t"` | Column to rank genes by |
| `species` | character | `"Mus musculus"` | Species for msigdbr |
| `db_species` | character | `NULL` | Species abbreviation ("MM", "HS") |
| `collection` | character | `"H"` | MSigDB collection code |
| `subcollection` | character | `""` | MSigDB subcollection code |
| `pvalue_cutoff` | numeric | `1` | P-value cutoff for result storage |
| `padj_method` | character | `"fdr"` | P-value adjustment method |
| `nperm` | integer | `100000` | Number of fgsea permutations |
| `seed` | integer | `123` | Random seed |

**Returns:** `gseaResult` object from clusterProfiler.

#### `run_gsea_analysis()`
**File:** `scripts/GSEA/GSEA_processing/run_gsea_analysis.R`

High-level pipeline that runs GSEA across multiple databases and generates plots automatically.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `de_table` | data.frame | (required) | DE results table |
| `analysis_name` | character | (required) | Name for plot labeling |
| `rank_metric` | character | `"t"` | Column to rank genes by |
| `species` | character | `"Mus musculus"` | Species |
| `n_pathways` | integer | `30` | Top pathways to display in plots |
| `padj_cutoff` | numeric | `0.05` | FDR cutoff for significance |
| `save_plots` | logical | `TRUE` | Whether to save generated plots |
| `output_dir` | character | `"./GSEA_Plots"` | Output directory for plots |
| `databases` | list | `NULL` | Custom database list (NULL = defaults) |
| `nperm` | integer | `100000` | Permutations for GSEA |
| `pvalue_cutoff` | numeric | `1` | Store all pathways |
| `sample_annotation` | data.frame | `NULL` | For sample-level heatmaps |
| `sample_order` | character | `NULL` | Custom sample ordering |
| `helper_root` | character | `NULL` | Root dir for helper scripts |

**Returns:** Named list of gseaResult objects (one per database).

**Auto-generated plots per database:** dotplot (up), dotplot (down), faceted dotplot, NES barplot, running sum (top 5), heatmap (if sample_annotation provided).

#### `run_pooled_gsea()`
**File:** `scripts/GSEA/GSEA_processing/run_pooled_gsea.R`

Cross-contrast GSEA aggregation. Runs GSEA per contrast, pools significant pathways, extracts leading edge genes, and calculates pathway scores.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `fit` | MArrayLM | (required) | limma fit object |
| `contrasts` | matrix | (required) | Contrast matrix |
| `DGEobject` | DGEList | (required) | For pathway score calculation |
| `species` | character | `"Mus musculus"` | Species |
| `top_n` | integer | `30` | Top pathways for gene extraction |
| `padj_cutoff` | numeric | `0.05` | Significance threshold |
| `gsea_pvalue_cutoff` | numeric | `1` | GSEA storage cutoff |
| `rank_metric` | character | `"t"` | Ranking column |
| `nperm` | integer | `100000` | Permutations |
| `databases` | list | `NULL` | Custom databases |
| `helper_root` | character | `NULL` | Script root directory |
| `cached_gsea_results` | list | `NULL` | Pre-computed results to skip recomputation |
| `verbose` | logical | `TRUE` | Print progress |
| `log_file` | character | `NULL` | Path for log output |

**Returns:** List with four elements:
- `gsea_results`: Nested list `contrast -> database -> gseaResult`
- `pools`: `database -> character vector of significant pathway IDs`
- `genes`: `database -> pathway_id -> character vector of genes`
- `scores`: `database -> matrix (samples x pathways)`

#### `normalize_gsea_results()`
**File:** `scripts/GSEA/GSEA_processing/normalize_gsea.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult or data.frame | (required) | GSEA result to normalize |
| `database` | character | (required) | Database label |
| `contrast` | character | (required) | Contrast label |
| `padj_cutoff` | numeric | `1` | Filter threshold (1 = keep all) |
| `format_names` | logical | `TRUE` | Apply smart name formatting |
| `max_name_length` | integer | `80` | Truncate names beyond this |

**Returns:** Tibble with standardized schema (see [Section 8.3](#83-master-gsea-table)).

#### `empty_gsea_tibble()`
Returns an empty tibble with the correct column schema for `normalize_gsea_results()`.

#### `summarize_gsea_results(gsea_df, padj_cutoff = 0.05)`
Returns a wide summary tibble with counts and mean NES per database/contrast/direction.

#### `clean_pathway_name_basic(names, max_length, strip_prefix)`
Fallback name cleaner: strips prefixes, replaces underscores, applies title case.

#### `get_significant_pathways()`
**File:** `scripts/GSEA/GSEA_processing/get_significant_pathways.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_results_list` | list | (required) | List of gseaResult objects |
| `padj_cutoff` | numeric | `0.05` | FDR cutoff |
| `verbose` | logical | `FALSE` | Print per-result stats |

**Returns:** Character vector of unique significant pathway IDs across all results.

#### `get_pathway_genes()`
**File:** `scripts/GSEA/GSEA_processing/get_pathway_genes.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult | (required) | Single GSEA result |
| `padj_cutoff` | numeric | `0.05` | FDR cutoff |
| `top` | integer or NULL | `NULL` | Limit to top N pathways |

**Returns:** Named list: `pathway_id -> character vector of core enrichment genes`.

#### `get_pathway_genes_all()`
**File:** `scripts/GSEA/GSEA_processing/get_pathway_genes_all.R`

Cross-contrast version. Finds pathways significant in any contrast, ranks by minimum p-value.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_results_list` | nested list | (required) | `contrast -> database -> gseaResult` |
| `database` | character | (required) | Database name to extract |
| `padj_cutoff` | numeric | `0.05` | FDR cutoff |
| `top` | integer or NULL | `NULL` | Limit to top N pathways |

**Returns:** Named list: `pathway_id -> character vector of genes`.

#### `calculate_pathway_scores()`
**File:** `scripts/GSEA/GSEA_processing/calculate_pathway_scores.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expression_data` | matrix | (required) | Normalized expression (genes x samples) |
| `pathway_genes` | named list | (required) | `pathway_id -> gene vector` |
| `method` | character | `"mean"` | Scoring method |
| `verbose` | logical | `FALSE` | Print progress |

**Returns:** Numeric matrix (samples x pathways) of pathway activity scores.

### 9.2 External Gene Set Parsers

**File:** `scripts/GSEA/GSEA_processing/parse_external_genesets.R`

All parsers return a list with: `T2G` (TERM2GENE data frame), `T2N` (TERM2NAME data frame), `stats`, `source`, `created`.

#### `parse_transportdb(file, prefix, group_by, org_db, id_type, min_size, max_size)`
Parses TransportDB2.0 CSV. Groups transporter genes by family. Returns TERM2GENE format for clusterProfiler.

#### `parse_geneset_file(file, gs_col, gene_col, desc_col, prefix, sep, header)`
Generic parser for any two-column or multi-column gene set file. Auto-detects separator.

#### `parse_gmt(file, prefix)`
Parses GMT (Gene Matrix Transposed) format: `gene_set_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...`

#### `parse_gmx(file, prefix, min_size, max_size)`
Parses GMX format (columns are gene sets, rows 1-2 are descriptions/names, rows 3+ are genes).

#### `parse_mitoxplorer(file, prefix, gene_col, process_col, min_size, max_size)`
Parses mitoXplorer3.0 `mouse_gene_function.txt`. Groups genes by `mito_process`.

#### `convert_geneset_ids(T2G, org_db, from_type, to_type, drop_unmapped)`
Converts gene IDs in a T2G data frame using an AnnotationDbi database. Supports Ensembl, Symbol, Entrez.

#### `convert_human_to_mouse(T2G, drop_unmapped, verbose)`
Uses the `homologene` package to convert human gene symbols to mouse orthologs.

### 9.3 Pathway Utility Functions

**File:** `scripts/GSEA/GSEA_processing/pathway_utils.R`

#### `filter_pathways_by_size(T2G, T2N, min_size, max_size, verbose)`
Filters gene sets to retain only those within size bounds. Returns list with filtered `T2G`, `T2N`, and counts.

#### `get_pathway_size_stats(T2G)`
Returns named vector: `n_pathways`, `min_size`, `q1_size`, `median_size`, `mean_size`, `q3_size`, `max_size`, `n_genes`.

#### `list_to_term2gene(geneset_list, gs_col, gene_col)`
Converts a named list of gene sets to a TERM2GENE data frame.

#### `term2gene_to_list(T2G)`
Inverse of above: converts TERM2GENE to a named list.

#### `create_term2name(T2G, format_names)`
Creates a TERM2NAME data frame from TERM2GENE, using gene set names as descriptions with optional formatting.

#### `export_to_gmx(T2G, T2N, output_file)`
Exports gene sets to GMX format for use with external GSEA tools.

### 9.4 GSEA Plotting Functions

#### `gsea_dotplot()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_dotplot.R`

Creates a dotplot with continuous NES gradient coloring (Blue-White-Orange). Separates pathway **selection** (what to show) from **highlighting** (significance).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult | (required) | GSEA result object |
| `filterBy` | character | `"p.adjust"` | Sort method: "p.adjust", "NES", "NES_positive", "NES_negative" |
| `sortBy` | character | `"GeneRatio"` | Display sort: "GeneRatio" or "p.adjust" |
| `showCategory` | integer | `10` | Number of pathways to show |
| `padj_cutoff` | numeric | `0.05` | Threshold for black outline highlighting |
| `title` | character | `"GSEA Dotplot"` | Plot title |
| `wrap_width` | integer | `50` | Text wrapping width |
| `neg_color` | character | `"#2166AC"` | Blue for negative NES |
| `mid_color` | character | `"#F7F7F7"` | White for zero NES |
| `pos_color` | character | `"#B35806"` | Orange for positive NES |
| `highlight_sig` | logical | `TRUE` | Add black outline to significant pathways |
| `highlight_threshold` | numeric | `NULL` | Override highlight cutoff |
| `strip_prefix` | logical | `TRUE` | Remove database prefixes |
| `use_gradient` | logical | `TRUE` | Continuous NES gradient vs binary colors |

**Returns:** ggplot2 object.

**Design pattern:** Base points have `stroke = 0` with `color = "transparent"`. Significant points get a second `geom_point` overlay with `stroke = 2, color = "black"`.

#### `gsea_dotplot_facet()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R`

Faceted version split by Up/Down direction. Shows top N pathways per direction.

Same parameters as `gsea_dotplot()` except `filterBy` is not available (both directions shown).

#### `gsea_barplot()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_barplot.R`

Horizontal barplot of NES values with continuous gradient fill.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult | (required) | GSEA result object |
| `padj_cutoff` | numeric | `0.05` | Only significant pathways |
| `top_n` | integer | `30` | Number to show |
| `title` | character | `"GSEA NES Barplot"` | Plot title |
| `nes_limits` | numeric(2) | `c(-3.5, 3.5)` | Color scale limits |

**Returns:** ggplot2 object.

#### `gsea_running_sum_plot()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R`

Three-panel running sum enrichment plot using `enrichplot::gseaplot2()` with patchwork layout.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult | (required) | GSEA result |
| `gene_set_ids` | integer or character | `NULL` | Pathway indices or IDs (NULL = top 5 by |NES|) |
| `palette` | character vector | `NULL` | Colors (NULL = 9-color default) |
| `labels` | named character | `NULL` | Custom legend labels |
| `legend_pos` | numeric(2) | `c(.98, .98)` | Legend position |
| `base_size` | numeric | `14` | Font size |
| `max_name_length` | integer | `40` | Truncate legend labels |
| `title` | character | `NULL` | Plot title |

**Returns:** patchwork object (three vertically stacked panels).

**Critical note:** The color palette vector must NOT be named. Named colors cause issues with `enrichplot::gseaplot2()` for custom databases.

#### `gsea_heatmap()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_heatmap.R`

Expression z-score heatmap for leading edge genes of significant pathways.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_obj` | gseaResult | (required) | GSEA result |
| `expr_data` | matrix | (required) | Expression matrix (genes x samples) |
| `sample_annotation` | data.frame | (required) | Sample metadata for column annotation |
| `ann_colors` | list | `NULL` | Annotation colors |
| `padj_cutoff` | numeric | `0.05` | Significance cutoff |
| `n_pathways` | integer | `20` | Number of pathways |
| `max_genes_per_pathway` | integer | `25` | Gene limit per pathway |
| `title` | character | `"GSEA Results Heatmap"` | Title |

**Returns:** Invisible list with `matrix`, `heatmap`, `genes`.

Falls back to `gsea_heatmap_nes()` (NES values replicated across samples) when expression data is missing.

#### `gsea_to_matrix()` / `gsea_heatmap_save()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_meta_heatmap.R`

Builds a NES matrix across databases and saves as pheatmap PDF.

`gsea_to_matrix(gsea_results, n_pathways, padj_cutoff)` returns a matrix where rows are pathways and columns are databases.

`gsea_heatmap_save(nes_matrix, file, annotation_col, ...)` saves the heatmap as PDF.

#### `format_pathway_name()`
**File:** `scripts/GSEA/GSEA_plotting/format_pathway_names.R`

Smart biological name formatting with a three-step process:
1. Replace underscores with spaces
2. Strip database prefixes (HALLMARK_, KEGG_, REACTOME_, etc.)
3. Apply smart capitalization preserving abbreviations (NF-kappaB, IL-2, STAT5, MHC, etc.)

Maintains dictionaries for:
- ~100 single-word biological abbreviations (`build_exception_dictionary()`)
- ~25 multi-word patterns (`build_multiword_patterns()`) like "nf kappa b" -> "NF-kappaB"
- Roman numerals (I through X)
- Greek letters (kept lowercase)
- Chemical prefixes (n-, cis-, trans-, etc.)
- Conjunctions kept lowercase (via, and, of, in, to, by, from, the)

### 9.5 Plotting Utilities

#### `save_gsea_plot()`
**File:** `scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R`

Saves a ggplot to PDF with dynamic font scaling based on plot dimensions.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `plot` | ggplot | (required) | Plot to save |
| `filename` | character | (required) | Output filename |
| `width` | numeric | (required) | Width in inches |
| `height` | numeric | (required) | Height in inches |
| `base_font_size` | numeric | `10` | Base font size |
| `dir` | character | `NULL` | Output directory |
| `dpi` | integer | `300` | Resolution |

Uses `on.exit(dev.off())` to guarantee the PDF device is closed.

#### `get_db_plot_params(db_name)`
Returns a list with `width`, `height`, `font_size` tuned per database type (hallmark, gobp, reactome, etc.).

#### `smart_wrap(text, width = 40)`
Splits text at word boundaries: two-part split for moderately long text, three-part split for very long text.

#### `save_gsea_log(gsea_obj, filename, padj_cutoff, dir)`
Writes a human-readable text log of GSEA results with summary statistics and per-pathway details.

### 9.6 DE Visualization Functions

#### `create_standard_volcano()`
**File:** `scripts/DE/plot_standard_volcano.R`

Standard volcano plot with `-log10(P.Value)` on y-axis. Significance decisions use either FDR or raw p-value.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `de_results` | data.frame | (required) | Must have `logFC`, `P.Value`, `adj.P.Val` columns |
| `decision_by` | character | `"fdr"` | "fdr" or "p" -- which column determines colors |
| `p_cutoff` | numeric | `0.05` | Significance threshold (FDR or p depending on `decision_by`) |
| `fc_cutoff` | numeric | `2` | Log2 fold-change threshold |
| `top_n` | integer | `5` | Top genes to label per side |
| `highlight_gene` | character | `NULL` | Extra genes to label (bold, black, never suppressed) |
| `label_method` | character | `"top"` | "top", "sig", "p", "log2fc" |
| `color_palette` | named character(4) | Grey/Blue/Green/Orange | Colors for NS, Log2FC, p-value, Both |
| `show_grid` | logical | `FALSE` | Show panel grid |
| `max.overlaps` | integer | `10` | ggrepel overlap limit |
| `subtitle` | character | `NULL` | Optional subtitle |

**Returns:** ggplot2 object.

**Key design:** When `decision_by = "fdr"`, the horizontal dashed line is placed at the raw p-value corresponding to the FDR boundary (the largest raw p among genes passing FDR). This avoids staircase artifacts from plotting -log10(FDR) directly.

#### `create_vertical_volcano()`
**File:** `scripts/DE/volcano_helpers.R`

90-degree rotated volcano: `-log10(p)` on x-axis, `logFC` on y-axis. Useful for stacked panel layouts.

Same parameters as `create_standard_volcano()`.

#### `combine_volcano_row()`
**File:** `scripts/DE/volcano_helpers.R`

Combines multiple vertical volcano plots into a row with unified scales and legend.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `volcano_list` | list of ggplots | (required) | From `create_vertical_volcano()` |
| `labels` | character | `names(volcano_list)` | Panel labels |
| `guide_position` | character | `"bottom"` | Legend position |
| `keep_first_caption` | logical | `FALSE` | Keep caption under first panel |

#### `create_pca_plot()`
**File:** `scripts/DE/plotPCA.R`

2D PCA plot from a DGEList object using ggplot2.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `DGE_object` | DGEList | (required) | Must have `group` and `organ` in `$samples` |
| `title` | character | `"PCA Plot"` | Title |
| `xlim_abs` | numeric | `NULL` | X-axis absolute limit |
| `ylim_abs` | numeric | `NULL` | Y-axis absolute limit |
| `point_size` | numeric | `5` | Point size |
| `label_size` | numeric | `4` | Label size |

**Returns:** ggplot2 object with `coord_fixed()`.

#### `annotate_genes_from_ensembl()`
**File:** `scripts/General/annotate_genes.R`

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ens_ids` | character | (required) | Ensembl gene IDs |
| `try_biomart` | logical | `TRUE` | Try biomaRt for biotype info |

**Returns:** Tibble with `Symbol`, `Ensembl`, `ENTREZID`, `gene_biotype`.

### 9.7 ORA Functions

#### `run_ora(gene_list, species, ont, pvalue_cutoff, qvalue_cutoff, min_gs_size, max_gs_size)`
**File:** `scripts/ORA/run_ora.R`

GO over-representation analysis using `clusterProfiler::enrichGO()`. Returns `enrichResult` object.

#### `run_ora_kegg(gene_list, species, pvalue_cutoff, qvalue_cutoff)`
KEGG ORA. Converts symbols to Entrez IDs internally.

#### `normalize_ora_results(ora_result, database, module, padj_cutoff, format_names, max_name_length)`
Converts `enrichResult` to standardized tibble with columns: `pathway_id`, `pathway_name`, `database`, `module`, `gene_ratio`, `pvalue`, `padj`, `gene_count`, `genes`.

#### `run_ora_all_ontologies(gene_list, species, ...)`
Convenience wrapper: runs GO BP, MF, CC and returns named list.

#### `ora_dotplot(ora_result, top_n, sort_by, padj_cutoff, ...)`
**File:** `scripts/ORA/ora_dotplot.R`

Dotplot for ORA results. Uses viridis color scale when available.

#### `ora_dotplot_facet(ora_list, top_n_per_module, padj_cutoff, ...)`
Multi-module faceted ORA dotplot.

### 9.8 Theme Function

#### `custom_minimal_theme_with_grid()`
**File:** `scripts/custom_minimal_theme.R`

Clean ggplot2 theme based on `theme_classic()` with:
- White background, no grid lines
- Black axis lines and ticks (0.5 linewidth)
- No panel border
- Centered plot title
- 10px margins on all sides

---

## 10. Color and Theme System

### Color Palette Summary

| Palette | Usage | Source |
|---------|-------|--------|
| Blue-White-Orange diverging | NES, logFC, activity scores | `pipeline.yaml` / `color_config.R` |
| Group colors (Blue/Orange) | NTC vs IL2RA_KO | `config.R` / `color_config.R` |
| Okabe-Ito (13 colors) | Database categories | `pipeline.yaml` / `color_config.R` |
| 9-color vibrant | Running sum plots | `color_config.R` |
| Blue tiers (3 levels) | GATOM robustness | `color_config.R` |

### Colorblind Safety

All palettes are designed for deuteranopia and protanopia accessibility:
- Diverging scale avoids red-green: uses blue (#2166AC) and orange (#B35806)
- Database palette uses Okabe-Ito scheme validated for color vision deficiency
- Group colors use blue (#1f77b4) and orange (#ff7f0e), distinguishable by luminance

### Where Colors Are Applied

| Visualization | Color Source | Function |
|---------------|-------------|----------|
| GSEA dotplot fill | `neg_color`/`mid_color`/`pos_color` params | `gsea_dotplot()` |
| GSEA barplot fill | Same gradient params | `gsea_barplot()` |
| ComplexHeatmap NES | `nes_color_scale()` | Returns `colorRamp2` |
| pheatmap diverging | `get_diverging_palette()` | Returns hex vector |
| ggplot2 NES color | `nes_ggplot_scale()` | Returns `scale_color_gradient2` |
| Volcano plot categories | `color_palette` param | `create_standard_volcano()` |
| Sample correlation | `get_correlation_palette()` | White-to-teal gradient |
| Z-score heatmap | `zscore_color_scale()` | Blue-white-orange |

### Loading Colors in Scripts

**R scripts:**
```r
source("02_analysis/config/color_config.R")
# Now DIVERGING_COLORS, GROUP_COLORS, DATABASE_COLORS, etc. are available
```

**Python scripts:** Read directly from `pipeline.yaml` using PyYAML.

---

## Appendix: Quick Resume Patterns

### Restart from Checkpoints

```r
source("02_analysis/config/config.R")

# Load latest normalized DGE
dge <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))

# Load fit object
fit_results <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_fit_object.rds"))

# Load all GSEA results
all_gsea <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_all_gsea_results.rds"))
```

### Force Recompute a Single Step

```r
de_results <- load_or_compute(
  checkpoint_file = CHECKPOINT_DE_RESULTS,
  description = "DE results (all contrasts)",
  force_recompute = TRUE,   # <-- Force fresh computation
  compute_fn = function() { ... }
)
```

### Add a New GSEA Database

1. Add entry to `MSIGDB_DATABASES` in `config.R`
2. Add display name to `DB_DISPLAY_NAMES` in `1.5.create_master_tables.R`
3. Add database color to `DATABASE_COLORS` in `color_config.R` and `pipeline.yaml`
4. Rerun the relevant GSEA checkpoint and master tables
