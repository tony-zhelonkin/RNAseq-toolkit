# MSigDB GSEA Pipeline Reference

## Overview

This document describes the complete Gene Set Enrichment Analysis (GSEA) pipeline used in the RNAseq-toolkit, covering gene set loading from MSigDB via the `msigdbr` R package, GSEA execution via `clusterProfiler::GSEA` (which delegates to the `fgsea` algorithm), result normalization into a standard tabular schema, checkpoint caching, and the master table output format. All code snippets reference actual source files.

**Key source files:**

| File | Purpose |
|------|---------|
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R` | Core single-database GSEA wrapper |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea_analysis.R` | Multi-database pipeline with auto-plotting |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R` | Result normalization to standard tibble |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/parse_external_genesets.R` | Parsers for non-MSigDB gene sets (TransportDB, GMT, GMX, mitoXplorer) |
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/pathway_utils.R` | Gene set utility functions (filtering, conversion, export) |
| `02_analysis/config/config.R` | Project-level GSEA configuration (databases, parameters) |
| `02_analysis/config/pipeline.yaml` | Shared YAML config (single source of truth) |
| `02_analysis/helpers/normalize_gsea.R` | Project-local normalization helper (used by master table script) |
| `02_analysis/1.1.core_pipeline.R` | Core pipeline that executes GSEA across all databases |
| `02_analysis/1.5.create_master_tables.R` | Aggregates GSEA results into master CSV tables |

---

## 1. MSigDB Gene Set Loading via msigdbr

### 1.1 The msigdbr Package

Gene sets are fetched at runtime from the Molecular Signatures Database (MSigDB) using the `msigdbr` R package. This package provides pre-compiled ortholog-mapped gene sets for multiple species, eliminating the need to download GMT files manually.

### 1.2 Species Configuration

The species is configured at the project level in `config.R`:

```r
# From 02_analysis/config/config.R
SPECIES <- "Mus musculus"
SPECIES_DB <- "MM"  # msigdbr database species code
```

`msigdbr` accepts the full species name (e.g., `"Mus musculus"`) for its `species` parameter. The `db_species` code (`"MM"`, `"HS"`) is used by the newer v8+ API.

### 1.3 Collections and Subcollections

MSigDB organizes gene sets into collections (major categories) and subcollections. The project configures which databases to analyze in `config.R`:

```r
# From 02_analysis/config/config.R
MSIGDB_DATABASES <- list(
  H = c("H", ""),                           # Hallmark gene sets
  C2_KEGG = c("C2", "CP:KEGG"),            # KEGG pathways
  C2_REACTOME = c("C2", "CP:REACTOME"),    # Reactome pathways
  C2_WIKIPATHWAYS = c("C2", "CP:WIKIPATHWAYS"), # WikiPathways
  C3_TF = c("C3", "TFT:GTRD"),             # TF targets (GTRD)
  C5_BP = c("C5", "GO:BP"),                # GO Biological Process
  C5_MF = c("C5", "GO:MF"),                # GO Molecular Function
  C5_CC = c("C5", "GO:CC")                 # GO Cellular Component
)
```

Each entry is a two-element vector: `c(category, subcategory)`. When the subcategory is `""`, it is treated as `NULL` (fetch the entire category). The core pipeline parses this as follows:

```r
# From 02_analysis/1.1.core_pipeline.R, lines 403-406
db_params <- MSIGDB_DATABASES[[db_name]]
category <- db_params[1]
subcategory <- if (nzchar(db_params[2])) db_params[2] else NULL
```

### 1.4 MSigDB Collection Reference

| Config Key | Category | Subcategory | Description | Typical Size |
|------------|----------|-------------|-------------|--------------|
| `H` | H | (none) | Hallmark gene sets - 50 well-defined biological states/processes | ~50 sets |
| `C2_KEGG` | C2 | CP:KEGG | KEGG pathway maps | ~180 sets |
| `C2_REACTOME` | C2 | CP:REACTOME | Reactome signaling and metabolic pathways | ~1600 sets |
| `C2_WIKIPATHWAYS` | C2 | CP:WIKIPATHWAYS | Community-curated WikiPathways | ~600 sets |
| `C3_TF` | C3 | TFT:GTRD | Transcription factor targets from GTRD | ~500 sets |
| `C5_BP` | C5 | GO:BP | Gene Ontology Biological Process | ~7500 sets |
| `C5_MF` | C5 | GO:MF | Gene Ontology Molecular Function | ~1700 sets |
| `C5_CC` | C5 | GO:CC | Gene Ontology Cellular Component | ~1000 sets |

### 1.5 TERM2GENE Data Frame Format

`clusterProfiler::GSEA()` requires gene sets in a two-column `TERM2GENE` data frame mapping gene set names to individual genes. The `run_gsea()` function constructs this from the `msigdbr` output:

```r
# From run_gsea.R, line 172
term2gene_df <- msigdb_df[, c("gs_name", "gene_symbol")]
```

The resulting data frame looks like:

| gs_name | gene_symbol |
|---------|-------------|
| HALLMARK_MYC_TARGETS_V1 | Gnl3 |
| HALLMARK_MYC_TARGETS_V1 | Ddx21 |
| HALLMARK_MYC_TARGETS_V1 | Kpna2 |
| HALLMARK_E2F_TARGETS | Brms1l |
| ... | ... |

This is a long-format table where each row represents one gene's membership in one gene set. Gene symbols (not Ensembl IDs) are used because the DE results have gene symbols as rownames.

---

## 2. Handling msigdbr v7.5 vs v8+ API Differences

The `msigdbr` package underwent a breaking API change between versions 7.5.x and 8+:

| | v7.5.x (Legacy) | v8+ (New) |
|---|---|---|
| **Category parameter** | `category` | `collection` |
| **Subcategory parameter** | `subcategory` | `subcollection` |
| **Species parameter** | `species` (full name) | `species` + `db_species` |

The `run_gsea()` function in the toolkit detects the installed version at runtime by inspecting the function signature:

```r
# From run_gsea.R, lines 109-111
msigdbr_params <- names(formals(msigdbr::msigdbr))
use_new_api <- "collection" %in% msigdbr_params
```

Then it dispatches accordingly:

```r
# From run_gsea.R, lines 112-142
if (use_new_api) {
  # New API (v8+)
  if (nzchar(subcollection)) {
    msigdb_df <- msigdbr(
      db_species    = db_species,
      species       = species,
      collection    = collection,
      subcollection = subcollection
    )
  } else {
    msigdb_df <- msigdbr(
      db_species    = db_species,
      species       = species,
      collection    = collection
    )
  }
} else {
  # Legacy API (v7.5.x): category/subcategory
  if (nzchar(subcollection)) {
    msigdb_df <- msigdbr(
      species       = species,
      category      = collection,
      subcategory   = subcollection
    )
  } else {
    msigdb_df <- msigdbr(
      species       = species,
      category      = collection
    )
  }
}
```

This dual-API support means the toolkit works regardless of which `msigdbr` version is installed. The `run_gsea()` function always accepts `collection`/`subcollection` as its own parameter names, and translates to the installed API internally.

**Note on the core pipeline:** The `1.1.core_pipeline.R` script calls `run_gsea()` using parameter names `category` and `subcategory` (lines 413-414), which are passed as extra arguments. This works because `run_gsea()` also accepts positional matching, but the canonical parameter names are `collection` and `subcollection`. Future code should use the canonical names.

---

## 3. Gene Ranking

### 3.1 Rank Metric: t-statistic

Genes are ranked by the moderated t-statistic from the limma-voom differential expression analysis. This is the default and recommended metric:

```r
# From 02_analysis/config/config.R
RANK_METRIC <- "t"  # t-statistic from limma
```

The t-statistic is arbitrarily preferred over logFC. t-statistic incorporates both the magnitude of the fold change and the precision of the estimate (standard error). This is expected to give more weight to genes with consistent effects across replicates.

### 3.2 Ranked Gene Vector Construction

Inside `run_gsea()`, the ranked gene list is built from the DE results table:

```r
# From run_gsea.R, lines 80-91
gene_vector <- DE_results[[rank_metric]]
# Check for NAs in ranking metric
if (any(is.na(gene_vector))) {
    warning("NA values found in rank_metric column '", rank_metric, "'. Removing corresponding genes.")
    valid_indices <- !is.na(gene_vector)
    gene_vector <- gene_vector[valid_indices]
    gene_names <- rownames(DE_results)[valid_indices]
} else {
    gene_names <- rownames(DE_results)
}
names(gene_vector) <- gene_names
ranked_genes <- sort(gene_vector, decreasing = TRUE)
```

The resulting `ranked_genes` is a named numeric vector sorted in descending order (highest t-statistic first). Genes with positive t-statistics (upregulated in the numerator group of the contrast) appear at the top; downregulated genes appear at the bottom.

Infinite values are also removed if present:

```r
# From run_gsea.R, lines 94-97
if (any(is.infinite(ranked_genes))) {
    warning("Infinite values found in ranked gene list. Removing them.")
    ranked_genes <- ranked_genes[!is.infinite(ranked_genes)]
}
```

### 3.3 DE Results Table Structure

The DE results table comes from `limma::topTable()` with gene symbols as rownames:

```r
# From 1.1.core_pipeline.R, lines 338-340
results <- topTable(fit_results$fit,
                   coef = contrast_name,
                   number = Inf,
                   sort.by = "none")
```

Key columns used by the GSEA pipeline:

| Column | Description |
|--------|-------------|
| `t` | Moderated t-statistic (default rank metric) |
| `logFC` | Log2 fold change |
| `P.Value` | Raw p-value |
| `adj.P.Val` | BH-adjusted p-value (FDR) |
| `AveExpr` | Average log2 expression |
| `B` | Log-odds of differential expression |

---

## 4. GSEA Execution

### 4.1 clusterProfiler::GSEA with fgsea Backend

The actual GSEA computation is performed by `clusterProfiler::GSEA()`, which internally delegates to the `fgsea` package. The call in `run_gsea()`:

```r
# From run_gsea.R, lines 175-185
set.seed(seed)
GSEA_result <- clusterProfiler::GSEA(
    geneList = ranked_genes,
    TERM2GENE = term2gene_df,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = padj_method,
    eps = 0,              # Recommended setting for fgsea
    by = "fgsea",         # Use fgsea implementation
    nPermSimple = nperm,  # Pass permutation number to fgsea
    verbose = FALSE       # Keep GSEA quiet
)
```

### 4.2 Parameter Details

| Parameter | Value | Explanation |
|-----------|-------|-------------|
| `eps` | `0` | Boundary for p-value estimation. Setting to 0 enables exact p-value calculation via the adaptive multi-level splitting algorithm in fgsea, rather than truncating at a minimum. This is critical for very significant pathways where p-values can be extremely small. |
| `nPermSimple` | `100000` | Number of permutations for the simple (non-adaptive) fgsea method. Higher values improve p-value precision at the cost of computation time. 100,000 is sufficient for most analyses. |
| `pvalueCutoff` | `1.0` | Set to 1.0 to retain ALL pathways in the result, not just significant ones. This is essential because filtering happens downstream (at the master table or visualization stage), and storing all results allows re-analysis at different thresholds. |
| `pAdjustMethod` | `"fdr"` | Benjamini-Hochberg false discovery rate correction for multiple testing. |
| `by` | `"fgsea"` | Use the fgsea algorithm (fast preranked GSEA). This is the only supported backend in current clusterProfiler versions. |
| `seed` | `123` | Random seed for reproducibility. Set via `set.seed()` before the GSEA call. |

### 4.3 The gseaResult Object

`clusterProfiler::GSEA()` returns a `gseaResult` S4 object. The key slot is `@result`, a data frame with these columns:

| Column | Description |
|--------|-------------|
| `ID` | Gene set identifier (e.g., `HALLMARK_MYC_TARGETS_V1`) |
| `Description` | Gene set description (same as ID for MSigDB) |
| `setSize` | Number of genes in the gene set that were found in the ranked list |
| `enrichmentScore` | Raw enrichment score (ES) |
| `NES` | Normalized enrichment score (accounts for gene set size) |
| `pvalue` | Nominal p-value from the permutation test |
| `p.adjust` | BH-adjusted p-value |
| `qvalue` | q-value (alternative FDR estimate) |
| `rank` | Position in the ranked list where the running sum reaches its maximum |
| `leading_edge` | String describing leading edge statistics |
| `core_enrichment` | Slash-separated list of leading edge genes |

### 4.4 How the Pipeline Executes Multiple Databases

The core pipeline (`1.1.core_pipeline.R`) runs GSEA in two batches, each cached as a separate checkpoint:

**Batch 1 (S1.5):** Hallmark + C2 + C3 collections
```r
# From 1.1.core_pipeline.R, lines 388-438
gsea_h_c2 <- load_or_compute(
  checkpoint_file = CHECKPOINT_GSEA_H_C2,   # "1.1_gsea_H_C2.rds"
  description = "GSEA results (H, C2)",
  compute_fn = function() {
    # Runs: H, C2_KEGG, C2_REACTOME, C2_WIKIPATHWAYS, C3_TF
    db_names <- c("H", "C2_KEGG", "C2_REACTOME", "C2_WIKIPATHWAYS", "C3_TF")
    # ... loops through each database calling run_gsea() ...
  }
)
```

**Batch 2 (S1.6):** GO terms (C5)
```r
# From 1.1.core_pipeline.R, lines 446-496
gsea_c5 <- load_or_compute(
  checkpoint_file = CHECKPOINT_GSEA_C5,      # "1.1_gsea_C5.rds"
  description = "GSEA results (C5 GO terms)",
  compute_fn = function() {
    # Runs: C5_BP, C5_MF, C5_CC
    db_names <- c("C5_BP", "C5_MF", "C5_CC")
    # ... loops through each database calling run_gsea() ...
  }
)
```

**Combination (S1.7):** All results merged into a single named list:
```r
# From 1.1.core_pipeline.R, lines 504-523
all_gsea <- load_or_compute(
  checkpoint_file = CHECKPOINT_ALL_GSEA,     # "1.1_all_gsea_results.rds"
  description = "All GSEA results combined",
  compute_fn = function() {
    combined <- c(gsea_h_c2, gsea_c5)
    return(combined)
  }
)
```

The result is a named list where each element is a `gseaResult` object:
```
all_gsea$H           → gseaResult (Hallmark)
all_gsea$C2_KEGG     → gseaResult (KEGG)
all_gsea$C2_REACTOME → gseaResult (Reactome)
all_gsea$C2_WIKIPATHWAYS → gseaResult (WikiPathways)
all_gsea$C3_TF       → gseaResult (TF targets)
all_gsea$C5_BP       → gseaResult (GO BP)
all_gsea$C5_MF       → gseaResult (GO MF)
all_gsea$C5_CC       → gseaResult (GO CC)
```

---

## 5. Result Normalization

### 5.1 The Normalize-Then-Visualize Pattern

Raw `gseaResult` objects are not suitable for direct aggregation across databases. The normalization step converts each `gseaResult` into a standardized tibble with a consistent schema. This is the bridge between the R compute layer and all downstream visualization:

```
gseaResult objects (per database)
    ↓ normalize_gsea_results()
Standardized tibbles (per database)
    ↓ bind_rows()
Master GSEA table (CSV)
    ↓
All visualizations read from CSV
```

### 5.2 normalize_gsea_results() Function

There are two versions of this function:

1. **Toolkit version** (`01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R`) -- the canonical, reusable implementation with `format_pathway_name()` integration and NES column capitalized as `NES`.
2. **Project-local version** (`02_analysis/helpers/normalize_gsea.R`) -- a slightly adapted version used by the master table script, with the NES column lowercase as `nes` and a project-specific `clean_pathway_name()` / `smart_capitalize()` function.

The project-local version is what actually produces the master table. Its signature:

```r
# From 02_analysis/helpers/normalize_gsea.R
normalize_gsea_results <- function(
    gsea_obj,
    database,
    contrast,
    padj_cutoff = 1,
    strip_prefix = TRUE,
    max_name_length = 60
)
```

### 5.3 Normalization Steps

The function performs these transformations:

1. **Extract result slot** from the `gseaResult` S4 object:
   ```r
   result_df <- as.data.frame(gsea_obj@result)
   ```

2. **Detect adjusted p-value column** (handles different GSEA implementations):
   ```r
   padj_col <- if ("qvalue" %in% colnames(result_df)) "qvalue" else "p.adjust"
   ```

3. **Calculate leading edge size** from `core_enrichment` (slash-separated gene list):
   ```r
   result_df$leading_edge_size <- sapply(result_df$core_enrichment, function(x) {
     if (is.na(x) || x == "") return(0L)
     length(strsplit(as.character(x), "/")[[1]])
   })
   ```

4. **Calculate gene ratio** (leading edge proportion):
   ```r
   result_df$gene_ratio <- result_df$leading_edge_size / result_df$setSize
   ```

5. **Clean pathway names** (strip database prefixes, title-case, smart capitalization):
   ```r
   result_df$pathway_name <- clean_pathway_name(
     result_df$Description,
     strip_prefix = strip_prefix,
     max_length = max_name_length
   )
   ```
   Prefixes removed: `HALLMARK_`, `KEGG_`, `REACTOME_`, `WP_`, `GOBP_`, `GOCC_`, `GOMF_`, `GO_`, `MITOPATHWAYS_`, `MITOXPLORER_`, `MITOCARTA_`.

   Smart capitalization preserves biological abbreviations (ATP, NADH, OXPHOS, IL, TNF, STAT, etc.) while title-casing other words.

6. **Build standardized tibble** with computed columns:
   ```r
   normalized <- tibble::tibble(
     pathway_id          = as.character(result_df$ID),
     pathway_name        = as.character(result_df$pathway_name),
     nes                 = as.numeric(result_df$NES),
     pvalue              = as.numeric(result_df$pvalue),
     padj                = as.numeric(result_df[[padj_col]]),
     set_size            = as.integer(result_df$setSize),
     leading_edge_size   = as.integer(result_df$leading_edge_size),
     gene_ratio          = as.numeric(result_df$gene_ratio),
     core_enrichment     = as.character(result_df$core_enrichment),
     database            = database,
     contrast            = contrast
   ) %>% dplyr::mutate(
     neg_log_padj = -log10(padj),
     neg_log_padj = ifelse(is.infinite(neg_log_padj), 16, neg_log_padj),
     direction    = ifelse(nes > 0, "Up", "Down")
   )
   ```

### 5.4 Database Display Name Mapping

During master table creation, the internal config keys are mapped to human-readable names:

```r
# From 02_analysis/1.5.create_master_tables.R, lines 76-85
DB_DISPLAY_NAMES <- c(
  H = "Hallmark",
  C2_KEGG = "KEGG",
  C2_REACTOME = "Reactome",
  C2_WIKIPATHWAYS = "WikiPathways",
  C5_BP = "GO_BP",
  C5_MF = "GO_MF",
  C5_CC = "GO_CC",
  Mitochondria = "Mitochondria"
)
```

---

## 6. Checkpoint Caching

### 6.1 The load_or_compute() Pattern

Every expensive GSEA computation is wrapped in `load_or_compute()`, which provides transparent checkpoint caching:

```r
# From 02_analysis/config/config.R
load_or_compute <- function(checkpoint_file,
                            compute_fn,
                            force_recompute = FALSE,
                            description = "Result") {
  # Resolve path relative to DIR_CHECKPOINTS
  if (!grepl("^/", checkpoint_file) && !grepl("^[A-Z]:", checkpoint_file)) {
    checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)
  } else {
    checkpoint_path <- checkpoint_file
  }

  # Load from cache if available
  if (file.exists(checkpoint_path) && !force_recompute) {
    message(sprintf("[CACHE] Loading %s from: %s", description, checkpoint_file))
    result <- readRDS(checkpoint_path)
    return(result)
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

### 6.2 GSEA Checkpoint Files

All checkpoints are stored as RDS files in `03_results/checkpoints/`:

| Checkpoint File | Contents | Dependencies |
|-----------------|----------|--------------|
| `1.1_gsea_H_C2.rds` | Named list of gseaResult objects for H, C2_KEGG, C2_REACTOME, C2_WIKIPATHWAYS, C3_TF | DE results |
| `1.1_gsea_C5.rds` | Named list of gseaResult objects for C5_BP, C5_MF, C5_CC | DE results |
| `1.1_all_gsea_results.rds` | Combined list (H_C2 + C5) | Both above |
| `1.3_gsea_mito.rds` | Mitochondrial GSEA (legacy) | DE results + custom gene sets |
| `1.4_gsea_mito_unified.rds` | Unified mitochondrial GSEA | DE results + unified mito gene sets |
| `1.4_gsea_mito_mitopathways.rds` | MitoPathways GSEA | DE results + MitoPathways gene sets |
| `1.4_gsea_mito_mitoxplorer.rds` | MitoXplorer GSEA | DE results + mitoXplorer gene sets |
| `1.4_all_mito_gsea.rds` | All mitochondrial results combined | Above three |
| `1.9_gsea_transportdb.rds` | TransportDB family GSEA | DE results + TransportDB gene sets |

### 6.3 Force Recomputation

To force recomputation (e.g., after changing parameters), set `force_recompute = TRUE`:

```r
gsea_h_c2 <- load_or_compute(
  checkpoint_file = CHECKPOINT_GSEA_H_C2,
  description = "GSEA results (H, C2)",
  force_recompute = TRUE,  # <-- bypass cache
  compute_fn = function() { ... }
)
```

---

## 7. Master Table Creation

### 7.1 Pipeline Flow

The `1.5.create_master_tables.R` script reads cached checkpoint files and produces CSV master tables:

```
03_results/checkpoints/1.1_all_gsea_results.rds   ─┐
03_results/checkpoints/1.4_gsea_mito_unified.rds   ─┤
                                                     ├──→  normalize_gsea_results()  ──→  bind_rows()
                                                     │         (per database)
                                                     └──→  03_results/tables/master_gsea_table.csv
```

The normalization loop:

```r
# From 02_analysis/1.5.create_master_tables.R, lines 164-187
for (db_name in names(all_gsea)) {
  gsea_obj <- all_gsea[[db_name]]
  db_display <- DB_DISPLAY_NAMES[db_name]

  normalized <- normalize_gsea_results(
    gsea_obj = gsea_obj,
    database = db_display,
    contrast = contrast_name,
    padj_cutoff = 1,          # Keep ALL pathways
    strip_prefix = TRUE
  )

  all_normalized[[db_name]] <- normalized
}

master_gsea <- bind_rows(all_normalized)
write_csv(master_gsea, file.path(DIR_TABLES, "master_gsea_table.csv"))
```

### 7.2 Master GSEA Table Schema

The `master_gsea_table.csv` has the following columns:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `pathway_id` | character | Original MSigDB gene set ID | `HALLMARK_MYC_TARGETS_V1` |
| `pathway_name` | character | Cleaned, human-readable name (prefix stripped, title-cased) | `Myc Targets V1` |
| `nes` | numeric | Normalized Enrichment Score. Positive = enriched in numerator group (e.g., IL2RA-KO); negative = enriched in denominator (e.g., NTC). | `3.16` |
| `pvalue` | numeric | Nominal p-value from permutation test | `3.44e-32` |
| `padj` | numeric | BH-adjusted p-value (FDR) | `3.44e-31` |
| `set_size` | integer | Number of genes from the gene set found in the ranked list | `194` |
| `leading_edge_size` | integer | Number of genes in the leading edge (core enrichment) | `140` |
| `gene_ratio` | numeric | `leading_edge_size / set_size` | `0.722` |
| `core_enrichment` | character | Slash-separated list of leading edge gene symbols | `Gnl3/Ddx21/Kpna2/...` |
| `database` | character | Source database display name | `Hallmark`, `KEGG`, `Reactome`, `GO_BP`, etc. |
| `contrast` | character | Name of the DE contrast | `IL2RAKO_vs_NTC` |
| `neg_log_padj` | numeric | `-log10(padj)`, capped at 16 for infinite values | `30.46` |
| `direction` | character | `"Up"` if NES > 0, `"Down"` if NES < 0 | `Up` |
| `module_theme` | character | (Optional) Thematic grouping for GATOM modules | `NA` for MSigDB |
| `n_metabolites` | integer | (Optional) GATOM-specific field | `NA` for MSigDB |
| `n_enzymes` | integer | (Optional) GATOM-specific field | `NA` for MSigDB |

The last three columns (`module_theme`, `n_metabolites`, `n_enzymes`) appear in the CSV because other analysis pipelines (e.g., GATOM) append results with additional columns. For standard MSigDB GSEA rows these are `NA`.

### 7.3 Output Files

| File | Description |
|------|-------------|
| `master_gsea_table.csv` | All GSEA results across all databases (no filtering) |
| `master_gsea_significant.csv` | Filtered to FDR < 0.05, sorted by database then padj |
| `gsea_summary_stats.csv` | Per-database/direction summary counts and mean NES |

### 7.4 Schema Validation

The `pipeline.yaml` defines the expected schema for validation:

```yaml
# From 02_analysis/config/pipeline.yaml
schemas:
  master_gsea_table:
    version: "1.0.0"
    required_columns:
      - pathway_id
      - pathway_name
      - database
      - nes
      - pvalue
      - padj
      - core_enrichment
    optional_columns:
      - contrast
      - direction
      - set_size
      - leading_edge_size
```

---

## 8. Custom (Non-MSigDB) Gene Sets

The pipeline also supports GSEA with gene sets from external databases. These use the same `clusterProfiler::GSEA()` engine but with custom TERM2GENE data frames instead of msigdbr.

### 8.1 External Gene Set Parsers

The toolkit provides parsers in `parse_external_genesets.R`:

| Function | Input Format | Database |
|----------|-------------|----------|
| `parse_transportdb()` | CSV with Family/Symbol columns | TransportDB 2.0 |
| `parse_gmt()` | GMT (Gene Matrix Transposed) | Any GMT file |
| `parse_gmx()` | GMX (column-oriented) | MitoPathways 3.0 |
| `parse_mitoxplorer()` | TSV with mito_process/MGI_symbol | mitoXplorer 3.0 |
| `parse_geneset_file()` | Generic TSV/CSV | Any tabular gene set file |

All parsers return a list with `T2G` (TERM2GENE data frame) and `T2N` (TERM2NAME data frame), which can be passed directly to `clusterProfiler::GSEA()`:

```r
transport_db <- parse_transportdb("00_data/references/TransportDB2.0.csv")
gsea_result <- clusterProfiler::GSEA(
    geneList = ranked_genes,
    TERM2GENE = transport_db$T2G,
    pvalueCutoff = 1,
    eps = 0,
    nPermSimple = 100000
)
```

### 8.2 Utility Functions

`pathway_utils.R` provides additional helpers:

- `filter_pathways_by_size(T2G, min_size, max_size)` -- Remove gene sets outside size bounds
- `list_to_term2gene(named_list)` -- Convert `list(pathway = c(genes))` to TERM2GENE format
- `term2gene_to_list(T2G)` -- Reverse: TERM2GENE to named list
- `convert_geneset_ids(T2G, org_db, from_type, to_type)` -- Convert gene IDs (e.g., Ensembl to Symbol)
- `convert_human_to_mouse(T2G)` -- Map human gene symbols to mouse orthologs via homologene
- `export_to_gmx(T2G, T2N, output_file)` -- Export to GMX format for external tools

---

## 9. The run_gsea_analysis() Multi-Database Pipeline

For standalone use (outside the 12868-EH project pipeline), the toolkit provides `run_gsea_analysis()`, which runs GSEA across multiple databases and auto-generates plots:

```r
# From run_gsea_analysis.R
results <- run_gsea_analysis(
    de_table = de_results,
    analysis_name = "MyAnalysis",
    rank_metric = "t",
    species = "Mus musculus",
    n_pathways = 30,
    padj_cutoff = 0.05,
    nperm = 100000,
    pvalue_cutoff = 1,       # Store ALL pathways
    save_plots = TRUE,
    output_dir = "./GSEA_Plots"
)
```

This function:
1. Auto-sources all required helper scripts from the toolkit
2. Defines a default database list (Hallmark, Canonical, GO BP/MF/CC, KEGG, Reactome, WikiPathways, CGP, GTRD)
3. Loops through each database calling `run_gsea()`
4. Generates dotplots (up/down/faceted), barplots, running sum plots, and heatmaps per database
5. Returns a named list of `gseaResult` objects

---

## 10. End-to-End Example

A complete GSEA run for a single database, from DE results to master table row:

```r
# 1. Load DE results (gene symbols as rownames, 't' column present)
de_results <- readRDS("03_results/checkpoints/1.1_de_results.rds")
de_table <- de_results[["IL2RA_KO - NTC"]]

# 2. Run GSEA for Hallmark
gsea_hallmark <- run_gsea(
    DE_results = de_table,
    rank_metric = "t",
    species = "Mus musculus",
    collection = "H",
    subcollection = "",
    pvalue_cutoff = 1,
    nperm = 100000,
    seed = 123
)

# 3. Normalize to standard tibble
source("02_analysis/helpers/normalize_gsea.R")
hallmark_df <- normalize_gsea_results(
    gsea_obj = gsea_hallmark,
    database = "Hallmark",
    contrast = "IL2RAKO_vs_NTC",
    padj_cutoff = 1,
    strip_prefix = TRUE
)

# 4. Combine with other databases and write master table
master_gsea <- bind_rows(hallmark_df, kegg_df, reactome_df, ...)
write_csv(master_gsea, "03_results/tables/master_gsea_table.csv")
```

---

## 11. Key Design Decisions and Rationale

1. **pvalueCutoff = 1.0**: All pathways are retained in the gseaResult object. Filtering happens at the visualization/export stage. This avoids re-running expensive GSEA when exploring different significance thresholds.

2. **eps = 0**: Enables exact p-value calculation via fgsea's adaptive multi-level splitting. Without this, extremely significant pathways get truncated p-values.

3. **Gene symbols (not Ensembl IDs)**: The ranked list uses gene symbols as names because msigdbr returns gene sets mapped to symbols. The Ensembl-to-symbol conversion happens during DGEList construction (S1.1 of the pipeline).

4. **Separate H/C2 and C5 checkpoints**: GO term analysis (C5) is the most computationally expensive batch (thousands of gene sets). Separating it allows the H/C2 checkpoint to be loaded quickly when only those results are needed.

5. **Two normalization functions**: The toolkit version (`normalize_gsea.R`) uses uppercase `NES`; the project-local version uses lowercase `nes`. The master table uses lowercase. This is a known inconsistency that should be unified in future toolkit updates.

6. **neg_log_padj capped at 16**: Very small p-values produce extremely large -log10 values. Capping at 16 prevents visualization scale distortion while still indicating extreme significance.

7. **Prefix stripping in pathway names**: MSigDB gene set IDs include database prefixes (e.g., `HALLMARK_MYC_TARGETS_V1`). These are stripped and the remaining text is title-cased for readability, while the original ID is preserved in the `pathway_id` column for programmatic access.
