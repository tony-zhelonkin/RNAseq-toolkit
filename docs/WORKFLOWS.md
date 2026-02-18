# RNAseq-toolkit Workflow Guide

> Quick-scan reference for plugging toolkit functions into analysis pipelines.
> Prevents slop. Token-efficient.

---

## 1. Sourcing Functions

```r
# === GSEA Processing ===
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/parse_external_genesets.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/pathway_utils.R")

# === GSEA Plotting ===
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_barplot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/format_pathway_names.R")

# === DE Visualization ===
source("{path}/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")
source("{path}/RNAseq-toolkit/scripts/DE/plotPCA.R")

# === Theme (required by plotting functions) ===
source("{path}/RNAseq-toolkit/scripts/custom_minimal_theme.R")
```

---

## 2. Core Function Signatures

### `run_gsea()` - MSigDB GSEA

**Input:** DE table with gene symbols as rownames, must have rank metric column.

```r
gsea_result <- run_gsea(
  DE_results   = de_table,          # data.frame, rownames = gene symbols
  rank_metric  = "t",               # Column to rank by (default: "t")
  species      = "Mus musculus",    # or "Homo sapiens"
  category     = "H",               # MSigDB: H, C2, C3, C5, etc.
  subcategory  = NULL,              # e.g., "CP:KEGG", "BP", "CP:REACTOME"
  pvalue_cutoff = 1,                # Filter at GSEA step (1 = no filter)
  nperm        = 100000,            # Permutations
  seed         = 123                # Reproducibility
)
```

**Returns:** `gseaResult` S4 object from clusterProfiler.

**Category/Subcategory Reference:**

| Database | category | subcategory |
|----------|----------|-------------|
| Hallmark | `"H"` | `NULL` |
| KEGG | `"C2"` | `"CP:KEGG"` |
| Reactome | `"C2"` | `"CP:REACTOME"` |
| WikiPathways | `"C2"` | `"CP:WIKIPATHWAYS"` |
| GO:BP | `"C5"` | `"GO:BP"` |
| GO:MF | `"C5"` | `"GO:MF"` |
| GO:CC | `"C5"` | `"GO:CC"` |
| TF Targets (GTRD) | `"C3"` | `"TFT:GTRD"` |

---

### `normalize_gsea_results()` - Standardize Output

Converts `gseaResult` → tibble for aggregation and cross-database comparison.

```r
df <- normalize_gsea_results(
  gsea_obj        = gsea_result,     # gseaResult or data.frame
  database        = "Hallmark",      # Source label
  contrast        = "A_vs_B",        # Contrast name
  padj_cutoff     = 1,               # Filter (1 = all)
  format_names    = TRUE,            # Clean pathway names
  max_name_length = 80
)
```

**Output Schema:**
```
pathway_id | pathway_name | database | contrast | NES | pvalue | padj |
set_size | leading_edge_size | gene_ratio | core_enrichment | direction | neg_log_padj
```

---

## 3. External Database Pattern (Custom Gene Sets)

For non-MSigDB databases (TransportDB, MitoDB, etc.), use `clusterProfiler::GSEA()` directly with TERM2GENE/TERM2NAME.

### A. Parse External Database

```r
# Option 1: TransportDB parser
transport_db <- parse_transportdb(
  file     = "00_data/references/TransportDB2.0.csv",
  prefix   = "TRANSPORTDB",
  min_size = 5,
  max_size = 500
)

# Option 2: GMT file
gene_sets <- parse_gmt(
  file   = "00_data/references/custom.gmt",
  prefix = "CUSTOM"
)

# Option 3: Generic TSV/CSV
gene_sets <- parse_geneset_file(
  file     = "myfile.tsv",
  gs_col   = 1,      # Gene set name column
  gene_col = 2,      # Gene symbol column
  prefix   = "MYDB"
)
```

**All return:** `list(T2G, T2N, stats, source)`
- `T2G`: `gs_name | gene_symbol`
- `T2N`: `gs_name | description`

### B. Run GSEA on External Database

```r
# Prepare TERM2GENE/TERM2NAME
T2G <- gene_sets$T2G %>% rename(term = gs_name, gene = gene_symbol)
T2N <- gene_sets$T2N %>% rename(term = gs_name, name = description)

# Create ranked gene list from DE results
ranked_genes <- de_table[["t"]]
names(ranked_genes) <- rownames(de_table)
ranked_genes <- sort(ranked_genes[!is.na(ranked_genes)], decreasing = TRUE)

# Filter T2G to genes in ranked list (IMPORTANT!)
T2G_filtered <- T2G %>% filter(gene %in% names(ranked_genes))

# Run GSEA
set.seed(123)
gsea_result <- clusterProfiler::GSEA(
  geneList      = ranked_genes,
  TERM2GENE     = T2G_filtered,
  TERM2NAME     = T2N,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  by            = "fgsea",
  nPermSimple   = 100000,
  verbose       = FALSE
)

# CRITICAL: Fix geneSets slot for running sum plots
T2G_list <- split(T2G_filtered$gene, T2G_filtered$term)
gsea_result@geneSets <- T2G_list
```

### C. Normalize External Results

Same function, just specify the database name:

```r
df <- normalize_gsea_results(
  gsea_obj = gsea_result,
  database = "TransportDB",
  contrast = "Obesity_effect"
)
```

---

## 4. Visualization Functions

### `gsea_dotplot()`

```r
p <- gsea_dotplot(
  gsea_obj,
  filterBy     = "p.adjust",       # "p.adjust", "NES_positive", "NES_negative", "NES"
  sortBy       = "GeneRatio",      # "GeneRatio" or "p.adjust"
  showCategory = 20,               # Max pathways
  padj_cutoff  = 0.05,
  title        = "Hallmark",
  use_gradient = TRUE,             # Continuous NES gradient
  highlight_sig = TRUE,            # Outline FDR < 0.005
  strip_prefix = TRUE              # Remove HALLMARK_, KEGG_, etc.
)
```

### `gsea_barplot()`

```r
p <- gsea_barplot(
  gsea_obj,
  padj_cutoff  = 0.05,
  top_n        = 30,
  title        = "GSEA NES",
  strip_prefix = TRUE
)
```

### `gsea_running_sum_plot()` - IMPORTANT

**The running sum plot requires the `@geneSets` slot to be populated.**

For MSigDB results via `run_gsea()`: Works automatically.
For external databases: **Must manually fix the slot** (see Section 3B).

```r
# By gene set ID (pathway name)
p <- gsea_running_sum_plot(
  gsea_obj,
  gene_set_ids = c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_ADIPOGENESIS"),
  palette      = c("#E41A1C", "#377EB8"),  # Colors per pathway
  labels       = NULL,                      # Auto from Description
  max_name_length = 45
)

# By row index (top 5 by |NES|)
p <- gsea_running_sum_plot(gsea_obj)  # Defaults to top 5
```

**External Database Gotcha:**

For external DBs, `Description` often equals `ID`. enrichplot uses `Description` as the color aesthetic. Two options:

```r
# Option A: Copy object, set Description = ID
gsea_plot <- gsea_result
gsea_plot@result$Description <- gsea_plot@result$ID
p <- gsea_running_sum_plot(gsea_plot, gene_set_ids = top_ids)

# Option B: Pass explicit labels (recommended)
original_desc <- gsea_result@result$Description
names(original_desc) <- gsea_result@result$ID

p <- gsea_running_sum_plot(
  gsea_result,
  gene_set_ids = top_ids,
  labels = original_desc[top_ids]  # Named vector
)
```

---

## 5. Volcano Plot

```r
p <- create_standard_volcano(
  de_results,                    # rownames = gene IDs
  decision_by   = "fdr",         # "fdr" or "p"
  p_cutoff      = 0.05,          # FDR threshold when decision_by="fdr"
  fc_cutoff     = 2,             # |log2FC| threshold
  top_n         = 5,             # Label top N per side
  highlight_gene = c("Setdb1"),  # Extra genes to label
  title         = "Contrast Name"
)
```

**Y-axis:** -log10(raw P.Value) for resolution
**Decision:** Based on `adj.P.Val` (FDR)
**Dashed line:** Raw p-value corresponding to FDR boundary

---

## 6. Pipeline Integration Checklist

### Phase 1: Compute

```r
# 1. Source required functions
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R")

# 2. Run GSEA per database
for (db in databases) {
  gsea_result <- run_gsea(de_table, category = db$cat, subcategory = db$subcat)
  all_results[[db$name]] <- gsea_result
}

# 3. Normalize to master table
combined <- tibble()
for (db in names(all_results)) {
  df <- normalize_gsea_results(all_results[[db]], database = db, contrast = contrast_name)
  combined <- bind_rows(combined, df)
}

# 4. Save checkpoint + master table
saveRDS(all_results, "checkpoints/gsea_combined.rds")
write_csv(combined, "tables/master_gsea.csv")
```

### Phase 2: Visualize (Read-Only)

```r
# 1. Source plotting functions
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("{path}/RNAseq-toolkit/scripts/custom_minimal_theme.R")

# 2. Load checkpoints (never recompute)
gsea_results <- readRDS("checkpoints/gsea_combined.rds")

# 3. Generate plots
for (contrast in contrasts) {
  for (db in databases) {
    gsea_obj <- gsea_results[[db]][[contrast]]

    # Dotplot
    p <- gsea_dotplot(gsea_obj, padj_cutoff = 0.05)
    ggsave(sprintf("plots/%s_%s_dotplot.pdf", contrast, db), p)

    # Running sum (top 5)
    if (sum(gsea_obj@result$p.adjust < 0.05) > 0) {
      p <- gsea_running_sum_plot(gsea_obj)
      ggsave(sprintf("plots/%s_%s_running_sum.pdf", contrast, db), p)
    }
  }
}
```

---

## 7. Common Pitfalls

| Problem | Cause | Fix |
|---------|-------|-----|
| `run_gsea()` fails | Wrong column name | Ensure `rank_metric` exists in DE table |
| Empty GSEA results | No overlapping genes | Check gene symbol format; filter T2G to ranked genes |
| Running sum plot fails | Empty `@geneSets` slot | For external DBs: `gsea_result@geneSets <- split(T2G$gene, T2G$term)` |
| Pathway names ugly | Not stripped | Ensure `format_pathway_names.R` sourced before plotting |
| Dotplot theme wrong | Missing theme | Source `custom_minimal_theme.R` |
| Colors wrong on running sum | Description mismatch | Pass explicit `labels` parameter for external DBs |

---

## 8. Directory Structure (Toolkit)

```
{path}/RNAseq-toolkit/
├── scripts/
│   ├── GSEA/
│   │   ├── GSEA_processing/
│   │   │   ├── run_gsea.R              # MSigDB GSEA wrapper
│   │   │   ├── normalize_gsea.R        # → Standard tibble
│   │   │   ├── parse_external_genesets.R  # TransportDB, GMT, CSV
│   │   │   └── pathway_utils.R         # Filtering, GMT export
│   │   └── GSEA_plotting/
│   │       ├── gsea_dotplot.R          # Main dotplot
│   │       ├── gsea_barplot.R          # NES barplot
│   │       ├── gsea_running_sum_plot.R # enrichplot wrapper
│   │       └── format_pathway_names.R  # Smart capitalization
│   ├── DE/
│   │   ├── plot_standard_volcano.R     # FDR-aware volcano
│   │   └── plotPCA.R                   # 2D PCA
│   └── custom_minimal_theme.R          # ggplot2 theme
└── CLAUDE.md                           # Full toolkit docs
```

---

## 9. Quick Reference: Minimal Working Example

```r
# === Setup ===
source("02_analysis/config/config.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_processing/normalize_gsea.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")
source("{path}/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("{path}/RNAseq-toolkit/scripts/custom_minimal_theme.R")

library(clusterProfiler)
library(msigdbr)

# === Run Hallmark GSEA ===
gsea_H <- run_gsea(de_results, category = "H", species = "Mus musculus")

# === Normalize ===
df <- normalize_gsea_results(gsea_H, database = "Hallmark", contrast = "Treatment_vs_Control")

# === Plot ===
p1 <- gsea_dotplot(gsea_H, padj_cutoff = 0.05, title = "Hallmark")
p2 <- gsea_running_sum_plot(gsea_H)  # Top 5 by |NES|

ggsave("hallmark_dotplot.pdf", p1, width = 10, height = 8)
ggsave("hallmark_running_sum.pdf", p2, width = 12, height = 10)
```

---

*Last updated: 2025-12-30*
