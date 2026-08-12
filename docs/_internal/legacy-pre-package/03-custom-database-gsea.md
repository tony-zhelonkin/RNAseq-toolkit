# Custom Database GSEA: From Raw Files to Master Table

This document describes how to integrate external (non-MSigDB) gene set databases
into the clusterProfiler GSEA pipeline and merge results into a unified master
table. It covers the T2G/T2N data model, database-specific parsing patterns,
gene ID mapping strategies, quality control, result normalization, and the
idempotent append pattern used for master table integration.

All code examples are drawn from the 12868-EH project (CAR-T IL2RA-KO analysis,
*Mus musculus*). The same patterns apply to human analyses with appropriate
organism substitutions.

---

## 0. Pre-packaged Reference Databases (Quick Start)

The RNAseq-toolkit ships bundled reference databases that can be loaded with a single function call — no manual file downloading required.

### Available Bundled Databases

| Database ID | Name | Species | Description |
|-------------|------|---------|-------------|
| `mitopathways` | MitoPathways 3.0 | Mouse (from human via homologene) | Hierarchical mitochondrial pathways |
| `mitoxplorer` | mitoXplorer 3.0 | Mouse (native MGI symbols) | Mitochondrial functional classifications |
| `mito_unified` | Unified Mito | Mouse | MitoPathways + mitoXplorer merged with Jaccard deduplication |
| `transportdb` | TransportDB 2.0 | Mouse | Membrane transporter families (TCDB classification) |

GATOM network files (~24MB) are not bundled. Use `download_gatom_references()` instead.

### Usage

```r
# Source the loader
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))

# Load a pre-built database (returns T2G/T2N list)
db <- load_reference_db("mito_unified")

# Use directly with clusterProfiler::GSEA()
gsea_result <- GSEA(ranked_genes, TERM2GENE = db$T2G, TERM2NAME = db$T2N,
                    eps = 0, nPermSimple = 100000, pvalueCutoff = 1.0)

# List all available databases
list_reference_dbs()

# Get citation info
info <- get_reference_db_info("mitopathways")
cat(info$citations_text)
```

See `data/references/README.md` for full details on bundled databases, rebuilding, and citations.

> **Note:** The sections below document how to parse databases from raw files (BYOD workflow). If using the pre-packaged databases above, you can skip directly to [Section 5: Running GSEA](#5-running-gsea-with-custom-databases).

---

## Table of Contents

1. [The T2G/T2N Data Model](#1-the-t2gt2n-data-model)
2. [Database-Specific Parsing Patterns](#2-database-specific-parsing-patterns)
   - [MitoPathways (GMX format, human-to-mouse)](#21-mitopathways-gmx-format-human-to-mouse-conversion)
   - [mitoXplorer (TSV format, native mouse)](#22-mitoxplorer-tsv-format-native-mouse)
   - [TransportDB (CSV with RefSeq protein IDs)](#23-transportdb-csv-with-refseq-protein-ids)
   - [GATOM (igraph module extraction)](#24-gatom-igraph-module-extraction)
3. [Gene ID Mapping Strategies](#3-gene-id-mapping-strategies)
4. [Quality Control and Filtering](#4-quality-control-and-filtering)
5. [Running GSEA with Custom Databases](#5-running-gsea-with-custom-databases)
6. [Result Normalization to Master Table Schema](#6-result-normalization-to-master-table-schema)
7. [Master Table Integration (Idempotent Append)](#7-master-table-integration-idempotent-append)
8. [Recipe: Adding a New Custom Database](#8-recipe-adding-a-new-custom-database)

---

## 1. The T2G/T2N Data Model

`clusterProfiler::GSEA()` accepts custom gene sets through two data frames
passed as the `TERM2GENE` and `TERM2NAME` arguments. Throughout this codebase
we refer to these as **T2G** and **T2N**.

### T2G (TERM2GENE)

A long-format data frame with exactly two columns. Each row maps one gene to
one gene set. A gene set with 50 genes produces 50 rows.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `gs_name` | character | Gene set identifier, prefixed by database (e.g. `MITOPATHWAYS_Oxidative_Phosphorylation`) |
| 2 | `gene_symbol` | character | Gene symbol matching the rownames of the DE results (e.g. `Cox5a`, `Ndufs1`) |

### T2N (TERM2NAME)

A lookup table mapping each gene set identifier to a human-readable description.
One row per gene set.

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `gs_name` | character | Must match exactly the values used in T2G |
| 2 | `description` | character | Human-readable pathway name |

### Naming convention

Gene set identifiers always carry a database prefix in SCREAMING_SNAKE_CASE:

```
{DATABASE}_{pathway_name}

MITOPATHWAYS_Oxidative_Phosphorylation.Complex_I
MITOXPLORER_TRANSLATION_INITIATION
TRANSPORTDB_AAAP
GATOM_GLUTAMINE_GLUTATHIONE_AXIS
```

This prefix serves two purposes: (1) it prevents name collisions when multiple
databases define identically-named pathways, and (2) it allows downstream code
to determine the source database from the pathway ID alone.

### Storage format

Parsed databases are persisted as RDS lists with the following structure:

```r
db <- list(
  T2G      = data.frame(gs_name, gene_symbol),
  T2N      = data.frame(gs_name, description),
  source   = "MitoPathways3.0",
  created  = Sys.time()
)
saveRDS(db, "00_data/references/mitochonria/mitocarta3.0/mito_mitopathways.rds")
```

Some databases also export a GMX file for use with external tools like the
Broad GSEA desktop application.

---

## 2. Database-Specific Parsing Patterns

### 2.1 MitoPathways (GMX format, human-to-mouse conversion)

**Source file:** `00_data/references/mitochonria/mitocarta3.0/MitoPathways3.0.gmx`

**Format:** GMX (Gene Matrix eXtended) -- the *transpose* of GMT:

```
Row 1:  Short names    (one per column)
Row 2:  Full names     (hierarchical: pathway.subpathway.detail)
Row 3+: Gene symbols   (one gene per row per column, variable length)
```

**Challenge:** MitoPathways ships human gene symbols. Mouse analysis requires
ortholog conversion.

**Parser:** `parse_mitopathways_gmx()` in `02_analysis/helpers/parse_mito_databases.R`

```r
mp_data <- parse_mitopathways_gmx(
  gmx_file          = "00_data/references/mitochonria/mitocarta3.0/MitoPathways3.0.gmx",
  prefix            = "MITOPATHWAYS",
  convert_to_mouse  = TRUE,
  org_db            = org.Mm.eg.db
)
```

The parser reads the GMX column-by-column, extracts genes from rows 3+,
optionally converts human symbols to mouse orthologs (see Section 3), and
emits T2G/T2N. The full hierarchical name from row 2 becomes the `gs_name`
(after prefixing), and the short name from row 1 becomes the `description`.

The toolkit also provides a generic `parse_gmx()` in
`scripts/GSEA/GSEA_processing/parse_external_genesets.R` that handles any
GMX file without species conversion.

### 2.2 mitoXplorer (TSV format, native mouse)

**Source file:** `00_data/references/mitochonria/mitoXplorer3.0/raw/mouse_gene_function.txt`

**Format:** Tab-separated with columns including `MGI_symbol` (gene) and
`mito_process` (pathway category).

**No species conversion needed** -- this database ships mouse symbols directly.

**Parser:** `parse_mitoxplorer()` in both the project helpers and the toolkit.

```r
mx_data <- parse_mitoxplorer(
  gene_function_file = "00_data/references/mitochonria/mitoXplorer3.0/raw/mouse_gene_function.txt",
  prefix             = "MITOXPLORER"
)
```

The parser groups genes by the `mito_process` column. Each unique process
becomes a gene set. Special characters in the process name are replaced with
underscores when constructing the `gs_name`:

```r
gs_name <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9_]", "_", proc)))
```

### 2.3 TransportDB (CSV with RefSeq protein IDs)

**Source file:** `00_data/references/TransportDB2.0.csv`

**Format:** Headerless CSV with 7 columns:

| Column | Content | Example |
|--------|---------|---------|
| 1 | RefSeq Protein ID | `NP_001009950` |
| 2 | Substrate description | `glucose` |
| 3 | (empty) | |
| 4 | Family abbreviation | `AAAP` |
| 5 | Family full name | `Amino Acid/Auxin Permease` |
| 6 | Transporter type | `Secondary Transporter` |
| 7 | TCDB family code | `2.A.18` |

**Challenge:** TransportDB uses RefSeq protein accessions (`NP_*`), not gene
symbols. These must be mapped to gene symbols through `org.Mm.eg.db`.

**Parser:** `parse_transportdb()` in `02_analysis/helpers/parse_transportdb.R`

```r
transport_db <- parse_transportdb(
  transportdb_file = "00_data/references/TransportDB2.0.csv",
  prefix           = "TRANSPORTDB",
  org_db           = org.Mm.eg.db,
  group_by         = "family"    # or "substrate"
)
```

The ID conversion is handled by `convert_refseq_to_symbols()`, which:

1. Strips version suffixes from protein IDs (`NP_001009950.1` -> `NP_001009950`).
2. Queries `org.Mm.eg.db` via `AnnotationDbi::select()` using `REFSEQ` keytype.
3. Falls back to a two-step Entrez intermediate if direct REFSEQ lookup fails.

```r
# Primary attempt
mapping <- AnnotationDbi::select(
  org_db,
  keys    = clean_ids,
  columns = "SYMBOL",
  keytype = "REFSEQ"
)

# Fallback: REFSEQ -> ENTREZID -> SYMBOL
entrez_map <- AnnotationDbi::select(org_db, keys = clean_ids,
                                     columns = "ENTREZID", keytype = "REFSEQ")
symbol_map <- AnnotationDbi::select(org_db, keys = entrez_ids,
                                     columns = c("SYMBOL", "REFSEQ"),
                                     keytype = "ENTREZID")
```

Gene sets are formed by grouping mapped gene symbols by the `family_abbrev`
column (column 4). An alternative `group_by = "substrate"` mode is available
but produces noisier, less interpretable gene sets.

### 2.4 GATOM (igraph module extraction)

**Source files:** `03_results/checkpoints/1.6_gatom_kegg.rds`, `1.6_gatom_combined.rds`

GATOM (Gene Active Topology-Oriented Module) analysis produces igraph objects
representing active metabolic subnetworks. These are *not* traditional gene
set databases -- they are topology-aware graphs where reaction order matters.

**The conversion to gene sets is lossy by design.** It exists solely to
inject GATOM modules into the pathway explorer's UMAP for cross-validation.
The script `1.10.gatom_to_gsea.R` documents two explicit caveats:

1. **Topology flattening:** Converting a graph to a gene list strips structural
   insight about reaction ordering.
2. **Metabolite invisibility:** The gene list contains only enzymes; metabolites
   that drove module selection are absent from downstream clustering.

**Extraction pattern:**

```r
module <- gatom_result$module  # igraph object

# Metabolites are nodes
node_df <- igraph::as_data_frame(module, "vertices")

# Enzymes are edges (reactions connecting metabolites)
edge_df <- igraph::as_data_frame(module, "edges")

# Extract gene symbols from edges
gene_col <- NULL
for (col in c("Symbol", "symbol", "gene", "gene_symbol", "label")) {
  if (col %in% colnames(edge_df)) {
    gene_col <- col
    break
  }
}
module_genes <- unique(na.omit(edge_df[[gene_col]]))
```

**Biological naming:** Rather than assigning generic names like "Active Module",
the `derive_module_name()` function analyzes metabolite composition and enzyme
membership to generate informative names:

```r
# Metabolite signature matching
tca_metabolites <- c("pyruvate", "oxaloacetate", "2-oxoglutarate", ...)
gsh_metabolites <- c("glutathione", "glutathione disulfide", ...)

# Enzyme signature matching
gsh_enzymes <- c("Gclc", "Gss", "Gpx4", ...)
tca_enzymes <- c("Idh1", "Idh2", "Sdha", ...)

# Composite naming
# -> "Glutamine-Dependent Redox & Anaplerosis Module (GATOM Combined)"
# -> pathway_id: GATOM_GLUTAMINE_DEPENDENT_REDOX___ANAPLEROSIS
```

**Pseudo-NES calculation:** Since GATOM modules do not go through GSEA, a
pseudo-NES is computed from the average log2FC of module edges, and a
pseudo-p-value from the geometric mean of edge p-values.

---

## 3. Gene ID Mapping Strategies

Three strategies are used depending on the source database:

### Strategy 1: Homologene (cross-species ortholog conversion)

Used for: **MitoPathways** (Human -> Mouse)

```r
# Primary: homologene package
h_res <- homologene::homologene(human_genes, inTax = 9606, outTax = 10090)
mouse_map <- setNames(h_res$`10090`, h_res$`9606`)

# Fallback: Title-case capitalization + org.Mm.eg.db validation
manual_mouse <- tools::toTitleCase(tolower(unmapped_genes))
valid_keys <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
manual_mouse <- manual_mouse[manual_mouse %in% valid_keys]
```

The hybrid approach (homologene first, then validated manual fallback) typically
achieves 85-95% mapping rates. The manual fallback catches genes like `MT-CO1`
that may be absent from the homologene database but have predictable
capitalization patterns between species.

### Strategy 2: AnnotationDbi (same-species ID type conversion)

Used for: **TransportDB** (RefSeq Protein -> Symbol)

```r
mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = refseq_ids,
  columns = "SYMBOL",
  keytype = "REFSEQ"
)
```

Also available as a generic utility in the toolkit:

```r
# From parse_external_genesets.R
T2G <- convert_geneset_ids(
  T2G,
  org_db    = org.Mm.eg.db,
  from_type = "ENSEMBL",    # or "REFSEQ", "ENTREZID"
  to_type   = "SYMBOL",
  drop_unmapped = TRUE
)
```

### Strategy 3: Direct (no conversion needed)

Used for: **mitoXplorer** (already provides MGI mouse symbols),
**GATOM** (symbols extracted from igraph edge attributes)

No conversion step is required. Gene symbols are used as-is after basic
cleaning (trimming whitespace, removing NAs/empty strings).

---

## 4. Quality Control and Filtering

### Gene set size filtering

All databases are filtered to retain only gene sets with between 5 and 500
genes. This is enforced by `filter_pathways_by_size()` from
`02_analysis/helpers/pathway_utils.R`:

```r
filter_pathways_by_size <- function(t2g, t2n, min_size = 5, max_size = 500) {
  sizes <- t2g %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop")

  valid_pathways <- sizes$gs_name[sizes$n >= min_size & sizes$n <= max_size]

  list(
    T2G = t2g[t2g$gs_name %in% valid_pathways, ],
    T2N = t2n[t2n$gs_name %in% valid_pathways, ]
  )
}
```

**Rationale:** Gene sets with fewer than 5 genes produce unreliable enrichment
statistics. Gene sets with more than 500 genes are too broad to be informative
(they test "is the transcriptome generally perturbed" rather than any specific
pathway).

### Cross-database deduplication (Jaccard merge)

When combining MitoPathways and mitoXplorer, nearly identical gene sets are
collapsed using Jaccard similarity:

```r
merged_res <- merge_redundant_pathways(
  combined_t2g,
  combined_t2n,
  similarity_cutoff = 0.99,
  preference_order  = c("MITOPATHWAYS")  # keep this name on collision
)
```

The `preference_order` parameter determines which database's name is kept
when two gene sets are merged. This prevents inflated significance from testing
the same biology twice under different names.

### Duplicate description resolution

After merging, multiple gene sets may share the same description string. The
preparation script appends the source prefix to disambiguate:

```r
final_t2n <- final_t2n %>%
  mutate(description = if_else(
    description %in% dup_descs,
    paste0(description, " (", gsub("_.*", "", gs_name), ")"),
    description
  ))
```

### Pre-GSEA filtering to ranked list

Before running GSEA, the T2G is filtered to include only genes present in
the ranked gene list. This avoids warnings from clusterProfiler and ensures
set sizes reflect actual testable genes:

```r
T2G_filtered <- T2G %>%
  dplyr::filter(gene %in% names(ranked_genes))
```

---

## 5. Running GSEA with Custom Databases

The call to `clusterProfiler::GSEA()` is identical for MSigDB and custom
databases. The only difference is that custom databases pass `TERM2GENE`
and `TERM2NAME` instead of relying on msigdbr's built-in gene sets.

```r
# Rename columns to match clusterProfiler expectations
T2G <- mito_db$T2G %>% dplyr::rename(term = gs_name, gene = gene_symbol)
T2N <- mito_db$T2N %>% dplyr::rename(term = gs_name, name = description)

# Filter to genes in ranked list
T2G <- T2G %>% dplyr::filter(gene %in% names(ranked_genes))

# Run GSEA
set.seed(GSEA_SEED)
gsea_result <- clusterProfiler::GSEA(
  geneList     = ranked_genes,        # Named numeric vector, sorted descending
  TERM2GENE    = T2G,
  TERM2NAME    = T2N,
  pvalueCutoff = 1,                   # Keep all results; filter later
  pAdjustMethod = "BH",
  eps          = 0,                   # Exact p-values for leading edge
  by           = "fgsea",
  nPermSimple  = 10000,
  verbose      = FALSE
)
```

**Important fix:** After running GSEA, the `@geneSets` slot of the result
object may be empty. This must be populated manually for downstream enrichplot
functions (e.g., running sum plots) to work:

```r
T2G_list <- split(T2G$gene, T2G$term)
gsea_result@geneSets <- T2G_list
```

**Checkpoint caching:** All GSEA runs are wrapped in `load_or_compute()` so
repeated executions skip the expensive permutation step:

```r
gsea_result <- load_or_compute(
  checkpoint_file = "1.9_gsea_transportdb.rds",
  description     = "TransportDB GSEA",
  compute_fn      = function() { ... }
)
```

---

## 6. Result Normalization to Master Table Schema

All GSEA results -- whether from MSigDB, MitoPathways, TransportDB, or GATOM
-- are normalized to a single schema before entering the master table.

### The master_gsea_table schema

| Column | Type | Description |
|--------|------|-------------|
| `pathway_id` | character | Original gene set ID (with database prefix) |
| `pathway_name` | character | Cleaned, human-readable name |
| `nes` | numeric | Normalized Enrichment Score |
| `pvalue` | numeric | Raw p-value |
| `padj` | numeric | FDR-adjusted p-value |
| `set_size` | integer | Number of genes in the gene set |
| `leading_edge_size` | integer | Number of genes in the leading edge |
| `gene_ratio` | numeric | `leading_edge_size / set_size` |
| `core_enrichment` | character | Leading edge genes, slash-separated |
| `database` | character | Source database name (e.g. "Hallmark", "TransportDB") |
| `contrast` | character | Contrast name (e.g. "IL2RAKO_vs_NTC") |
| `neg_log_padj` | numeric | `-log10(padj)`, capped at 16 |
| `direction` | character | "Up" or "Down" based on NES sign |

### Normalization for clusterProfiler results

The `normalize_gsea_results()` function in `02_analysis/helpers/normalize_gsea.R`
converts a `gseaResult` object to this schema:

```r
normalized <- normalize_gsea_results(
  gsea_obj     = gsea_result,
  database     = "TransportDB",
  contrast     = "IL2RAKO_vs_NTC",
  padj_cutoff  = 1,           # Keep all results
  strip_prefix = TRUE
)
```

This function:
1. Extracts the `@result` slot from the gseaResult object.
2. Computes `leading_edge_size` by counting slash-separated genes in `core_enrichment`.
3. Computes `gene_ratio = leading_edge_size / setSize`.
4. Cleans pathway names via `clean_pathway_name()` (strips prefixes, title-cases, truncates to 60 chars).
5. Adds `neg_log_padj` and `direction` columns.

### Normalization for non-GSEA results (GATOM)

GATOM modules do not produce gseaResult objects. The `1.10.gatom_to_gsea.R`
script manually constructs rows matching the schema:

```r
entry <- data.frame(
  pathway_id        = naming_result$pathway_id,
  pathway_name      = naming_result$pathway_name,
  nes               = ifelse(is.na(avg_logfc), 0, avg_logfc),  # pseudo-NES
  pvalue            = geom_mean_pval,
  padj              = geom_mean_pval,
  set_size          = length(module_genes),
  leading_edge_size = length(module_genes),
  gene_ratio        = NA_real_,
  core_enrichment   = paste(module_genes, collapse = "/"),
  database          = "GATOM",
  contrast          = contrast_name,
  neg_log_padj      = -log10(max(geom_mean_pval, 1e-50)),
  direction         = ifelse(avg_logfc > 0, "Up", "Down")
)
```

---

## 7. Master Table Integration (Idempotent Append)

### Initial creation (1.5.create_master_tables.R)

The master GSEA table is first created by `1.5.create_master_tables.R`, which
normalizes all MSigDB and mitochondrial GSEA results and writes
`master_gsea_table.csv`:

```r
all_normalized <- list()
for (db_name in names(all_gsea)) {
  all_normalized[[db_name]] <- normalize_gsea_results(gsea_obj, database, contrast)
}
master_gsea <- bind_rows(all_normalized)
write_csv(master_gsea, "03_results/tables/master_gsea_table.csv")
```

### Subsequent appends (idempotent pattern)

Scripts that run after master table creation (e.g., `1.9.transportdb_gsea.R`,
`1.10.gatom_to_gsea.R`) use an **idempotent append** pattern: they delete
any existing rows from their database before inserting new ones. This makes
the script safe to re-run without creating duplicates:

```r
master_file <- file.path(DIR_TABLES, "master_gsea_table.csv")

if (file.exists(master_file)) {
  master_df <- read_csv(master_file, show_col_types = FALSE)

  # Remove any existing entries from this database
  master_df <- master_df %>%
    dplyr::filter(database != "TransportDB")

  # Append new results
  master_df <- bind_rows(master_df, export_df)
  write_csv(master_df, master_file)
} else {
  write_csv(export_df, master_file)
}
```

The key operation is `filter(database != "TransportDB")` -- this removes all
prior TransportDB rows before appending the fresh results. The `database`
column acts as the partitioning key.

### Execution order

The intended execution order is:

```
1.1 core_pipeline.R        -> DE results + MSigDB GSEA checkpoints
1.3 mito_db_prepare.R      -> Parsed mito databases (RDS)
1.4 mito_gsea.R            -> Mito GSEA checkpoints
1.5 create_master_tables.R  -> master_gsea_table.csv (MSigDB + Mito)
1.9 transportdb_gsea.R     -> Appends TransportDB rows
1.10 gatom_to_gsea.R        -> Appends GATOM rows
```

Scripts 1.9 and 1.10 can run in any order after 1.5.

---

## 8. Recipe: Adding a New Custom Database

Follow this template to integrate any new gene set database. The example
assumes a hypothetical "MetaboAtlas" database.

### Step 1: Place raw data

```
00_data/references/metaboatlas/MetaboAtlas_v2.tsv
```

### Step 2: Write a parser (or use generic parsers)

If your file is a simple two-column (gene_set, gene) TSV/CSV, use the
toolkit's generic parser:

```r
source("01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/parse_external_genesets.R")

db <- parse_geneset_file(
  file     = "00_data/references/metaboatlas/MetaboAtlas_v2.tsv",
  gs_col   = "pathway",      # Column with gene set names
  gene_col = "gene_symbol",  # Column with gene symbols
  desc_col = "description",  # Optional: column with descriptions
  prefix   = "METABOATLAS"
)
```

For GMT files:
```r
db <- parse_gmt("path/to/file.gmt", prefix = "METABOATLAS")
```

For GMX files:
```r
db <- parse_gmx("path/to/file.gmx", prefix = "METABOATLAS")
```

If you need a custom parser, write a function that returns:
```r
list(
  T2G = data.frame(gs_name = character(), gene_symbol = character()),
  T2N = data.frame(gs_name = character(), description = character()),
  source = "MetaboAtlas v2",
  created = Sys.time()
)
```

### Step 3: Handle gene ID conversion (if needed)

If the database provides a non-symbol ID type:

```r
db$T2G <- convert_geneset_ids(
  db$T2G,
  org_db    = org.Mm.eg.db,
  from_type = "ENSEMBL",     # or "REFSEQ", "ENTREZID", etc.
  to_type   = "SYMBOL",
  drop_unmapped = TRUE
)
```

If the database is human and the analysis is mouse:

```r
db$T2G <- convert_human_to_mouse(db$T2G, drop_unmapped = TRUE)
```

### Step 4: Filter by gene set size

```r
source("02_analysis/helpers/pathway_utils.R")

filtered <- filter_pathways_by_size(db$T2G, db$T2N, min_size = 5, max_size = 500)
db$T2G <- filtered$T2G
db$T2N <- filtered$T2N
```

### Step 5: Save the processed database

```r
saveRDS(db, "00_data/references/metaboatlas/metaboatlas_genesets.rds")
```

### Step 6: Create the GSEA script

Create `02_analysis/1.XX.metaboatlas_gsea.R` following this skeleton:

```r
#!/usr/bin/env Rscript
source("02_analysis/config/config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(fgsea)
})
source_toolkit(verbose = FALSE)
source("02_analysis/helpers/pathway_utils.R")

# --- Load database ---
db <- readRDS("00_data/references/metaboatlas/metaboatlas_genesets.rds")

# --- Load DE results and create ranked list ---
de_results <- readRDS(file.path(DIR_CHECKPOINTS, CHECKPOINT_DE_RESULTS))
de_table <- de_results[[1]]
contrast_name <- names(CONTRASTS)[1]

ranked_genes <- de_table[[RANK_METRIC]]
names(ranked_genes) <- rownames(de_table)
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# --- Prepare TERM2GENE ---
T2G <- db$T2G %>%
  dplyr::rename(term = gs_name, gene = gene_symbol) %>%
  dplyr::filter(gene %in% names(ranked_genes))

T2N <- db$T2N %>%
  dplyr::rename(term = gs_name, name = description)

# --- Run GSEA ---
gsea_result <- load_or_compute(
  checkpoint_file = "1.XX_gsea_metaboatlas.rds",
  description     = "MetaboAtlas GSEA",
  compute_fn      = function() {
    set.seed(GSEA_SEED)
    clusterProfiler::GSEA(
      geneList      = ranked_genes,
      TERM2GENE     = T2G,
      TERM2NAME     = T2N,
      pvalueCutoff  = GSEA_PVALUE_CUTOFF,
      pAdjustMethod = GSEA_PADJ_METHOD,
      eps           = 0,
      by            = "fgsea",
      nPermSimple   = GSEA_NPERM,
      verbose       = FALSE
    )
  }
)

# Fix geneSets slot
gsea_result@geneSets <- split(T2G$gene, T2G$term)

# --- Normalize to master table schema ---
result_df <- gsea_result@result
if (nrow(result_df) > 0) {
  export_df <- result_df %>%
    dplyr::mutate(
      pathway_id         = ID,
      pathway_name       = Description,
      database           = "MetaboAtlas",
      contrast           = contrast_name,
      nes                = NES,
      pvalue             = pvalue,
      padj               = p.adjust,
      set_size           = setSize,
      leading_edge_size  = as.integer(gsub(".*size=(\\d+).*", "\\1", leading_edge)),
      gene_ratio         = NA_real_,
      core_enrichment    = core_enrichment,
      direction          = ifelse(NES > 0, "Up", "Down"),
      neg_log_padj       = -log10(p.adjust)
    ) %>%
    dplyr::select(
      pathway_id, pathway_name, nes, pvalue, padj, set_size,
      leading_edge_size, gene_ratio, core_enrichment, database,
      contrast, neg_log_padj, direction
    )

  # --- Idempotent append to master table ---
  master_file <- file.path(DIR_TABLES, "master_gsea_table.csv")
  if (file.exists(master_file)) {
    master_df <- read_csv(master_file, show_col_types = FALSE)
    master_df <- master_df %>% dplyr::filter(database != "MetaboAtlas")
    master_df <- bind_rows(master_df, export_df)
    write_csv(master_df, master_file)
  } else {
    write_csv(export_df, master_file)
  }
}
```

### Step 7: Update documentation

Add the new database to `CLAUDE.md`, `plan.md`, and the pipeline status table.
Register the checkpoint file name in `config.R` if it will be referenced by
other scripts.

### Checklist

- [ ] Raw data placed under `00_data/references/`
- [ ] Parser produces T2G/T2N with correct column names (`gs_name`, `gene_symbol`, `description`)
- [ ] Gene symbols match species used in DE analysis
- [ ] Gene set names carry a unique database prefix
- [ ] Size filtering applied (5-500 genes)
- [ ] Processed database saved as RDS
- [ ] GSEA script uses `load_or_compute()` for checkpoint caching
- [ ] `@geneSets` slot populated after GSEA
- [ ] Results normalized to master table schema (13 columns)
- [ ] Master table append uses idempotent filter-then-bind pattern
- [ ] Script is safe to re-run without creating duplicate rows

---

## Appendix: File Inventory

### Toolkit (reusable across projects)

| File | Functions |
|------|-----------|
| `01_scripts/RNAseq-toolkit/scripts/GSEA/GSEA_processing/parse_external_genesets.R` | `parse_transportdb()`, `parse_geneset_file()`, `parse_gmt()`, `parse_gmx()`, `parse_mitoxplorer()`, `convert_human_to_mouse()`, `convert_geneset_ids()` |

### Project-specific helpers

| File | Functions |
|------|-----------|
| `02_analysis/helpers/parse_mito_databases.R` | `parse_mitopathways_gmx()`, `parse_mitoxplorer()`, `convert_human_to_mouse_symbols()` |
| `02_analysis/helpers/parse_transportdb.R` | `parse_transportdb()`, `convert_refseq_to_symbols()`, `parse_transportdb_all()` |
| `02_analysis/helpers/pathway_utils.R` | `filter_pathways_by_size()`, `get_pathway_sizes()`, `calculate_pathway_similarity()`, `merge_redundant_pathways()` |
| `02_analysis/helpers/normalize_gsea.R` | `normalize_gsea_results()`, `empty_gsea_tibble()`, `clean_pathway_name()`, `smart_capitalize()` |

### Pipeline scripts

| File | Purpose |
|------|---------|
| `02_analysis/1.3.mito_db_prepare.R` | Parse and merge MitoPathways + mitoXplorer |
| `02_analysis/1.4.mito_gsea.R` | Run GSEA on 3 mito database variants |
| `02_analysis/1.5.create_master_tables.R` | Create master_gsea_table.csv (MSigDB + Mito) |
| `02_analysis/1.9.transportdb_gsea.R` | TransportDB GSEA + master table append |
| `02_analysis/1.10.gatom_to_gsea.R` | GATOM module extraction + master table append |

### Reference data

```
00_data/references/
  mitochonria/
    MitoPathways3.0.gmx              # Human GMX (raw input)
    mitoXplorer3.0/
      mouse_gene_function.txt         # Mouse TSV (raw input)
    mito_mitopathways.rds             # Parsed MitoPathways (T2G/T2N)
    mito_mitoxplorer.rds              # Parsed mitoXplorer (T2G/T2N)
    unified_mito_pathways.rds         # Merged, deduplicated (T2G/T2N)
    unified_mito_pathways.gmx         # Merged, GMX export
  TransportDB2.0.csv                  # Raw input (headerless CSV)
  transportdb_genesets.rds            # Parsed TransportDB (T2G/T2N)
  transportdb_genesets.gmx            # Parsed TransportDB, GMX export
  gatom/
    network.kegg.rds                  # GATOM network definitions
    network.combined.rds
    org.Mm.eg.gatom.anno.rds          # GATOM annotation
```
