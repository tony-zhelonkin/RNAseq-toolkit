# 04 - Output Artifacts and Visualization Pipeline

This document describes how GSEA output artifacts are produced by the R visualization
scripts (`2.1.visualizations.R`, `2.2.gsea_viz.R`) and consumed by the interactive
pathway explorer (`3.1.pathway_explorer.py`). It covers file formats, directory layout,
CSV schemas, color systems, and extension patterns.

---

## Table of Contents

1. [Pipeline Flow](#pipeline-flow)
2. [Output Artifact Inventory](#output-artifact-inventory)
3. [CSV Schema Documentation](#csv-schema-documentation)
4. [Visualization Script Details](#visualization-script-details)
5. [Pathway Explorer Architecture](#pathway-explorer-architecture)
6. [Theming and Color Palettes](#theming-and-color-palettes)
7. [Adding a New Database](#adding-a-new-database)

---

## Pipeline Flow

```
Phase 1: Core Compute (R)
    1.1.core_pipeline.R
        |
        v
    checkpoints/*.rds          <-- gseaResult objects per database
        |                          (clusterProfiler S4 objects)
        |
        +---> 1.4.mito_gsea.R ----------> 1.4_*.rds  (Mito variants)
        +---> 1.9.transportdb_gsea.R ----> 1.9_*.rds  (TransportDB)
        +---> 1.5/1.6 master tables -----> tables/master_*.csv

Phase 2: Master Tables (R/Python bridge)
    tables/
        master_gsea_table.csv          <-- All GSEA results (all databases)
        master_de_table.csv            <-- DE results with annotations
        master_tf_activities.csv       <-- DecoupleR TF activity scores
        master_progeny_activities.csv  <-- PROGENy signaling scores
        master_gsea_significant.csv    <-- FDR-filtered subset

Phase 3: R Visualizations
    2.1.visualizations.R   (basic volcano, dotplot, heatmap, running sum)
        |
    2.2.gsea_viz.R         (comprehensive per-database GSEA plots)
        |
        v
    plots/GSEA/{Contrast}/{Database}/*.pdf

Phase 4-5: Python Interactive
    3.1.pathway_explorer.py  -->  pathway_explorer/ package
        |
        reads: tables/master_gsea_table.csv
               tables/master_tf_activities.csv
               tables/master_progeny_activities.csv
               tables/master_de_table.csv
        |
        v
    interactive/pathway_explorer.html
```

### Data Flow Diagram (Text)

```
 [RDS checkpoints]
       |
       |  readRDS()              [master CSV tables]
       v                                |
 +------------------+                   |  pd.read_csv()
 | 2.1 / 2.2 R Viz  |                  v
 |   (per-database)  |          +-------------------+
 +--------+---------+          | 3.1 Pathway        |
          |                    | Explorer (Python)   |
          v                    +--------+------------+
   plots/GSEA/                          |
   {Contrast}/                          v
   {Database}/               interactive/
    *.pdf                    pathway_explorer.html
```

Key insight: R visualization scripts consume **RDS checkpoint objects** directly,
while the Python pathway explorer consumes **CSV master tables**. The master tables
serve as the bridge between the two languages.

---

## Output Artifact Inventory

### Directory: `03_results/plots/`

```
plots/
+-- GSEA/
|   +-- IL2RAKO_vs_NTC/                    # Per-contrast root
|       +-- gsea_summary_barplot.pdf        # Cross-database summary (from 2.1)
|       +-- H/                             # Hallmark
|       |   +-- IL2RAKO_vs_NTC_H_dotplot.pdf
|       |   +-- IL2RAKO_vs_NTC_H_up_dot.pdf
|       |   +-- IL2RAKO_vs_NTC_H_down_dot.pdf
|       |   +-- IL2RAKO_vs_NTC_H_nes_bar.pdf
|       |   +-- IL2RAKO_vs_NTC_H_running_sum.pdf    # Top-5 stacked
|       |   +-- IL2RAKO_vs_NTC_H_results.txt         # Text summary
|       |   +-- dotplot_H.pdf                         # From 2.1 (legacy)
|       |   +-- dotplot_H_up.pdf                      # From 2.1 (legacy)
|       |   +-- dotplot_H_down.pdf                    # From 2.1 (legacy)
|       |   +-- running_sum/                          # Individual plots
|       |       +-- running_HALLMARK_OXIDATIVE_PHOSPHORYLATION.pdf
|       |       +-- running_HALLMARK_GLYCOLYSIS.pdf
|       |       +-- ...
|       +-- C2_KEGG/                       # Same structure as H/
|       +-- C2_REACTOME/
|       +-- C2_WIKIPATHWAYS/
|       +-- C3_TF/
|       +-- C5_BP/
|       +-- C5_MF/
|       +-- C5_CC/
|       +-- MitoPathways/                  # Custom databases
|       +-- MitoXplorer/
|       +-- Mito_Unified/
|       +-- TransportDB/
|       +-- cross_database_pooled/
|           +-- pooled_IL2RAKO_vs_NTC.pdf         # All databases
|           +-- focused_top5_IL2RAKO_vs_NTC.pdf   # Key databases, top 5
|           +-- focused_top10_IL2RAKO_vs_NTC.pdf  # Key databases, top 10
|           +-- pooled_pathways_data.csv           # Underlying data
+-- Volcano/
|   +-- volcano_IL2RA_KO - NTC.pdf
|   +-- volcano_IL2RA_KO - NTC.png
+-- Heatmap/
|   +-- top_de_genes_heatmap.pdf
+-- GATOM/
|   +-- combined_module_enzymes.{pdf,png,svg}
|   +-- combined_module_ggraph.{png,svg}
|   +-- combined_module_robustness.{pdf,png,svg}
|   +-- kegg_module_*.{pdf,png,svg}
|   +-- sensitivity_panel.{pdf,png,svg}
|   +-- transport_network.{pdf,png,svg}
+-- PROGENy/
|   +-- progeny_activity_barplot.{pdf,png}
|   +-- progeny_activity_volcano.{pdf,png}
|   +-- progeny_targets_*.{pdf,png}
+-- TF/
|   +-- tf_activity_barplot.pdf
|   +-- tf_activity_distribution.pdf
|   +-- tf_activity_volcano.pdf
|   +-- tf_direction_summary.pdf
+-- QC/
|   +-- library_sizes.pdf
|   +-- mean_variance_plot.pdf
|   +-- PCA_plot.pdf
|   +-- pvalue_histogram.pdf
|   +-- sample_correlation_heatmap.pdf
|   +-- voom_mean_variance.pdf
+-- Il2ra/
    +-- il2ra_context_barplot.pdf
    +-- il2ra_counts_barplot.pdf
    +-- il2ra_expression_boxplot.pdf
    +-- il2ra_expression_samples.pdf
    +-- il2ra_sashimi_custom.pdf
```

### Directory: `03_results/interactive/`

```
interactive/
+-- pathway_explorer.html      # Self-contained interactive dashboard (~7.8 MB)
+-- docs/                      # Planning documents
```

### Directory: `03_results/tables/`

| File | Size | Purpose |
|------|------|---------|
| `master_gsea_table.csv` | ~2.7 MB | All GSEA results across all databases |
| `master_de_table.csv` | ~3.3 MB | Full DE results with Ensembl annotations |
| `master_gsea_significant.csv` | ~784 KB | FDR-filtered GSEA subset |
| `master_tf_activities.csv` | ~283 KB | DecoupleR/CollecTRI TF activities |
| `master_progeny_activities.csv` | ~35 KB | PROGENy pathway activities |
| `significant_de_genes.csv` | ~1.6 MB | FDR-significant DE genes only |
| `de_summary.csv` | ~83 B | One-row summary (counts of sig genes) |
| `gsea_summary_stats.csv` | ~1 KB | Per-database summary statistics |

### Directory: `03_results/tables/GSEA_results/`

Per-database GSEA exports (used for supplementary analysis, not by the main pipeline):

| File | Description |
|------|-------------|
| `mito_gsea_results.csv` | All mitochondrial GSEA results |
| `mito_gsea_significant.csv` | Significant mitochondrial results (FDR < 0.05) |
| `mito_gsea_unified.csv` | Unified Mito database results |
| `mito_gsea_unified_significant.csv` | Significant unified Mito results |
| `mito_gsea_mitopathways.csv` | MitoPathways-specific results |
| `mito_gsea_mitopathways_significant.csv` | Significant MitoPathways results |
| `mito_gsea_mitoxplorer.csv` | MitoXplorer-specific results |
| `mito_gsea_mitoxplorer_significant.csv` | Significant MitoXplorer results |
| `mito_gsea_combined_significant.csv` | All significant mito results combined |
| `transportdb_gsea.csv` | TransportDB family GSEA results |
| `gatom_module_genesets.csv` | Gene sets derived from GATOM active modules |

### Directory: `03_results/tables/DE_results/`

| File | Description |
|------|-------------|
| `IL2RAKO_vs_NTC_annotated.csv` | Annotated DE table for the primary contrast |

### Directory: `03_results/tables/GATOM/`

| File | Description |
|------|-------------|
| `active_module_genes.csv` | Genes in GATOM active modules |
| `network_comparison.csv` | Network topology comparison across k-values |
| `k_sensitivity_summary.csv` | K-parameter sensitivity summary |
| `core_genes_robust.csv` | Robust core genes (found across all k-values) |
| `gene_robustness_scorecard.csv` | Per-gene robustness scoring |
| `k_parameter_comparison.csv` | Full k-parameter comparison |
| `transport_genes.csv` | Transport-related genes from GATOM |

### Per-Database File Naming Convention

Files generated by `2.2.gsea_viz.R` follow this pattern:

```
{contrast}_{database}_{plot_type}.pdf

Examples:
  IL2RAKO_vs_NTC_H_dotplot.pdf
  IL2RAKO_vs_NTC_C2_KEGG_up_dot.pdf
  IL2RAKO_vs_NTC_MitoPathways_nes_bar.pdf
  IL2RAKO_vs_NTC_TransportDB_running_sum.pdf
```

Plot types per database:
- `_dotplot.pdf` -- Combined (up + down) dotplot
- `_up_dot.pdf` -- Upregulated-only dotplot (NES > 0)
- `_down_dot.pdf` -- Downregulated-only dotplot (NES < 0)
- `_nes_bar.pdf` -- NES barplot (horizontal bars colored by direction)
- `_running_sum.pdf` -- Stacked running sum (top 5 by |NES|)
- `_results.txt` -- Text summary with counts and top pathways
- `running_sum/running_{PATHWAY_ID}.pdf` -- Individual running sum (top 10)

---

## CSV Schema Documentation

### master_gsea_table.csv

The primary GSEA results table consumed by the pathway explorer.

| Column | Type | Description |
|--------|------|-------------|
| `pathway_id` | string | Unique pathway identifier (e.g., `HALLMARK_GLYCOLYSIS`) |
| `pathway_name` | string | Human-readable name (same as ID for MSigDB; descriptive for custom) |
| `nes` | float | Normalized Enrichment Score |
| `pvalue` | float | Raw p-value from GSEA |
| `padj` | float | FDR-adjusted p-value |
| `set_size` | int | Number of genes in the gene set |
| `leading_edge_size` | int | Number of genes in the leading edge |
| `gene_ratio` | float | leading_edge_size / set_size |
| `core_enrichment` | string | Slash-separated gene symbols (`Hk2/Ldha/Slc2a1/...`) |
| `database` | string | Source database (e.g., `Hallmark`, `KEGG`, `MitoPathways`, `GATOM`) |
| `contrast` | string | Contrast name (e.g., `IL2RAKO_vs_NTC`) |
| `neg_log_padj` | float | -log10(padj) for visualization |
| `direction` | string | `Up` or `Down` based on NES sign |
| `module_theme` | string | (GATOM only) Module theme annotation |
| `n_metabolites` | int | (GATOM only) Number of metabolites in module |
| `n_enzymes` | int | (GATOM only) Number of enzymes in module |

Schema version: 1.0.0 (defined in `config/pipeline.yaml`).

### master_de_table.csv

| Column | Type | Description |
|--------|------|-------------|
| `contrast` | string | Contrast name |
| `gene_symbol` | string | Gene symbol (e.g., `Il2ra`) |
| `ensembl_id` | string | Ensembl gene ID with version |
| `ensembl_id_base` | string | Ensembl gene ID without version |
| `logFC` | float | Log2 fold change |
| `AveExpr` | float | Average expression |
| `t` | float | Moderated t-statistic (limma) |
| `P.Value` | float | Raw p-value |
| `adj.P.Val` | float | FDR-adjusted p-value |
| `B` | float | B-statistic (log-odds of DE) |

### master_tf_activities.csv / master_progeny_activities.csv

These share an identical schema (validated against `master_tf_activities` schema):

| Column | Type | Description |
|--------|------|-------------|
| `pathway_id` | string | TF name or signaling pathway |
| `pathway_name` | string | Display name |
| `database` | string | `CollecTRI` (TFs) or `PROGENy` (signaling) |
| `contrast` | string | Contrast name |
| `nes` | float | Activity score (z-score-like) |
| `pvalue` | float | Raw p-value |
| `padj` | float | FDR-adjusted p-value |
| `set_size` | int | Number of target genes |
| `leading_edge_size` | int | Number of significant targets |
| `direction` | string | `Up` or `Down` |
| `core_enrichment` | string | Slash-separated target genes |

### gsea_summary_stats.csv

| Column | Type | Description |
|--------|------|-------------|
| `database` | string | Database name |
| `contrast` | string | Contrast name |
| `direction` | string | `Up` or `Down` |
| `n_total` | int | Total pathways tested |
| `n_significant` | int | FDR-significant pathways |
| `mean_abs_nes` | float | Mean |NES| of significant pathways |

### de_summary.csv

| Column | Type | Description |
|--------|------|-------------|
| `contrast` | string | Contrast name |
| `total_genes` | int | Total genes tested |
| `sig_genes` | int | FDR-significant genes |
| `sig_up` | int | Upregulated |
| `sig_down` | int | Downregulated |

---

## Visualization Script Details

### 2.1.visualizations.R (Basic Visualizations)

**Inputs:** RDS checkpoints loaded via `readRDS()`:
- `1.1_de_results.rds` -- DE results (named list of data frames per contrast)
- `1.1_all_gsea_results.rds` -- MSigDB GSEA results (named list of gseaResult objects)
- `1.1_dge_normalized.rds` -- Normalized DGEList for heatmap expression

**Plot types generated:**

1. **Volcano plots** (`plots/Volcano/volcano_{contrast}.{pdf,png}`)
   - Uses toolkit function `create_standard_volcano()`
   - Highlights biologically relevant genes (IL-2/STAT5, glycolysis, mitochondrial)
   - Dual format: PDF (vector) + PNG (300 DPI)

2. **GSEA dotplots** (`plots/GSEA/{contrast}/{db}/dotplot_{db}[_{direction}].pdf`)
   - Uses toolkit function `gsea_dotplot()`
   - Three variants: combined, upregulated-only (`_up`), downregulated-only (`_down`)
   - Filtered by `GSEA_FDR_CUTOFF` (0.05)
   - Shows up to `N_TOP_PATHWAYS` (20)

3. **GSEA summary barplot** (`plots/GSEA/{contrast}/gsea_summary_barplot.pdf`)
   - Bidirectional bar chart: upregulated counts right, downregulated left
   - Colors: `#fc8d59` (up/orange), `#91bfdb` (down/blue)
   - One bar per database

4. **Running sum plots** (`plots/GSEA/{contrast}/H/running_sum/running_{pathway}.pdf`)
   - Only for Hallmark (H) database in 2.1
   - Uses toolkit function `gsea_running_sum_plot()`
   - Key pathways hardcoded: OXPHOS, Glycolysis, MYC, mTORC1, IL2-STAT5, Apoptosis

5. **Top DE genes heatmap** (`plots/Heatmap/top_de_genes_heatmap.pdf`)
   - Top 25 up + top 25 down by FDR
   - Z-score normalized expression
   - Color: `colorRampPalette(c("#91bfdb", "white", "#fc8d59"))`
   - Sample annotation bar: NTC = `#1f77b4`, IL2RA_KO = `#ff7f0e`

### 2.2.gsea_viz.R (Comprehensive GSEA Visualization)

**Inputs:** RDS checkpoints loaded via `readRDS()`:
- `1.1_all_gsea_results.rds` -- MSigDB databases (H, C2_KEGG, C2_REACTOME, etc.)
- `1.4_all_mito_gsea.rds` -- All mitochondrial variants (unified, MitoPathways, MitoXplorer)
- `1.9_gsea_transportdb.rds` -- TransportDB GSEA results

**Processing pattern:** Identical 7-step procedure applied to every database:

```r
# For each database:
# 1. Combined dotplot         ->  {prefix}_dotplot.pdf
# 2. Upregulated dotplot      ->  {prefix}_up_dot.pdf
# 3. Downregulated dotplot    ->  {prefix}_down_dot.pdf
# 4. NES barplot              ->  {prefix}_nes_bar.pdf
# 5. Stacked running sum      ->  {prefix}_running_sum.pdf  (top 5 by |NES|)
# 6. Individual running sums  ->  running_sum/running_{safe_name}.pdf  (top 10)
# 7. Text summary             ->  {prefix}_results.txt
```

**Database-specific handling:**

- **MSigDB** (H, C2_*, C3_TF, C5_*): Direct processing from `all_gsea` list.
  Pathway IDs used directly for running sum plots.
- **Mitochondrial** (MitoPathways, MitoXplorer, Mito_Unified): Requires a
  workaround where `Description` is temporarily set to `ID` for running sum
  compatibility with `enrichplot::gseaplot2()`. Original descriptions preserved
  via `labels` parameter.
- **TransportDB**: Same workaround as mitochondrial databases.

**Cross-database pooled dotplot:**
Collects top 10 significant pathways from each database, formats names, and
creates a unified dotplot with:
- X-axis: Gene Ratio (leading edge / set size)
- Y-axis: Pathway name (prefixed with database)
- Fill color: NES gradient (Blue-White-Orange via `nes_scale()`)
- Size: -log10(FDR)
- Outline: Black border on pathways with FDR < 0.01

Three variants saved:
- `pooled_{contrast}.pdf` -- All databases
- `focused_top5_{contrast}.pdf` -- Key databases, top 5 per db
- `focused_top10_{contrast}.pdf` -- Key databases, top 10 per db

---

## Pathway Explorer Architecture

The interactive pathway explorer (`3.1.pathway_explorer.py`) is a thin wrapper
around the `02_analysis/pathway_explorer/` package, which consists of:

### Module Structure

```
pathway_explorer/
+-- __init__.py          # Public API exports
+-- __main__.py          # CLI entry point
+-- main.py              # Orchestration: load -> similarity -> embed -> HTML
+-- config.py            # Paths, YAML loading, schema validation, constants
+-- data_loader.py       # CSV loading, score standardization
+-- similarity.py        # Jaccard/Overlap similarity matrices
+-- embedding.py         # UMAP/t-SNE/PCA dimensionality reduction
+-- html_generator.py    # Self-contained HTML dashboard generation
```

### Data Ingestion Pipeline

1. **Load master tables** (`data_loader.py`):
   - `master_gsea_table.csv` -- all pathway GSEA results
   - `master_tf_activities.csv` -- CollecTRI TF activities (concatenated)
   - `master_progeny_activities.csv` -- PROGENy signaling (concatenated)
   - `master_de_table.csv` -- gene rankings for running sum display

2. **Schema validation** (`config.py`):
   Validates required columns from `pipeline.yaml` schema definitions:
   ```python
   validate_schema(df, 'master_gsea_table', SCHEMAS)
   # Checks: pathway_id, pathway_name, database, nes, pvalue, padj, core_enrichment
   ```

3. **Database reclassification** (`data_loader.py`):
   Generic `Mitochondria` database entries are reclassified based on ID prefix:
   - `MITOPATHWAYS_*` --> `MitoPathways`
   - `MITOXPLORER_*` --> `MitoXplorer`

4. **Entity type assignment** (`data_loader.py`):
   - `CollecTRI` --> `TF`
   - `PROGENy` --> `PROGENy`
   - Everything else --> `Pathway`

5. **Score standardization** (`data_loader.py`):
   Creates `signed_sig = -log10(padj) * sign(nes)`, capped at +/-50.
   This unifies TF activity scores and GSEA NES onto a comparable scale.

6. **Gene set parsing**:
   `core_enrichment` column split on `/` to create Python sets for similarity.

### Similarity Computation (`similarity.py`)

Uses a **hybrid similarity strategy**:
- **Same-type** (Pathway-Pathway, TF-TF, PROGENy-PROGENy): Jaccard index
  `|A intersection B| / |A union B|`
- **TF cross-type** (TF-Pathway, TF-PROGENy): Overlap coefficient
  `|A intersection B| / min(|A|, |B|)`
- **PROGENy-Pathway**: Jaccard (both are pathway-like)

Implementation uses sparse binary gene-entity matrix for efficient computation.

Neighbors extracted: top-5 per pathway with minimum similarity >= 0.15.

### Embedding (`embedding.py`)

Converts similarity matrix to 2D coordinates:
- Preferred: **UMAP** (metric=precomputed, n_neighbors=15, min_dist=0.1)
- Fallback chain: UMAP > PCA > random projection
- Output normalized to [0, 1] range

### HTML Dashboard (`html_generator.py`)

Generates a self-contained HTML file (~7.8 MB) with embedded:
- Plotly.js (loaded from CDN)
- Pathway data as JSON
- Gene rankings for running sum computation

Interactive features:
- 2D scatter plot with UMAP-embedded pathways
- Color encoding: NES (diverging blue-white-orange)
- Size encoding: -log10(FDR)
- Shape encoding by entity type (circle=Pathway, diamond=TF, square=PROGENy)
- Edge overlay: on-click shows connected pathways (high gene overlap)
- Sidebar filters: database checkboxes, FDR slider, NES threshold, search
- Rich tooltips: pathway name, database, NES, FDR, set size, leading edge genes
- Running sum enrichment plot (computed client-side from gene rankings)
- Gene table for selected pathway

---

## Theming and Color Palettes

### Single Source of Truth

Colors are defined in **two locations** that stay synchronized:

1. **`config/pipeline.yaml`** -- shared YAML read by both R and Python
2. **`config/color_config.R`** -- R-specific helper functions that load from YAML

The Python config (`pathway_explorer/config.py`) also reads `pipeline.yaml`.

### Diverging Color Scale (NES / logFC)

Used across all NES-colored plots, heatmaps, and the pathway explorer.

| Role | Hex | Visual |
|------|-----|--------|
| Negative (downregulated) | `#2166AC` | Blue |
| Neutral | `#F7F7F7` | White |
| Positive (upregulated) | `#B35806` | Orange/Brown |

5-point gradient for heatmaps:
- `#2166AC` -- `#92c5de` -- `#F7F7F7` -- `#f4a582` -- `#B35806`

R helper functions:
```r
nes_color_scale(limits = c(-3.5, 3.5))       # ComplexHeatmap colorRamp2
nes_ggplot_scale(limits, name = "NES")        # ggplot2 scale_color_gradient2
nes_ggplot_fill_scale(limits, name = "NES")   # ggplot2 scale_fill_gradient2
get_diverging_palette(n = 100)                # pheatmap-compatible vector
```

### Experimental Group Colors

| Group | Hex | Visual |
|-------|-----|--------|
| NTC | `#1f77b4` | Blue |
| IL2RA_KO | `#ff7f0e` | Orange |

### Database Colors (Okabe-Ito Colorblind-Safe)

Used for categorical database coloring in cross-database plots and the pathway explorer.

| Database | Hex | Visual |
|----------|-----|--------|
| Hallmark | `#E69F00` | Orange |
| KEGG | `#56B4E9` | Sky Blue |
| Reactome | `#009E73` | Bluish Green |
| WikiPathways | `#F0E442` | Yellow |
| GO_BP | `#0072B2` | Blue |
| GO_MF | `#D55E00` | Vermillion |
| GO_CC | `#CC79A7` | Reddish Purple |
| MitoPathways | `#332288` (R) / `#5E4FA2` (Py) | Indigo / Deep Violet |
| MitoXplorer | `#999999` (R) / `#9E9AC8` (Py) | Gray / Periwinkle |
| CollecTRI | `#882255` | Wine |
| PROGENy | `#44AA99` | Teal |
| TransportDB | `#DDCC77` | Sand |
| GATOM | `#117733` | Forest |

Note: MitoPathways and MitoXplorer have slightly different hex values between
R (`color_config.R`) and Python (`pathway_explorer/config.py`). The Python
explorer uses violet-family hues to signal mito-database relatedness.

### Running Sum Plot Palette

9-color palette for multi-pathway running sum overlay plots:

```r
RUNNING_SUM_PALETTE <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#FFFF33",  # Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#999999"   # Gray
)
```

### Robustness Tier Colors (GATOM)

| Tier | Hex | Width | Alpha |
|------|-----|-------|-------|
| core | `#2166AC` | 3.0 | 1.0 |
| moderate | `#92C5DE` | 2.0 | 0.8 |
| exploratory | `#D1E5F0` | 1.0 | 0.5 |

### Entity Type Shapes (Pathway Explorer)

Defined in `pipeline.yaml`:
- Pathway: circle
- TF: diamond
- PROGENy: square
- GATOM: hexagon (reserved for future use)

### Custom Minimal Theme (`custom_minimal_theme_with_grid`)

The RNAseq-toolkit provides a shared ggplot2 theme sourced by all visualization scripts:

```r
custom_minimal_theme_with_grid(base_size = 12, base_family = "")
```

Properties:
- Based on `theme_classic()`
- White background with no grid lines (major or minor)
- Black axis lines and ticks at 0.5pt width
- No panel border
- Centered plot title with bottom margin
- 10px margins on all sides
- Axis title margins: x-title 10pt above, y-title 10pt right

### Output Format Conventions

| Property | Standard | Set in |
|----------|----------|--------|
| Primary static format | PDF (vector) | All viz scripts |
| Secondary raster format | PNG at 300 DPI | `PLOT_DPI` in config |
| Tertiary format | SVG (GATOM network plots only) | GATOM scripts |
| Default plot size | 10 x 8 inches | `PLOT_WIDTH`, `PLOT_HEIGHT` |
| Dotplot size | 12 x 10 inches | Hardcoded in 2.2 |
| Heatmap size | 10 x 12 inches | Hardcoded in 2.1 |
| Volcano size | 10 x 8 inches | Hardcoded in 2.1 |
| Running sum size | 12 x 10 (stacked), 10 x 8 (individual) | Hardcoded in 2.2 |

Plots are saved with `pdf()` / `png()` (in 2.1) or `ggsave()` (in 2.2). The `ggsave()`
approach is preferred as it handles device management automatically.

---

## Adding a New Database

To add a new database (e.g., `MyCustomDB`) to the full visualization pipeline:

### Step 1: Create GSEA checkpoint

Run GSEA and save as RDS checkpoint:

```r
mycustomdb_gsea <- run_gsea(
  DE_results = de_table,
  rank_metric = "t",
  species = "Mus musculus",
  custom_gene_sets = my_gene_sets  # Named list of character vectors
)
saveRDS(mycustomdb_gsea, file.path(DIR_CHECKPOINTS, "1.X_gsea_mycustomdb.rds"))
```

### Step 2: Add to master_gsea_table.csv

Ensure results are appended to the master GSEA table with the correct schema:

```r
# Use normalize_gsea_results() from toolkit
new_rows <- normalize_gsea_results(mycustomdb_gsea,
                                    database = "MyCustomDB",
                                    contrast = contrast_name)
# Append to existing master table
master <- read_csv("03_results/tables/master_gsea_table.csv")
master <- bind_rows(master, new_rows)
write_csv(master, "03_results/tables/master_gsea_table.csv")
```

### Step 3: Add to 2.2.gsea_viz.R

Add a new processing block (following the TransportDB pattern):

```r
# Load checkpoint
mycustomdb_gsea <- tryCatch({
  readRDS(file.path(DIR_CHECKPOINTS, "1.X_gsea_mycustomdb.rds"))
}, error = function(e) NULL)

# Process (same 7-step pattern as other custom databases)
if (!is.null(mycustomdb_gsea) && nrow(mycustomdb_gsea@result) > 0) {
  db_display <- "MyCustomDB"
  db_dir <- file.path(dir_contrast, db_display)
  # ... (follow TransportDB block pattern)

  # IMPORTANT: For custom databases, use the Description<->ID workaround
  # for running sum plots:
  mycustomdb_plot <- mycustomdb_gsea
  original_desc <- mycustomdb_gsea@result$Description
  names(original_desc) <- mycustomdb_gsea@result$ID
  mycustomdb_plot@result$Description <- mycustomdb_plot@result$ID

  p_running <- gsea_running_sum_plot(
    mycustomdb_plot,
    gene_set_ids = top_ids,
    palette = RUNNING_SUM_PALETTE[seq_along(top_ids)],
    labels = original_desc[top_ids]
  )
}
```

### Step 4: Register display name and color

In `2.2.gsea_viz.R`:
```r
DB_DISPLAY_NAMES <- c(
  ...,
  MyCustomDB = "MyCustomDB"
)
```

In `config/color_config.R` and `config/pipeline.yaml`:
```r
DATABASE_COLORS <- c(
  ...,
  MyCustomDB = "#AABBCC"
)
```

In `pathway_explorer/config.py`:
```python
DB_COLORS = {
    ...,
    'MyCustomDB': '#AABBCC',
}
```

### Step 5: Add to cross-database pooled plot

In the `[V7]` section of `2.2.gsea_viz.R`, add the new database to the
pooled pathway collection (following the TransportDB pattern).

Optionally add to `key_dbs` for the focused variants:
```r
key_dbs <- c("Hallmark", "KEGG", ..., "MyCustomDB")
```

### Step 6: Pathway explorer picks it up automatically

The pathway explorer reads all rows from `master_gsea_table.csv` and groups
by the `database` column. As long as Step 2 is complete, the new database
will appear in the explorer with its registered color. If the database name
does not match any key in `DB_COLORS`, it will receive a default gray color.

### Step 7: Update pipeline.yaml schema (optional)

If the new database introduces additional columns, update the schema in
`config/pipeline.yaml` under `schemas.master_gsea_table.optional_columns`.

---

## Toolkit Functions Used

The visualization scripts depend on these RNAseq-toolkit functions:

| Function | File | Purpose |
|----------|------|---------|
| `create_standard_volcano()` | `DE/plot_standard_volcano.R` | Volcano plot with FDR decision |
| `gsea_dotplot()` | `GSEA/GSEA_plotting/gsea_dotplot.R` | NES-gradient dotplot |
| `gsea_barplot()` | `GSEA/GSEA_plotting/gsea_barplot.R` | Horizontal NES barplot |
| `gsea_running_sum_plot()` | `GSEA/GSEA_plotting/gsea_running_sum_plot.R` | Three-panel running sum |
| `format_pathway_name()` | `GSEA/GSEA_plotting/format_pathway_names.R` | Smart name formatting |

### Key Pattern: gsea_dotplot() Selection vs Highlighting

```r
# filterBy + showCategory = what to SHOW
# padj_cutoff + highlight_sig = what gets BLACK OUTLINE
gsea_dotplot(
  gsea_obj,
  filterBy = "p.adjust",      # Sort criterion
  showCategory = 20,           # How many to display
  padj_cutoff = 0.05,          # FDR threshold for outline
  highlight_sig = TRUE,        # Enable outline
  use_gradient = TRUE          # Continuous NES color gradient
)
```

### Key Pattern: Running Sum for Custom Databases

Custom databases (MitoPathways, MitoXplorer, TransportDB) require a
workaround because `enrichplot::gseaplot2()` uses `Description` for color
mapping, but `gsea_running_sum_plot()` names palette entries by `ID`:

```r
# 1. Save original descriptions
original_desc <- gsea_obj@result$Description
names(original_desc) <- gsea_obj@result$ID

# 2. Set Description = ID for internal color mapping
plot_obj <- gsea_obj
plot_obj@result$Description <- plot_obj@result$ID

# 3. Pass original descriptions via labels parameter
p <- gsea_running_sum_plot(
  plot_obj,
  gene_set_ids = top_ids,
  labels = original_desc[top_ids]
)
```

---

## Checkpoint File Reference

| Checkpoint | Script | Contents |
|------------|--------|----------|
| `1.1_dge_normalized.rds` | 1.1 | Normalized DGEList |
| `1.1_fit_object.rds` | 1.1 | limma model fit |
| `1.1_de_results.rds` | 1.1 | Named list of DE data frames |
| `1.1_all_gsea_results.rds` | 1.1 | Named list of gseaResult (MSigDB) |
| `1.4_gsea_mito_unified.rds` | 1.4 | Unified mito GSEA |
| `1.4_gsea_mito_mitopathways.rds` | 1.4 | MitoPathways GSEA |
| `1.4_gsea_mito_mitoxplorer.rds` | 1.4 | MitoXplorer GSEA |
| `1.4_all_mito_gsea.rds` | 1.4 | Combined mito (all variants) |
| `1.6_gatom_combined.rds` | 1.6 | GATOM metabolic modules |
| `1.9_gsea_transportdb.rds` | 1.9 | TransportDB family GSEA |

---

## Toolkit Source Loading Order

The `source_toolkit()` function (defined in `config/config.R`) sources scripts in a
specific order to satisfy dependencies:

```r
# 1. GSEA processing functions (no plotting dependencies)
"GSEA/GSEA_processing/run_gsea.R"
"GSEA/GSEA_processing/get_significant_pathways.R"
"GSEA/GSEA_processing/get_pathway_genes.R"

# 2. GSEA plotting functions (format_pathway_names MUST come first)
"GSEA/GSEA_plotting/format_pathway_names.R"     # <-- Required by gsea_dotplot
"GSEA/GSEA_plotting/gsea_dotplot.R"
"GSEA/GSEA_plotting/gsea_barplot.R"
"GSEA/GSEA_plotting/gsea_running_sum_plot.R"

# 3. DE visualization
"DE/plot_standard_volcano.R"
"DE/plotPCA.R"

# 4. Shared theme
"custom_minimal_theme.R"
```

All paths are relative to `01_scripts/RNAseq-toolkit/scripts/`. The toolkit root
is set via `DIR_TOOLKIT` in `config.R`.

Additional toolkit scripts available but not auto-sourced (use manually as needed):

| Script | Purpose |
|--------|---------|
| `GSEA/GSEA_processing/run_gsea_analysis.R` | Multi-database pipeline with auto-plotting |
| `GSEA/GSEA_processing/run_pooled_gsea.R` | Cross-contrast GSEA aggregation |
| `GSEA/GSEA_processing/normalize_gsea.R` | Convert gseaResult to standard tibble |
| `GSEA/GSEA_processing/calculate_pathway_scores.R` | Pathway activity scores |
| `GSEA/GSEA_processing/parse_external_genesets.R` | Load custom gene set formats |
| `GSEA/GSEA_processing/pathway_utils.R` | Pathway utility functions |
| `GSEA/GSEA_plotting/gsea_dotplot_facet.R` | Up/Down faceted dotplot |
| `GSEA/GSEA_plotting/gsea_dotplot_compare.R` | Multi-contrast comparison dotplot |
| `GSEA/GSEA_plotting/gsea_heatmap.R` | Pathway x sample heatmaps |
| `GSEA/GSEA_plotting/gsea_meta_heatmap.R` | Cross-database meta-heatmap |
| `GSEA/GSEA_plotting/gsea_nes_comparison.R` | NES comparison across contrasts |
| `GSEA/GSEA_plotting/gsea_path_scores_heatmap.R` | Pathway score heatmap |
| `GSEA/GSEA_plotting/gsea_path_genes_heatmap.R` | Leading edge gene heatmap |
| `GSEA/GSEA_plotting/plot_pooled_contrast_dotplot.R` | Pooled contrast visualization |
| `GSEA/GSEA_plotting/plot_all_gsea_results.R` | Batch GSEA plotting |
| `General/annotate_genes.R` | Ensembl to Symbol/ENTREZID annotation |
| `General/dge_helpers.R` | DGEList construction helpers |
| `General/io_helpers.R` | File I/O utilities |
| `ORA/run_ora.R` | Over-representation analysis |
| `ORA/ora_dotplot.R` | ORA dotplot visualization |
