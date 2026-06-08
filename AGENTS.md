# Repository Guidelines

> Single source of truth for agent/contributor guidance. `CLAUDE.md` and `GEMINI.md`
> intentionally contain only `@AGENTS.md` so every assistant reads the same instructions.

## Overview

RNAseq-toolkit is a modular R **script library** (not a formal R package) for bulk
RNA-seq differential expression (DE) and Gene Set Enrichment Analysis (GSEA). It wraps
`limma`/`edgeR` for DE and `clusterProfiler`/`msigdbr`/`fgsea` for GSEA, providing
consistent, publication-ready visualization functions. It is designed to be consumed as
a git **submodule** inside larger analysis projects.

**Key philosophy:**
- **Modular** — `source()` only what you need; there is no `library(RNAseqToolkit)`.
- **Standardized** — consistent themes/colors across the parent project.
- **Wrapper-based** — simplifies complex calls to `clusterProfiler`, `limma`, `ggplot2`.

## Project Structure & Module Organization

```
scripts/
├── General/                    # Core utilities
│   ├── annotate_genes.R        # Ensembl → Symbol/ENTREZID annotation
│   ├── dge_helpers.R           # DGEList construction (build_dge)
│   └── io_helpers.R            # File I/O utilities
├── DE/                         # Differential expression visuals
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
│   │   ├── calculate_pathway_scores.R
│   │   ├── load_reference_db.R # Load bundled reference databases
│   │   └── build_reference_databases.R # Rebuild processed RDS from raw
│   ├── GSEA_plotting/          # R visualization functions
│   │   ├── gsea_dotplot.R      # Standard dotplot
│   │   ├── gsea_dotplot_facet.R # Up/Down faceted dotplot
│   │   ├── gsea_barplot.R      # NES barplot
│   │   ├── gsea_running_sum_plot.R # Running sum enrichment plot
│   │   ├── gsea_heatmap.R      # Pathway × sample heatmaps
│   │   └── format_pathway_names.R # Clean MSigDB pathway names
│   └── GSEA_plotting_python/   # Python dotplot renderer
├── custom_minimal_theme.R      # Shared ggplot2 theme
└── utils_plotting.R / famd_plotting_rnaseq.R # Shared styling helpers

data/
└── references/                 # Bundled reference databases for GSEA
    ├── METADATA.yaml           # Database registry
    ├── mitocarta3.0/           # MitoCarta + MitoPathways (raw + processed)
    ├── mitoxplorer3.0/         # mitoXplorer (raw + processed)
    ├── mitochondria_unified/   # Merged mito databases (processed only)
    └── transportdb/            # TransportDB 2.0 (raw + processed)

examples/example_analysis.R     # Runnable reference for sourcing + running analyses
tests/                          # testthat volcano tests + visual-inspection PDFs
README.md / MIGRATION.md        # Overview/usage and upgrade notes
```

## Architecture & Key Patterns

**1. GSEA Normalize-Then-Visualize**
```r
source("scripts/GSEA/GSEA_processing/normalize_gsea.R")
df <- normalize_gsea_results(gsea_obj, database = "Hallmark", contrast = "A_vs_B")
# Returns tibble with: pathway_id, pathway_name, NES, padj, direction, etc.
```

**2. Pipeline functions auto-source dependencies**
`run_gsea_analysis()` and `run_pooled_gsea()` source their helpers automatically. Pass
`helper_root` if sourcing from a non-standard location.

**3. Decision-by-FDR volcano plots**
The volcano uses `-log10(P.Value)` on the y-axis but decides significance by `adj.P.Val`.
The dashed line is placed at the raw p-value corresponding to the FDR boundary. Supply
`fixed_p_boundary` to pin the line to a known p-value (e.g. a pre-filtered dataset)
instead of deriving it from the significant gene set.

**4. GSEA dotplot "show all, highlight significant"**
`gsea_dotplot()` separates **selection** from **highlighting**:
- `filterBy` + `showCategory` → which pathways are displayed (top N by NES, p.adjust, …)
- `padj_cutoff` + `highlight_threshold` → which pathways get a black outline (FDR < threshold)

```r
gsea_dotplot(
  gsea_obj,
  filterBy = "NES",       # sort by |NES| magnitude
  showCategory = 20,      # show top 20 pathways
  padj_cutoff = 0.10,     # black outline for FDR < 0.10
  highlight_sig = TRUE,
  use_gradient = TRUE
)
```

## Golden Path: Usage in a Parent Project

Functions are `source()`'d relative to the project root (adjust the submodule path).

```r
# Differential expression
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")
source("01_modules/RNAseq-toolkit/scripts/custom_minimal_theme.R")
volcano <- create_standard_volcano(
  de_results = topTable_results,  # needs logFC, P.Value, adj.P.Val, gene symbols
  decision_by = "fdr",
  p_cutoff = 0.05,
  fc_cutoff = 1                   # log2 scale
)
ggsave("03_results/plots/volcano_AvsB.pdf", volcano, width = 8, height = 6)

# GSEA
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_processing/run_gsea.R")
source("01_modules/RNAseq-toolkit/scripts/GSEA/GSEA_plotting/gsea_dotplot.R")
gsea_res <- run_gsea(
  DE_results = de_table,          # rownames = gene symbols, col = rank_metric
  rank_metric = "t",
  species = "Mus musculus",
  category = "H"                  # Hallmark
)
dotplot <- gsea_dotplot(gsea_res, showCategory = 20, padj_cutoff = 0.05)
```

## Build, Test, and Development Commands

- Run volcano test suite (auto + visual): `Rscript tests/test_volcano_plots.R`
  (writes PDFs to `tests/output/` if enabled).
- Other suites: `Rscript tests/test_gsea_dotplot.R`, `Rscript tests/test_pathway_formatting.R`.
- Quick example workflow: `Rscript examples/example_analysis.R` (run from repo root;
  adjust input paths inside the script).
- Ad-hoc sourcing: `source("scripts/GSEA/GSEA_processing/run_gsea_analysis.R")`, then call
  helpers per README examples.

## Key Functions

| Function | Purpose | Input → Output |
|----------|---------|----------------|
| `run_gsea()` | Single GSEA analysis | DE table → gseaResult |
| `run_gsea_analysis()` | Multi-database GSEA + plots | DE table → list of gseaResult |
| `normalize_gsea_results()` | Standardize GSEA output | gseaResult → tibble |
| `create_standard_volcano()` | DE volcano plot | DE table → ggplot |
| `gsea_dotplot()` | GSEA dotplot | gseaResult → ggplot |
| `build_dge()` | Construct DGEList | count matrix + metadata → DGEList |
| `load_reference_db()` | Load bundled reference database | database name → T2G/T2N list |
| `list_reference_dbs()` | List available reference databases | → data frame |
| `download_gatom_references()` | Download GATOM network files | dest_dir → file paths |

## Coding Style & Naming Conventions

- Language: R. Indent 2 spaces; no tabs. Prefer tidyverse style (pipes / `|>`),
  consistent spacing around `=` in args.
- Functions/objects: `snake_case`; exported helpers start lowercase (`run_gsea`,
  `create_volcano_plot`). File names mirror the main function.
- Keep plotting side effects optional: return plot objects; gate file I/O behind explicit
  flags/paths.
- Add lightweight comments for non-obvious wrangling/plotting; keep docstrings minimal but precise.

## Testing Guidelines

- Primary coverage: volcano alignment and edge cases via `testthat` in
  `tests/test_volcano_plots.R`. Expect PASS and PDFs whose dashed thresholds align with
  color boundaries.
- Treat plotting changes as **visual regression**: inspect PDFs in `tests/output/`.
- When adding DE/GSEA helpers, add a minimal reproducible fixture to `tests/` and extend
  the volcano test script or create a sibling `test_*.R`.
- Clean up large artifacts; keep generated outputs in `tests/output/` (or a user-specified
  directory), not versioned.

## Required R Packages

- **Core:** `limma`, `edgeR`, `dplyr`, `tibble`, `ggplot2`
- **GSEA:** `clusterProfiler`, `msigdbr`, `enrichplot`, `fgsea`
- **Annotation:** `org.Mm.eg.db` (mouse), `org.Hs.eg.db` (human), `biomaRt`
- **Visualization:** `ggrepel`, `pheatmap`, `plotly`, `scales`

Do not assume libraries are pre-loaded — load them explicitly in scripts or ensure the
parent environment provides them. Document any new package imports.

## Git Branching, Versioning & Submodule Pinning

The toolkit uses a **two-branch model**:
- **`main`** — stable, released state.
- **`dev`** — integration branch for in-progress work.

Releases are marked with annotated **tags** (semver, e.g. `v0.2.0`). Project-specific
`dev-{project}` branches are **no longer used**; per-project work lands on `dev` and is
released via a tag.

**Pinning in a parent project:** parent repos consume the toolkit as a submodule pinned to
a **tag** (recorded as the submodule's gitlink commit), independent of which branch the tag
sits on. To bump a project to a new release:

```bash
cd path/to/RNAseq-toolkit        # the submodule
git fetch --tags
git checkout v0.2.0              # detached at the tag
cd -                             # back to parent repo
git add path/to/RNAseq-toolkit  # stage the new gitlink
git commit -m "Bump RNAseq-toolkit to v0.2.0"
```

**Commits & PRs:** concise imperative subjects (`Add pooled GSEA cache`); rationale in the
body when non-trivial. Group related changes; avoid formatting-only noise. PRs should
describe intent, key changes, and test evidence (`Rscript tests/test_volcano_plots.R`
output), note new dependencies, and attach representative PDFs/PNGs when plots change.

## Known Issues & Troubleshooting

- **ggplot2 4.0+:** use `color = "transparent"` instead of `color = NA` for shape-21 points
  (`NA` drops points as "missing values").
- **"Function not found":** you likely missed sourcing the specific script — there is no
  `library(RNAseqToolkit)`.
- **"Pathway name mismatch":** check `species` ("Mus musculus" vs "Homo sapiens") in `run_gsea`.
- **Theme errors:** source `custom_minimal_theme.R` before plotting if using toolkit defaults.

## Security & Configuration Tips

- Ensure R dependencies are installed in your library; add new imports explicitly and document them.
- Prefer relative paths from repo root; avoid writing outside project directories by default.
- Gate caching/output directories via parameters (`cache_dir`, `output_dir`) and create them if missing.
