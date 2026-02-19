# Repository Guidelines

## Project Structure & Module Organization
- `scripts/DE/` – Differential expression visuals (volcano, PCA 2D/3D); core plotting entry points like `plot_standard_volcano.R`, `plotPCA.R`, `plotPCA3d.R`.
- `scripts/GSEA/` – GSEA pipelines: processing (`GSEA_processing/*.R`) and plotting (`GSEA_plotting/*.R`, including pooled and cross-database dotplots).
- `scripts/custom_minimal_theme.R` plus `utils_plotting.R`/`famd_plotting_rnaseq.R` – shared styling helpers.
- `examples/example_analysis.R` – runnable reference for sourcing scripts and running analyses.
- `tests/` – `testthat`-driven volcano plot tests, PDFs for visual inspection, and `VOLCANO_FIX_SUMMARY.md` notes.
- `README.md` – overview, usage, and branch model; `MIGRATION.md` – upgrade notes.

## Build, Test, and Development Commands
- Run volcano test suite (auto + visual):  
  `Rscript tests/test_volcano_plots.R` (writes PDFs to `tests/output/` if enabled).
- Quick example workflow:  
  `Rscript examples/example_analysis.R` (safest from repo root; adjust input paths inside the script).
- Ad-hoc sourcing in R:  
  `source("scripts/GSEA/GSEA_processing/run_gsea_analysis.R")` then call helpers per README examples.

## Coding Style & Naming Conventions
- Language: R. Indent 2 spaces; avoid tabs. Prefer tidyverse style (pipes or |>), consistent spacing around `=` in args.
- Functions/objects: `snake_case`; exported helpers start lowercase (`run_gsea`, `create_volcano_plot`). File names mirror main function.
- Keep plotting side effects optional: return plot objects; gate file I/O behind explicit flags/paths.
- Add lightweight comments for non-obvious data wrangling or plotting logic; keep docstrings minimal but precise.

## Testing Guidelines
- Primary coverage: volcano alignment and edge cases via `testthat` in `tests/test_volcano_plots.R`. Expect PASS and PDFs that align dashed thresholds with color boundaries.
- When adding DE/GSEA plotting helpers, add a minimal reproducible fixture to `tests/` and extend the volcano test script or create a sibling `test_*.R`.
- Clean up large artifacts; keep generated outputs in `tests/output/` or a user-specified directory, not versioned.

## Commit & Pull Request Guidelines
- Branch model: `main` (stable), `dev` (integration), `dev-{project}` (project-specific). Develop on `dev-{project}`, merge to `dev`, then to `main` when validated.
- Commits: concise, imperative subject (`Add pooled GSEA cache`), include rationale in body if non-trivial. Group related changes; avoid formatting-only noise unless stated.
- PRs: describe intent, key changes, test evidence (`Rscript tests/test_volcano_plots.R` output), and note any new dependencies. Link issue/analysis context and, when plot changes matter, attach representative PDFs/PNGs.

## Security & Configuration Tips
- R dependencies: `clusterProfiler`, `msigdbr`, `limma`, `edgeR`, `ggplot2`, `ggrepel`, `pheatmap`, `plotly`, `testthat`; ensure they are installed in your R library. Add any new package imports explicitly and document them.
- File paths: prefer relative paths from repo root; avoid writing outside project directories by default. Gate caching/output directories via parameters (`cache_dir`, `output_dir`) and create them if missing.
