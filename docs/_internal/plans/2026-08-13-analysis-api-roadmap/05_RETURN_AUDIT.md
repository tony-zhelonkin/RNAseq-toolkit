# S5 return-type and error-message audit

**Date:** 2026-08-19 · **Status:** implemented; R gates pending
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

## Measured surface

`bulkirna_api()` contains 79 exports: 46 stable, 12 experimental, and 21
deprecated. The live surface is therefore **58**, as the S5 brief states. The
57-live count in older roadmap prose predates `gsdb_coresh()`.

The return audit below is derived from the non-deprecated rows of
`bulkirna_api()`. Forty-three exports have direct, source-tree fixture probes in
`test-return-audit.R`. Fifteen dependency-, network-, or external-tool-bound
exports have individual reasons and existing focused tests that assert their
return classes. This work was performed in a sandbox with no R or Docker, so
the table records the implemented and already-tested contracts; a fresh runtime
observation remains part of the unverified gates.

## Return types, one row per live export

“Holds” means the broad S5 category applies. “Exception” means the function's
established contract is narrower and should not be changed to make the slogan
literal. A `gs_result` is a tibble. A `gs_db` is intentionally a named list and
a `gs_matrix` is intentionally a numeric matrix under the frozen architecture
in `CONVENTIONS.md` §11a and the existing object contract.

| Export | Category | Actual return on the audit fixture | Finding |
|---|---|---|---|
| `annotate_genes()` | Tabular compute | `tbl_df` | Holds; optional annotation packages are covered by its focused test. |
| `build_dge()` | Domain constructor | `DGEList` | Exception: edgeR consumers require the established object. |
| `bulkirna_api()` | Metadata table | `tbl_df` | Holds. |
| `bulkirna_check_deps()` | Metadata table | `tbl_df` | Holds; it is not an analysis compute function. |
| `bulkirna_stochastic()` | Metadata table | `tbl_df` | Holds. |
| `coresh_chunks()` | Tabular compute | `tbl_df` | Holds. |
| `coresh_convergence()` | Tabular compute | `tbl_df` | Holds. |
| `coresh_loadings()` | Tabular compute | `tbl_df` | Holds. |
| `coresh_match()` | Tabular compute | `tbl_df` | Holds. |
| `coresh_search()` | Tabular compute | `tbl_df` | Holds; BiocParallel execution is covered by its focused tests. |
| `coresh_sets()` | Domain constructor | `gs_db` (named list) | Exception required by the `gs_db` provider contract. |
| `coresh_validate()` | Validation report | `tbl_df`, invisibly | Holds. |
| `de_bfc_plot()` | Renderer | `ggplot` | Holds. |
| `de_md_plot()` | Renderer | `ggplot` | Holds. |
| `de_pca()` | Renderer | `ggplot` | Holds; optional edgeR path is covered by its focused test. |
| `de_pca_3d()` | Interactive renderer | `plotly` | Exception: interactivity is the documented purpose and stable return. |
| `de_volcano()` | Renderer | `ggplot` | Holds. |
| `de_volcano_grid()` | Renderer | `patchwork`, inheriting `ggplot` | Holds; composition adds the patchwork subclass. |
| `ensure_dir()` | Writer utility | character path, invisibly | Holds; it is neither compute nor a renderer. |
| `entrez_to_gene()` | Identifier mapper | character vector | Exception: a vector preserves one output per input identifier. |
| `filter_confounder_genes()` | Gene preparation | character vector | Exception: filtering a vector returns a vector. |
| `format_pathway_name()` | Presentation helper | character vector | Exception: formatting preserves vector shape. |
| `gatom_de()` | Domain constructor | `gatom_de` / `data.frame` | Exception: this stable table subclass is GATOM's input boundary. |
| `gatom_download_refs()` | Writer | character paths, invisibly | Holds. |
| `gatom_genes()` | Domain extractor | character vector | Exception: the result is a set of gene symbols, not a result table. |
| `gatom_module()` | Domain compute | `igraph` | Exception: graph topology and attributes are the result. |
| `gatom_refs()` | Domain constructor | `gatom_refs` / list | Exception: it holds three heterogeneous reference objects. |
| `gatom_save_html()` | Writer | character path, invisibly | Holds. |
| `gene_to_entrez()` | Identifier mapper | integer vector | Exception: a vector preserves one mapping result per input symbol. |
| `gs_coregulation()` | Tabular compute | `gs_result` / `tbl_df` | Holds. |
| `gs_filter()` | Tabular compute | `gs_result` / `tbl_df` | Holds. |
| `gs_leading_edge()` | Domain extractor | named list by default | Exception: each pathway maps to a gene vector. |
| `gs_plot_bar()` | Renderer | `ggplot` | Holds. |
| `gs_plot_dot()` | Renderer | `ggplot` | Holds. |
| `gs_plot_heatmap()` | Renderer | `ggplot` | Holds. |
| `gs_plot_running()` | Renderer | `ggplot` | Holds. |
| `gs_ranks()` | Domain preparation | named numeric vector | Exception: fgsea consumes this exact representation. |
| `gs_read()` | Reader | `gs_result` / `tbl_df` | Holds. |
| `gs_save()` | Writer | character paths, invisibly | Holds. |
| `gs_score()` | Domain compute | `gs_matrix` / numeric matrix | Exception required by the `gs_matrix` score contract. |
| `gs_split()` | Domain extractor | named list of `gs_result` objects | Exception: the names encode split groups. |
| `gs_stat_types()` | Package metadata | named character vector | Exception: this is a closed vocabulary lookup. |
| `gs_test()` | Tabular compute | `gs_result` / `tbl_df` | Holds. |
| `gs_to_master()` | Tabular compute | `tbl_df` | Holds. |
| `gs_top()` | Tabular compute | `gs_result` / `tbl_df` | Holds. |
| `gs_validate_master()` | Validator | input tibble invisibly, or report tibble with `error = FALSE` | Holds; both branches remain tabular and are deliberate. |
| `gs_write()` | Writer | character directory with manifest attribute, invisibly | Defect fixed: it had returned the documented value visibly. |
| `gsdb_coresh()` | Domain constructor | `gs_db` (named list) | Exception required by the provider contract. |
| `gsdb_from_file()` | Domain constructor | `gs_db` (named list) | Exception required by the provider contract. |
| `gsdb_info()` | Metadata | list | Exception: citations and provenance are heterogeneous. |
| `gsdb_list()` | Metadata | `data.frame` | Exception: stable metadata API, not compute. |
| `gsdb_load()` | Domain constructor | `gs_db` for one key; named list for several | Exception required by the provider and multi-load contracts. |
| `gsdb_msigdb()` | Domain constructor | `gs_db` (named list) | Exception required by the provider contract; the network boundary is mocked in its focused test. |
| `gsdb_register()` | Domain constructor | `gs_db` (named list) | Exception required by the provider contract. |
| `read_counts_matrix()` | Reader | numeric matrix | Exception: the matrix is the natural count representation. |
| `read_metadata()` | Reader | `data.frame` | Exception: stable I/O return, not compute. |
| `theme_bulki()` | Renderer component | ggplot2 `theme` | Exception: a theme is composed into a plot and is not itself a plot. |
| `write_session_provenance()` | Writer | character path, invisibly | Holds. |

There are no live exported predicates, so the `logical(1)` category is empty.
The enforcement code still has an explicit category field; adding a predicate
without categorising it fails the live-export comparison.

## Error-message measurement

The source-tree audit parses every `stop()` call in `R/`. It uses exact
function-local call ordinals for its reasoned exceptions, so one broad reason
cannot silently excuse a new call. The baseline was measured before edits.

| Check across 311 `stop()` calls | Baseline | After S5 | Decision |
|---|---:|---:|---|
| `call. = FALSE` | 304 | 311 | Seven legacy/internal paths were fixed. |
| Argument named in backticks | 268 | 297 | Fourteen remaining calls have individual reasons because no caller argument exists. |
| Free of option-dependent `sQuote()` | 281 | 311 | Thirty error calls now use `encodeString(..., quote = "\"")`; all 39 uses in `R/`, including warnings and deferred validation text, were removed. |
| Explicit correction language | 231 | 242 | Conservative count: acceptable values/types, a function or argument to use, an install command, or a rebuild/refresh action must be stated literally. |

The last row is deliberately conservative. I reviewed the calls excluded by
that text rubric rather than treating a keyword as proof. Many exclusions are
still adequate because they show the required/allowed columns dynamically or
are private assertions about malformed package or third-party data. The worst
baseline validation messages were:

| Offender | Finding | Disposition |
|---|---|---|
| Deprecated `list_to_term2gene()` | Said only `geneset_list must be a named list`; omitted backticks, remedy, and `call. = FALSE`. | Now names `geneset_list`, says to name elements with pathway IDs, and suppresses the call. |
| Deprecated GMX and mitoXplorer parsers | Reported only “file not found”; omitted the `file` argument, remedy, and `call. = FALSE`. | Now name `file` and request an existing format-appropriate path. |
| Deprecated `run_gsea()` | Three validation branches leaked calls; the missing-rank message used locale-dependent single quotes and gave no available columns. | Calls are suppressed; observed values use stable quotes and available columns are listed. |
| CoReSh loading/set inputs | Several messages named a GSE in prose but not `gse_id`, `gpl`, `query`, or the indexed field. | Caller-controlled inputs are named and the matching GSE/GPL/query remedy is explicit. |
| Object validators | “Not a gs_matrix/gs_result” named neither `x` nor its constructor. | They now name `x` and the public constructor/compute function that creates it. |

The 14 no-argument error exceptions fall into four narrow groups: corrupt
deserialized CoReSh payload fields, malformed GESECA return objects, two
transfer failures caught and re-emitted as file-specific warnings, and private
package-schema/parser assertions. Each call has its own one-line reason in
`test-error-audit.R`.

## What changed

- `gs_write()` now wraps its existing directory-plus-manifest value in
  `invisible()`. No value, attribute, path layout, or signature changed.
- Seven `stop()` calls gained `call. = FALSE`.
- Caller-validation messages that named no argument now name it in backticks
  and, where the old text stopped at the problem, name the correction.
- Option-dependent `sQuote()` was removed from errors, warnings, and deferred
  validation details in favour of stable ASCII quoting.

No compute return class, frozen signature, exported name, generated manual,
fixture, golden baseline, RNG mutation site, or consumer project changed.

## Enforcement and failure conditions

`test-return-audit.R` derives the live names from `bulkirna_api()` and the
source `NAMESPACE`. Every name must appear exactly once in the category table
and exactly once among an executable fixture probe or a per-name exception.
The probe asserts the class and asserts invisibility for writers. It fails when
a new live export is uncategorised, a name is duplicated, a probed return class
changes, a writer becomes visible, or an exception lacks its own reason.

`test-error-audit.R` parses every function in `R/`. It fails when any new or
existing `stop()` omits `call. = FALSE`, uses `sQuote()`, or lacks an argument
name in backticks without matching one exact reasoned exception. Its allowlist
helper normalises `NULL` to `character(0)`, and both audit files exercise the
empty-allowlist case so the success shape remains stable.

Both tests skip under an installed-package test run because `R/` and the source
`NAMESPACE` are unavailable there. Their skip messages say that they are
source-tree-only and why. They execute under `devtools::test()`.

## Gates not run

No R command and no Docker command was executed in this sandbox. Therefore
`devtools::document()`, `devtools::test()`, the golden verifier, and
`rcmdcheck::rcmdcheck()` remain unverified. `NAMESPACE` and `man/` were not
edited.
