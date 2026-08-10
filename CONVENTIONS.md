# bulkiRNA — internal conventions

Repo-internal. `.Rbuildignore`d, not shipped. Read this before writing a line of
`R/`; it is the contract the parallel module agents (B1–B5) share.

The design record lives outside this repo in
`/data1/users/antonz/pipeline/sciagent-rna/docs/` (`07_api-design.md` is the API
contract). Live state is `docs/_internal/plans/2026-08-10-bulkirna-package/00_INDEX.md`.

---

## 1. Ownership

| File | Owner |
|---|---|
| `DESCRIPTION` | Step A + integrator |
| `NAMESPACE` | **roxygen only** — never hand-edited, never committed by hand |
| `R/utils.R` | Step A + integrator |
| `R/gs-result.R` | Step A + integrator |
| `R/gs-matrix.R` | Step A + integrator |
| `R/bulkiRNA-package.R` | Step A + integrator |
| `tests/golden/`, `tests/fixtures/` | Step 0 — **do not touch** |
| everything else in `R/` and `tests/testthat/` | the module agent that created it |

Need a new dependency, a new shared helper, or a change to `gs_result`? **Put it
in your handback report.** Do not edit the file. Five agents editing one file is
the failure this table exists to prevent.

## 2. No new shared helpers

There is exactly one place shared helpers live: `R/utils.R`, and it currently
holds three things — `ensure_dir()`, `ensure_parent_dir()`, `%||%`. If you find
yourself writing a fourth, request it in the handback. This is the rule that
stops the `ensure_dir` collision from happening a second time.

Module-private helpers stay in the module file, prefixed `.` and documented
`@keywords internal`.

## 3. Roxygen

- Markdown roxygen (`Roxygen: list(markdown = TRUE)`). Use `[fn()]`, backticks,
  and `#'` throughout — no `\code{}`.
- **Every** function gets `@param` for each formal and an `@return`. Exported
  functions get a runnable `@examples`; internal ones do not need them.
- Non-exported functions: `@keywords internal` and **no `@export`**. They are
  still callable from other files in the package.
- Exported functions: `@export`. The export list is the ~30 names in
  `07_api-design.md` §7. Do not export anything not on it.
- S3 methods: `@export` on the method (roxygen emits `S3method()`), plus
  `@importFrom` for a generic owned by another package (e.g.
  `@importFrom tibble as_tibble`).
- Package imports go in the roxygen block of the function that needs them via
  `@importFrom pkg fn`. **`@import ggplot2` appears in exactly one place —
  `R/bulkiRNA-package.R`** — so no other file needs `ggplot2::` prefixes and no
  two files fight over it. Everything else uses `pkg::fn()` or `@importFrom`.
- `Suggests` packages are reached only behind a guard:

  ```r
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("`method = \"gsva\"` requires the GSVA package. Install it with ",
         "BiocManager::install(\"GSVA\").", call. = FALSE)
  }
  ```

## 4. Naming

- `snake_case` for functions and objects; `.` prefix for module-private helpers.
- Prefixes are the namespace: `gsdb_*` providers, `gs_*` compute/result/persist,
  `gs_plot_*` renderers, `de_*` differential expression, `theme_bulki`.
- File names mirror content: `R/gsdb-load.R`, `R/gs-test.R`, `R/gs-plot-dot.R`.
  Hyphens in file names, underscores in function names.
- Test files mirror the R file: `R/gs-test.R` → `tests/testthat/test-gs-test.R`.

## 5. Errors and messages

- `stop(..., call. = FALSE)` — the internal call stack is noise to a user.
- Say what was wrong, name the offending argument in backticks, and where useful
  say what to do: `` "`species` must be one of \"Mus musculus\", \"Homo sapiens\"; got \"mouse\"." ``
- Never `warning()` for something the user must act on, and never silently
  return an empty result where an error is the truth. An empty `gs_result` with
  zero rows is a valid answer; a `NULL` is not.
- No `print()`/`cat()` for progress in compute functions. Use `message()`, gated
  behind a `verbose` argument that defaults to `FALSE`.

## 6. Compute never plots; viz never computes

Structural, not aspirational:

- **Compute** (`gs_ranks`, `gs_test`, `gs_score`, `gs_leading_edge`) returns
  objects. No file I/O. No plotting. No `ggsave`.
- **Persist** (`gs_write`, `gs_read`, `gs_save`) is the *only* code that touches
  disk.
- **Render** (`gs_plot_*`, `de_*`) takes an object, returns a ggplot. No
  statistics computed inside a renderer beyond what the geom needs — if a plot
  needs a number the object does not have, the object's contract is wrong.
- `gs_save(plot, path)` writes the figure **and** its same-stem source table,
  PDF + PNG from one plot object. That craft rule is now the function's job, not
  an agent's memory.

## 7. The `stat_type` vocabulary

`gs_result$stat` is meaningless without `gs_result$stat_type`. The closed set,
from `gs_stat_types()`:

| `stat_type` | Method | Axis label |
|---|---|---|
| `NES` | fgsea | `NES` |
| `t` | GSVA → limma | `t statistic` |
| `log2_fold_enrichment` | ORA (`fgsea::fora`, `foldEnrichment`) | `log2 fold enrichment` |
| `signed_log10p` | reserved (camera/roast, no effect size) | `signed -log10 p` |

Rules:

1. Never write `stat` without `stat_type`.
2. Renderers get the axis label from `gs_stat_label(res)` — **never** a literal
   `"NES"`. A GSVA t-statistic mislabelled NES is exactly the class of quiet
   error this column exists to kill.
3. Adding a value means adding it to `gs_stat_types()` in `R/gs-result.R` — a
   Step A file, so it goes through the handback.
4. `direction` is `"up"` / `"down"` / `"ns"`, derived from `sign(stat)` by
   `gs_direction()`. Not a factor, not `"UP"`, not `1`/`-1`.

`gs_matrix` scores are *not* a `stat_type`; they carry `score_type` (`"gsva"`,
`"ssgsea"`, `"zscore"`, `"plage"`), and become `stat_type = "t"` only once
`gs_test()` runs limma over them.

## 8. Style

- Two spaces, no tabs. 80-character soft limit. Native pipe `|>`, not `%>%`.
- `TRUE`/`FALSE`, never `T`/`F`. Explicit `return()` only for early returns.
- ggplot2 is **4.0.3**: use `colour = "transparent"`, never `colour = NA` — `NA`
  drops shape-21 points as missing values.
- `seq_len()` / `seq_along()`, never `1:n`.
- Vectorized `vapply(..., FUN.VALUE=)` over `sapply()`.

## 9. Testing

- `tests/testthat/`, edition 3, one file per `R/` file.
- Test the contract, not the implementation: shape, types, error messages,
  boundary cases. Plot tests assert on `ggplot_build(p)$data`, never pixels.
- `tests/testthat/helper-gs.R` holds the small builders (`fake_gs_result()`,
  `fake_gs_matrix()`). Extend it rather than re-rolling a fixture per file.
- **`tests/golden/` and `tests/fixtures/` are off limits.** They are the
  behavioural baseline for the whole refactor; changing them destroys the only
  evidence that the refactor preserved results.
- Before handback: `devtools::load_all()` and your own tests pass.

## 10. Running R

There is no R on the host. Everything runs in a throwaway container, and both
`--user` and `HOME` are mandatory (`saveRDS` permissions; msigdbr's runtime
cache):

```bash
DK() { docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 "$@"; }

DK Rscript -e 'devtools::document()'
DK Rscript -e 'devtools::load_all(); testthat::test_local()'
DK Rscript tests/golden/verify_golden.R      # the Step C gate — must exit 0
```

R 4.5.3 · ggplot2 4.0.3 · fgsea 1.36.2 · msigdbr 26.1.0 (downloads MSigDB at
runtime to `$HOME/.cache/R/msigdbr`; needs network on a cold cache) · GSVA 2.4.9
(param-object API: `gsvaParam()` then `gsva(param)`).

## 11. Reference data

Processed databases ship in `inst/extdata/<dir>/processed/<Species>/*.rds` and
resolve **only** via `system.file("extdata", ..., package = "bulkiRNA")`. Never
walk up from a script's own location; `.resolve_toolkit_dir()` did that with
`sys.frame()$ofile` and worked only during sourcing.

Raw source files stay in `data/references/` in the git checkout and are
`.Rbuildignore`d. Anything requiring them (`rebuild = TRUE`) must fail with an
explicit "source checkout required" message rather than a missing-file error.

## 12. Deletion is the job

The surface goes from 123 definitions to ~30 exports. If a function has no
consumer in `data/used-functions.tsv` and no internal caller, **delete it** — do
not port it. Report what you deleted.
