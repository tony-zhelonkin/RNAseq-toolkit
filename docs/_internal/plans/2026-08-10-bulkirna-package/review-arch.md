# Architecture / deprecation-shim review — bulkiRNA @ fc016ff

Read-only review. Every behavioural claim below was executed in the container
(`scdock-r-dev:v0.5.11`, `devtools::load_all(".")`); scripts in `/tmp/arch/s1..s7.R`.

### 1. `gsea_barplot()` errors outright whenever no pathway passes `padj_cutoff`  [severity: high]
**R/deprecated-plot.R:~365** (`if (nrow(p$data) > 0L)`) with **R/gs-plot-utils.R:267** (`.gs_empty_plot()`)
— `gs_plot_bar()` applies `padj_max` as a hard filter and, when nothing survives,
returns `.gs_empty_plot(title)`: a bare `ggplot()` whose `$data` is an empty
`waiver`/`list`, so `nrow(p$data)` is `NULL` and `if (NULL > 0L)` throws
"argument is of length zero". The default `padj_cutoff = 0.05` makes this the
*normal* outcome for any run with no significant pathways.

**Evidence** (container):
```
### shim gsea_barplot on gs_result: ERR: argument is of length zero
### empty -> gsea_barplot : ERR: argument is of length zero
### gs_read -> gsea_barplot : ERR: argument is of length zero
Error in if (nrow(p$data) > 0L) { : argument is of length zero
3: gsea_barplot(res)
```
(the `gs_result` above had all `padj = 0.69`; `gsea_dotplot`, `gsea_dotplot_facet`,
`gsea_running_sum_plot` all returned OK on the same object)

**Failure scenario** — `res <- run_gsea(de); gsea_barplot(res)` on a dataset with
no FDR<0.05 pathway → hard error instead of an empty plot. Also
`gsea_barplot(empty_gsea_tibble())`.

**Suggested fix** — guard with `NROW(p$data)` (or `!is.null(p$data) && is.data.frame(p$data)`)
before the re-sort block in `gsea_barplot()`.

### 2. Six shims name unexported internals as their migration target  [severity: medium]
**R/deprecated-gs.R:110, 234, 279; R/deprecated-plot.R:~500, ~560** — `.Deprecated()`
targets and the `@description`/`Superseded by` prose point at
`the internal filter_by_size(...)`, `.gsdb_parse_mitoxplorer()`,
`.gsdb_human_to_mouse()`, `.gs_plot_all()`, `.gs_write_log()`, and `gs_db()` —
none of which a user can call. `gs_db()`, `gs_result()`, `gs_matrix()` and
`is_gs_db()` are all `@keywords internal` and absent from NAMESPACE.

**Evidence**:
```
$ for n in gs_db gs_result gs_matrix is_gs_db; do grep -c "^export($n)$" NAMESPACE; done
0 0 0 0
R/deprecated-gs.R:  .Deprecated("the internal filter_by_size(db, min_size, max_size) in R/gs-db.R")
R/deprecated-gs.R:  .Deprecated("the internal .gsdb_parse_mitoxplorer() in R/gsdb-rebuild.R")
R/deprecated-plot.R: .Deprecated(".gs_plot_all") / .Deprecated(".gs_write_log")
```
Corroborating: `R/gs-ops.R:177-182` builds a `gs_db` by hand with `structure(...)`
in its own `@examples` — the package's own docs work around the missing export.

**Failure scenario** — none at runtime; a user following the deprecation warning
gets `object '.gs_plot_all' not found` / must reach for `:::`.

**Suggested fix** — re-point each target at an exported successor
(`list_to_term2gene` → `gsdb_register()`, `parse_mitoxplorer` → `gsdb_from_file()`,
`plot_all_gsea_results` → `gs_plot_dot()`+`gs_save()`, `save_gsea_log` → `gs_save()`),
or export `gs_db()`/`gs_result()` as the public object vocabulary.

### 3. The `filter_by_size` shim's docs and warning describe a bug that was already fixed, and name a function that no longer exists  [severity: medium]
**R/deprecated-gs.R:75-110** — the roxygen block states at length that "the
collision itself ... is a package-level defect that this shim cannot fix", and
`.Deprecated()` sends users to "the internal `filter_by_size(db, ...)` in
R/gs-db.R". The integrator already renamed that internal to `.gs_filter_size()`.

**Evidence**:
```
$ grep -rn 'filter_by_size\|\.gs_filter_size' R/
R/gs-db.R:224:.gs_filter_size <- function(db, min_size = NULL, max_size = NULL, ...
R/deprecated-gs.R:110: .Deprecated("the internal filter_by_size(db, min_size, max_size) in R/gs-db.R")
# no top-level `filter_by_size` outside deprecated-gs.R
```
**Failure scenario** — none; documentation is now actively wrong (a reader will
look for a collision that isn't there, or call a name that isn't defined).

**Suggested fix** — rewrite the block to reference `.gs_filter_size()` (past
tense on the collision) and shorten the warning to a real successor.

### 4. `normalize_gsea_results()` tells users to call `gs_leading_edge()` with a `gs_db` — that signature does not exist  [severity: medium]
**R/deprecated-gs.R (`atlas_universe` warning) + R/gs-ops.R:187** —
the runtime warning is "Use `gs_leading_edge()` with the original `gs_db`", but
`gs_leading_edge(x, padj = NULL, top_n = NULL, unique_genes = FALSE)` has no `db`
argument, and `padj` is unvalidated, so following the advice produces a cryptic
coercion error.

**Evidence** (container):
```
### gs_leading_edge: ERR: 'list' object cannot be coerced to type 'double'
# from gs_leading_edge(res, db)  -> db lands in `padj`
$ grep -n 'gs_leading_edge <- function' R/gs-ops.R
187:gs_leading_edge <- function(x, padj = NULL, top_n = NULL, unique_genes = FALSE)
```
**Failure scenario** — `gs_leading_edge(res, db)` → `'list' object cannot be
coerced to type 'double'` rather than "unused/invalid argument".

**Suggested fix** — drop `gs_db` from the warning text, and validate `padj`/`top_n`
as numeric(1) in `gs_leading_edge()`.

### 5. `gs_test.gs_matrix()` requires `db` and then ignores it completely  [severity: medium]
**R/gs-test.R:107-120** — the generic is `gs_test(x, db, ...)` with no default for
`db`, and the `gs_matrix` method's body never references `db`; the `database`
column is taken from `attr(x, "database")`. Passing a mismatched db is silent.

**Evidence** (container):
```
### db ignored by gs_test.gs_matrix -> database col = td      # db was database="WRONG-DB"
### gs_test(sc, design=~group) without db: OK td              # db not actually required
```
**Failure scenario** — `gs_test(scores, gsdb_msigdb(...), design = ~group)` labels
every row with the *scoring* database, not the one passed; a user combining two
`gs_test()` calls by `rbind()` gets duplicate `database` keys with no warning.

**Suggested fix** — either error when `db` is supplied and its `database` attr
differs from `attr(x, "database")`, or give the method `db = NULL` and document
that scores already carry their database.

### 6. `gs_matrix` has no persistence path — the IO layer only handles `gs_result`  [severity: medium]
**R/gs-io.R:38,100** — `gs_write()`/`gs_read()` accept only `gs_result`; a `gs_matrix`
(the second first-class object in the contract, with `[`, `print`, `summary`,
`as_tibble` and a renderer) can only be saved by dropping to `saveRDS()`.

**Evidence** (container):
```
### gs_write(gs_matrix): ERR: `x` must be a `gs_result`, as returned by `gs_test()`; got 'gs_matrix'/'matrix'/'array'.
```
**Failure scenario** — none; a GSVA workflow simply has no supported way to
persist scores, and `gs_save()` covers only the plot.

**Suggested fix** — add a `gs_matrix` branch to `gs_write()`/`gs_read()` (scores
TSV + `sample_data` TSV + attrs), or state the asymmetry in `gs_write()`'s docs.

### 7. `test-namespace-hygiene.R` guards only the exact bug that was found  [severity: medium]
**tests/testthat/test-namespace-hygiene.R:19-38** — `top_level_functions()` matches
only `name <- function(...)` / `name = function(...)` at top level. It therefore
does **not** catch: (a) `assign("f", function(...))` or `` `f` <- `` via a string
LHS, (b) an S3 method whose generic is missing or is registered for a class with
no constructor, (c) a name exported *and* documented `@keywords internal`,
(d) `.Deprecated()`/`@seealso` targets that resolve to nothing, (e) functions
re-bound at load time. Items (c)/(d) are exactly findings 2 and 3, which this
test is silent about.

**Evidence** — read of the parser above; and both checks pass today while
`.Deprecated("the internal filter_by_size(...)")` names a nonexistent function.
Also `grep -Hn 'assign(' R/*.R` → no hits today, so (a) is a latent gap only.

**Failure scenario** — none; the class of bug the test was written for can recur
in a form it does not see.

**Suggested fix** — add two cheap checks: every `.Deprecated()` string argument
that looks like an identifier must exist in the namespace, and every
`S3method(g, c)` in NAMESPACE must have `g` available (base/imported/defined).

### 8. No shim test exercises the "nothing significant" path — which is why finding 1 shipped  [severity: medium]
**tests/testthat/test-deprecated-plot.R:90, 142** — both `gsea_barplot()` tests
pass a cutoff chosen to keep rows (`padj_cutoff = 1` and `padj_cutoff = 0.6`).
The default-cutoff / empty-result path is untested for all four plot shims.

**Evidence**:
```
90:  expect_warning(p <- gsea_barplot(res, padj_cutoff = 1, top_n = 3),
142:    p <- gsea_barplot(g, padj_cutoff = 0.6, top_n = 20),
```
`tests/testthat/test-gs-plot-bar.R:47` does test the empty case — but only on the
*new* renderer, never through the shim. Classic each-side-tests-its-own-half seam.

**Failure scenario** — see finding 1.

**Suggested fix** — one test per plot shim asserting a ggplot (not an error) when
every `padj` is above the cutoff, plus `gsea_barplot(empty_gsea_tibble())`.

### 9. `.Deprecated()` is called without `old =`, so the warning names the call site, not the function  [severity: low]
**R/deprecated-gs.R / R/deprecated-plot.R (all 20 shims)** — `.Deprecated(new)`
derives the old name from `sys.call()`, so indirect calls mis-report.

**Evidence** (container, all 20 shims invoked via `get(s)()`):
```
run_gsea      | deprecatedWarning,warning,condition | 'get(s)' is deprecated.
gsea_barplot  | deprecatedWarning,warning,condition | 'get(s)' is deprecated.
... (20/20 identical pattern)
```
Positive result worth recording: **all 20 shims do warn**, with a single
consistent condition class `deprecatedWarning`.

**Failure scenario** — legacy code using `do.call("run_gsea", args)` or
`lapply(dbs, run_gsea)` gets a warning naming the wrapper, so users cannot grep
their code for the deprecated name.

**Suggested fix** — pass `old` explicitly, e.g. `.Deprecated("gs_plot_bar", old = "gsea_barplot")`.

### 10. `rbind.gs_result()` silently merges incompatible `stat_type`s  [severity: low]
**R/gs-result.R (`rbind.gs_result`)**, consumed by
**R/deprecated-plot.R (`plot_all_gsea_results`)** — an fgsea result (`NES`) and an
ORA result (`log2_fold_enrichment`) rbind without complaint, and the renderers
then place both on one numeric axis.

**Evidence** (container):
```
### rbind fgsea+ora stat_types: NES+log2_fold_enrichment
### rbind warns? no warning
### mixed axis label: statistic | n rows: 5
```
Mitigating: the axis label degrades to the generic "statistic" rather than
claiming "NES", so the figure is vague rather than false.

**Failure scenario** — `plot_all_gsea_results(list(ora = ora_res, gsea = fgsea_res), ...)`
draws one dotplot mixing two unrelated statistics.

**Suggested fix** — warn in `rbind.gs_result()` when `unique(stat_type)` has
length > 1 (do not error: cross-database pooling is legitimate within one method).

---

## Seams enumerated and exercised (for the record)

Confirmed working in the container: `gs_ranks -> gs_test -> {gs_plot_dot, bar,
running, heatmap, gs_save, gs_write/gs_read}`; `gsdb_from_file -> gs_test`;
`gs_score -> gs_matrix -> gs_plot_heatmap / gs_save`; `gs_test.character (ORA) ->
gs_plot_dot/bar/.gs_plot_all`; `gs_filter/gs_top/gs_split` all preserve
`gs_result` class; `gs_read -> renderers` (`leading_edge` survives as a list);
`run_gsea()-shaped gs_result -> all four plot shims` (except finding 1);
`plot_all_gsea_results(list(gs_result))` with 1 and 2 elements (12 / 24 files);
`normalize_gsea_results -> gsea_dotplot / gsea_running_sum_plot`;
`empty_gsea_tibble -> gs_plot_dot/gs_plot_bar`.

Seams that error *by design*, with a clear message (not reported as findings):
ORA `gs_result -> gs_plot_running` / `gsea_running_sum_plot` ("`ranks` is
required"); `gs_test(ranks, <T2G data frame>)` (tells you to use `gsdb_*`);
`.dep_gsea_to_gs_result()` on a non-`gseaResult`.

Unexercised-anywhere pairs I could not break but which no test covers:
`gs_score -> gs_test.gs_matrix -> gs_plot_*` end-to-end (only the halves are
tested), and `gs_read()` output into any shim.

## Architecture verdict

**The four-layer contract holds at the level it was written for.** Mechanically:
`grep -l 'ggplot(|geom_|aes('` over `R/*.R` returns *only* the `gs-plot-*`,
`de-*`, `*theme*` files — no compute file plots; and `fgsea::`/`limma::`/`GSVA::`
appear in `gs-plot-running.R` solely as `plotEnrichmentData()` (curve geometry,
not testing) and in `de-md.R` as `limma::topTable` on a fit the caller supplied.
No duplicate top-level definition survives (`uniq -d` over all
`name <- function` LHS across the 34 files: empty). The 34-file split is the
right granularity: file names map onto layer prefixes (`gsdb-*`, `gs-*`,
`gs-plot-*`, `de-*`), the largest non-shim file is 544 lines, and I found no
function living in the wrong file.

**Weakest point: the layer boundaries are enforced in code but not in the public
surface.** The vocabulary of layers 1 and 3 — `gs_db()`, `gs_result()`,
`gs_matrix()`, `is_gs_db()` — is entirely unexported, so every "superseded by"
pointer into those layers is unfollowable (finding 2), the package's own examples
hand-build a `gs_db` with `structure()`, and `gs_matrix` has no IO (finding 6).
Secondary weakness: the `gs_matrix` branch of `gs_test()` takes a `db` it ignores
(finding 5) — the one place the layer signature and the layer behaviour disagree.

**Single highest-value structural change:** export the object vocabulary —
`gs_db()`, `gs_result()`, `gs_matrix()`, `is_gs_db()` (roxygen `@export`, drop
`@keywords internal`; NAMESPACE regenerates) — and then sweep every
`.Deprecated()` target and `Superseded by` line so it names an exported function.
That one move makes the contract described in DESCRIPTION checkable by a
newcomer, unblocks all six broken migration pointers, and removes the `:::` and
`structure()` workarounds already present in the package's own docs. It changes
no frozen name and no behaviour.
