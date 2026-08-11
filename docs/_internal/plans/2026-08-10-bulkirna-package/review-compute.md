# Compute-layer review — bulkiRNA @ fc016ff (read-only)

Scope: R/gs-test.R, gs-score.R, gs-ranks.R, gs-ops.R, gs-result.R, gs-matrix.R + their tests.
All claims below were executed in the container (`scdock-r-dev:v0.5.11`, `devtools::load_all("/pkg")`).

### 1. `gs_test()` silently clobbers the caller's global RNG stream  [severity: medium]
**R/gs-test.R:236** (`.gs_fgsea`: `if (!is.null(seed)) set.seed(seed)`)
`set.seed()` is called on the global stream and never restored. Any random draw after a
`gs_test()` call is taken from the seed-123 stream, not the user's. Because the seed is
re-applied per database, this also happens once per element of a `db` list.

**Evidence**
```r
set.seed(999); before <- runif(3)
set.seed(999); gs_test(ranks, db, min_size=5, max_size=50); after <- runif(3)
before -> 0.38907138 0.58306072 0.09466569
after  -> 0.0455565  0.5281055  0.8924190     # identical: FALSE
```
**Failure scenario** A script does `set.seed(42)` at the top, runs GSEA, then bootstraps or
permutes; the bootstrap is silently drawn from seed 123 and is *not* reproducible from the
script's own seed. Two scripts with different top-level seeds produce identical "random" draws.
**Suggested fix** Save and restore: `old <- .Random.seed; on.exit(assign(".Random.seed", old, .GlobalEnv))`,
or use `withr::with_seed(seed, fgsea::fgseaMultilevel(...))`. (The old `run_gsea()` had the same
leak, so this is inherited, not introduced — but it is a package now, and packages should not
mutate the user's RNG.)

### 2. A partially named `db` list produces an empty `database` label  [severity: medium]
**R/gs-test.R:135** — `nms <- names(db) %||% vapply(...)`. `%||%` only fires when `names(db)`
is `NULL`. For a partially named list `names()` returns `""` for the unnamed slots, so the
fallback never runs, and the inner `nm %||% attr(d, "database")` is likewise dead code because
`""` is not `NULL`.

**Evidence**
```r
gs_test(ranks, list(fake_gs_db(database="alpha"), beta = fake_gs_db(database="betaDB")), ...)
unique(res$database)  ->  ""  "beta"
```
**Failure scenario** `gs_test(ranks, list(gsdb_msigdb("H"), kegg = ...))` → half the rows carry
`database = ""`; `gs_split()`/`gs_top(per="database")`/facetted plots then group under a blank
label instead of "Hallmark".
**Suggested fix** Replace with a per-element fallback:
`nms <- names(db); nms[is.null(nms) | !nzchar(nms)] <- vapply(db[...], function(d) attr(d,"database"), "")`.

### 3. `gs_test()` on a `gs_matrix` hard-codes `method = "gsva"` regardless of the scoring method  [severity: medium]
**R/gs-test.R:118** — `gs_result(out, ..., method = "gsva", stat_type = "t")`, ignoring
`attr(x, "method")` / `attr(x, "score_type")` that `gs_score()` faithfully recorded.

**Evidence**
```r
gm <- gs_matrix(m, database="d", method="ssgsea", score_type="ssgsea", sample_data=meta)
r  <- gs_test(gm, design = ~0+group, contrast = "groupKO-groupWT")
attr(gm,"method") -> "ssgsea"      unique(r$method) -> "gsva"
```
**Failure scenario** An ssgsea/zscore/plage run is recorded and exported (`gs_save`, methods
section of a paper) as GSVA. Provenance in the result object is wrong for three of the four
`gs_score()` methods.
**Suggested fix** `method = attr(x, "method") %||% "gsva"`, and widen the documented `method`
vocabulary in `gs_result-class` accordingly (or keep a `score_type` column).

### 4. Type and vocabulary validation is bypassed by `[` and every dplyr verb  [severity: medium]
**R/gs-result.R:349 and :373** — both `[.gs_result` and `dplyr_reconstruct.gs_result` call
`new_gs_result()` (the unchecked constructor) after testing only *column presence*.
`validate_gs_result()` is therefore only ever reached from `gs_result()` itself.

**Evidence**
```r
bad <- r |> dplyr::mutate(direction = "Up", stat_type = "bogus", padj = "oops")
class(bad)[1]                      -> "gs_result"
validate_gs_result(bad)            -> Error: column(s) have the wrong type: padj
gs_filter(bad, direction = "up")   -> 0 rows      # silently, no error
```
**Failure scenario** A user (or a future renderer) writes `mutate(res, padj = p.adjust(p_value))`
returning a character/factor column, or restores the old `"Up"`/`"Down"` capitalisation. The
object still claims to be a `gs_result`; `gs_filter(direction = "up")` then silently returns zero
rows and the plot comes out empty with no error. This is exactly the class of silent-wrong the
contract exists to prevent.
**Suggested fix** Call `validate_gs_result()` (not just the column-presence check) in
`dplyr_reconstruct.gs_result`, downgrading to a plain tibble rather than erroring when it fails.

### 5. Subsetting a `gs_matrix` to zero rows or zero columns is a hard error  [severity: medium]
**R/gs-matrix.R:262–283** — `[.gs_matrix` funnels through the full `gs_matrix()` constructor,
which rejects a matrix whose `rownames()`/`colnames()` are `character(0)` via the
`is.null(rownames(x))` check (it also rejects duplicated indices).

**Evidence**
```r
gm[integer(0), ]      -> Error: `x` must have both row names (pathways) and column names (samples).
gm[, character(0)]    -> Error: (same)
gm[c(1,1), ]          -> Error: `x` has duplicated pathway row names.
```
Compare `gs_result`, where "an empty result is a valid answer; `NULL` is not"
(`.gs_empty_core`, R/gs-test.R:373) — the two object contracts disagree on emptiness.
**Failure scenario** `gm[rownames(gm) %in% keep, ]` where `keep` matches nothing (a pathway
filter with no hits) aborts the script instead of yielding an empty matrix. Selecting by an
unmatched name (`gm[c("p1","nope"), ]`) errors with a raw `subscript out of bounds` rather than
a package-level message.
**Suggested fix** Allow zero-extent matrices in the constructor (test `is.null()`, not
emptiness), and drop the duplicate check for the subset path (or use `new_gs_matrix()` there,
since `[` already produces a well-formed object).

### 6. An empty fgsea/ORA result loses the optional columns, so `rbind` and `gs_leading_edge` behave differently for it  [severity: low]
**R/gs-test.R:369–380** — `.gs_empty_core()` returns only the 12 core columns; no
`leading_edge`, `es`, `fold_enrichment`, `overlap`.

**Evidence**
```r
e <- gs_test(ranks, fake_gs_db(list(TOO_BIG = paste0("G",1:100))), min_size=5, max_size=10)
names(e) -> pathway_id … padj      # 12 cols, no leading_edge
gs_test(ranks, list(a = full_db, b = empty_db), min_size=5, max_size=10)  # -> 0 x 12, no leading_edge
rbind(e, nonempty)  # OK: bind_rows keeps the list column (verified, class "list")
```
**Failure scenario** A loop over contrasts where one contrast yields no pathways: the pooled
result is fine, but a per-contrast branch calling `gs_leading_edge()` hits the
"carries no `leading_edge` column" error instead of getting an empty list — an empty answer is
reported as a usage error.
**Suggested fix** Give `.gs_empty_core()` an optional-column argument (or have `.gs_fgsea`/`.gs_ora`
attach a zero-length `leading_edge` list before returning), so the shape is method-stable.

### 7. The two adapters adjust p-values over different sets, and it is undocumented  [severity: low]
**R/gs-test.R:250 vs :299** — the fgsea path recomputes `stats::p.adjust(..., "BH")` over the
rows it returns; the ORA path takes `fora()`'s `padj`, which was computed over all size-passing
pathways *including* the zero-overlap ones that line 292 then drops.

**Evidence**
```r
# 5 sets pass minSize; 3 have zero overlap and are dropped
   pathway_id      p_value         padj overlap  BH_over_returned_rows
1           A 3.292663e-23 1.646331e-22      20           6.585325e-23
2           B 8.398315e-02 2.099579e-01       5           8.398315e-02
# padj/p_value ratio = 5 (all tested sets), not 2 (returned rows)
```
**Failure scenario** None numerically wrong — arguably ORA's scope is the *better* one — but a
reader cannot reproduce `padj` from the returned table, and `summary.gs_result()`'s `n_pathways`
understates the multiplicity that produced it.
**Suggested fix** State the adjustment scope for each method in `gs_result-class`, and consider
keeping the zero-overlap rows (with `stat = NA`) rather than dropping them.

---

## Categories with nothing to report

- **Old-vs-new equivalence, ranking**: `gs_ranks()` reproduces `run_gsea()`'s preamble faithfully
  (`git show ff80de2^:scripts/GSEA/GSEA_processing/run_gsea.R`): NA-drop with warning, decreasing
  sort, infinite-drop with warning, `minSize=10`/`maxSize=500`/`eps=0`/`nPermSimple=1e5`/`seed=123`
  all match `clusterProfiler::GSEA`'s forwarded defaults. NaN is covered by `is.na()`; the
  sort-before-Inf-drop order is harmless. Verified by execution.
- **Constructor bypass**: no function in `R/` fabricates a `gs_result` or `gs_matrix` with a bare
  `structure(...)`/`class<-`. `grep -rn "structure(" R/` hits are all roxygen examples building the
  frozen `gs_db` shape, plus `gatom.R` and the legitimate `new_gs_matrix()`. The `gs_db`
  examples do bypass a constructor, but that is B1's layer, not mine.
- **Namespace hygiene / layer separation**: no ggplot2 call in any of the six compute/object
  files; no `fgsea::`/`GSVA::`/`limma::`/`p.adjust` in `gs-plot-*.R` other than
  `fgsea::plotEnrichmentData()` in `gs-plot-running.R`, which is fgsea's own plotting helper.
  No duplicate top-level definitions (test-namespace-hygiene.R covers this and I found no
  adjacent shadowing in scope).
- **`gs_ops` / attribute preservation**: `gs_filter`/`gs_top`/`gs_split`/`gs_leading_edge` all
  gate on `.gs_check_result()`, and class survives `[`, `rbind` and `dplyr` verbs (finding 4 is
  that it survives *too* readily). `gs_top(by_direction=TRUE)` and NA-padj push-to-bottom both
  behave as documented — verified by execution.
