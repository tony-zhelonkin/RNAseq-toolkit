# Review — gene-set database provider layer (`R/gs-db.R`, `R/gsdb-*.R`, `R/format-pathway-name.R`)

Branch `feat/bulkirna-package`, HEAD fc016ff. Read-only; nothing in the repo was modified.
All R was run in the `scdock-r-dev:v0.5.11` container; scratch scripts in `/tmp/gsdbrev/`.

Baseline: the seven scoped test files all pass
(`testthat::test_dir("tests/testthat", filter="gs-db|gsdb-|format-pathway-name|namespace")`
→ 0 failures, 0 warnings).

---

### 1. A GMT file with uniform column counts and a non-`.gmt` extension is silently parsed as GMX, producing garbage sets  [severity: medium-high]

**R/gsdb-file.R:118-129** (`.gsdb_sniff_format`) — the GMX heuristic only requires
`widths[1] == widths[2] && widths[1] >= 2 && all(widths[-c(1,2)] <= widths[1])`.
Any GMT whose first six lines happen to have the same number of tab fields (very common
for small/hand-made or programmatically emitted GMTs) satisfies it. The file is then run
through `.gsdb_parse_gmx()`, which takes row 2 as the set-name row and row 1 as
descriptions — so set ids become gene symbols and the description field becomes a gene.
No warning is emitted.

**Evidence** — ran in the container:
```r
l <- c("S1\tdesc\tg1\tg2","S2\tdesc\tg3\tg4","S3\tdesc\tg5\tg6")
f <- tempfile(fileext=".txt"); writeLines(l,f)
bulkiRNA:::.gsdb_sniff_format(f,l)   #> "gmx"
names(gsdb_from_file(f))             #> "S2" "desc" "g3" "g4"
```
Expected `c("S1","S2","S3")`. The old code had no sniffing at all — `parse_gmt()` and
`parse_gmx()` were separate explicit entry points — so this failure mode is new.

**Failure scenario** — a user calls `gsdb_from_file("my_sets.txt")` on a legitimate GMT
where every set has the same number of genes. They get four bogus "pathways" named after
the second record's fields, and every downstream `gs_test()` result is meaningless.

**Suggested fix** — tighten the heuristic (require row 1 to contain no value that also
appears in rows 3+, or require ≥3 data rows *and* that row 2's fields look like ids), and
at minimum `message()`/`warning()` which format was inferred when the extension did not
decide. A `format = c("auto","gmt","gmx")` argument would let callers opt out entirely.

---

### 2. `pathway_descriptions` is written by `gsdb_msigdb()`, destroyed by any subsetting, and never read anywhere  [severity: medium]

**R/gsdb-msigdb.R:86-90** and **R/gs-db.R:203-214** (`[.gs_db`) — the MSigDB provider goes
out of its way to attach `pathway_descriptions` *after* the size filter (with a comment
explaining that `[.gs_db` would otherwise misalign it), but `[.gs_db` reconstructs via
`gs_db()` passing only the five known attributes, so the attribute is dropped by the very
next subset. Grep shows no reader for it anywhere in `R/` or `tests/`.

**Evidence** — container run, plus grep:
```r
db <- bulkiRNA:::gs_db(list(A="a",B="b"), database="d", species="Mus musculus")
attr(db,"pathway_descriptions") <- c(A="desc A", B="desc B")
!is.null(attr(db,"pathway_descriptions"))        #> TRUE
!is.null(attr(db["A"],"pathway_descriptions"))   #> FALSE
```
```
$ grep -rn "pathway_descriptions" R/ tests/
R/gsdb-msigdb.R:89:    attr(db, "pathway_descriptions") <- descriptions[names(db)]
```

**Failure scenario** — a user does `db <- gsdb_msigdb("Mus musculus","H"); db <- db[1:10]`
and the descriptions vanish with no error. More generally: the `gs_db` contract documented
at gs-db.R:9-19 lists five attributes; a sixth exists in practice and is not part of the
subsetting contract, which is exactly the class of drift `[.gs_db` was written to prevent.

**Suggested fix** — either make `pathway_descriptions` a first-class optional argument of
`gs_db()` that `[.gs_db` carries through (and document it in the contract block), or drop
it until something reads it. Half-supported attributes are worse than absent ones.

---

### 3. `gsdb_list()`'s yaml-less fallback returns a different column order — and a different set of rows — than its documented shape  [severity: medium-low]

**R/gsdb-load.R:87-99** — the fallback branch `return()`s early with columns
`database, name, bundled, description, species`, while the yaml branch ends with an
explicit `out[c("database","name","bundled","species","description")]`. The roxygen
`@return` documents the latter order. The fallback also silently omits non-bundled
entries (`gatom`) and carries a stray `names` attribute on the `name` column
(`vapply(reg, [[, ..., "label")` keeps registry names).

**Evidence** — container run reproducing both branches:
```
== yaml branch ==     "database" "name" "bundled" "species" "description"
== fallback branch == "database" "name" "bundled" "description" "species"
```

**Failure scenario** — any caller that positionally indexes (`gsdb_list()[[4]]`) or
`rbind()`s results across machines gets `species` in one environment and `description` in
another, depending only on whether `yaml` happens to be installed. Silent, environment-
dependent wrong column.

**Suggested fix** — build both branches into the same `out[c(...)]` reordering statement
(one shared tail), and `unname()` the `name` column.

---

### 4. `gsdb_from_file()` has no `database_label` argument, so the rebuild driver mutates the attribute behind the constructor's back  [severity: medium, maintainability]

**R/gsdb-file.R:28-34** (formals) and **R/gsdb-rebuild.R:60-63**. Every other provider
(`gsdb_register()`, `gsdb_load()`, `.gsdb_parse_mitoxplorer()`, `.gsdb_parse_transportdb()`)
sets `database_label` through `gs_db()`. `gsdb_from_file()` instead hard-wires
`database_label = database %||% basename(path)`, so when `database` is supplied the display
label *is* the machine key — the one thing gs-db.R:14-19 says the two attributes must not
be. `.gsdb_rebuild()` then works around it:
```r
db <- gsdb_from_file(gmx, database = "mitopathways", ...)
attr(db, "database_label") <- "MitoPathways 3.0"   # bypasses the constructor
```

**Evidence** — `names(formals(gsdb_from_file))` in the container:
`"path" "database" "species" "prefix" "min_size" "max_size" "verbose"` — no
`database_label`. Source quoted above.

**Failure scenario** — none at runtime today; maintainability only. But it is a live trap:
the next attribute added to `gs_db()` will be initialised by the constructor and then
silently clobbered/skipped by this direct `attr<-`, and any user calling
`gsdb_from_file(f, database = "my_key")` gets `"my_key"` as their plot's database label.

**Suggested fix** — add `database_label = NULL` to `gsdb_from_file()`, default it to
`basename(path)` independent of `database`, and delete the `attr<-` in `.gsdb_rebuild()`.

---

### 5. Subsetting a `gs_db` by an unknown set name errors with a message about the constructor's internals  [severity: low]

**R/gs-db.R:203-214** — `db["nope"]` produces `unclass(x)["nope"]`, a one-element list with
an `NA` name, which `gs_db()` rejects with `"`sets` must have non-missing, non-empty
names."` The message names the wrong thing: the caller supplied an index, not `sets`.

**Evidence** — container run:
```r
db <- bulkiRNA:::gs_db(list(A="a"), database="d", species="Mus musculus")
try(db["nope"])
#> Error : `sets` must have non-missing, non-empty names.
```

**Failure scenario** — a user typos a pathway id in `db[my_ids]` and gets an error that
sends them looking at the database file rather than at their id vector.

**Suggested fix** — in `[.gs_db`, validate character `i` against `names(x)` first and error
with the offending ids (`"`i` selects sets not in this database: \"nope\""`).

---

### 6. `.gs_wrap_label()` is a plot-layout helper living in the name-formatting file  [severity: low]

**R/format-pathway-name.R:294-335** — the file's stated job is turning database ids into
human labels. `.gs_wrap_label()` inserts `\n` for axis rendering and is called only from
the renderer layer.

**Evidence** — grep for callers:
```
R/gs-plot-utils.R:169:  df$label <- .gs_wrap_label(
R/gs-plot-heatmap.R:157:      .gs_wrap_label(format_pathway_name(pretty, ...),
```
No provider or compute-layer caller.

**Failure scenario** — none; maintainability only. It is the one thing in this file that a
renderer-only change would touch, which makes the file do two jobs.

**Suggested fix** — move it to `R/gs-plot-utils.R` next to its only callers (and move the
two `.gs_wrap_label` tests alongside).

---

### 7. Rebuild writes a different `source` provenance string than the old builder  [severity: low]

**R/gsdb-rebuild.R:99** — `legacy$source <- attr(out[[key]], "database_label")` yields
`"MitoPathways 3.0"` / `"Unified Mitochondrial Pathways"` / `"TransportDB 2.0"`.

**Evidence** — old vs new source. Old `build_reference_databases.R`:
```r
mp_result$source <- "MitoPathways3.0"
... source = "Unified_Mitochondrial_Pathways",
tdb_result$source <- "TransportDB2.0"
```
New: the display label, with spaces.

**Failure scenario** — nothing in the package reads `source` (`gsdb_load()` only touches
`raw$T2G`), so no user-visible breakage; but a parent project that recorded provenance by
string-matching `readRDS(...)$source == "TransportDB2.0"` would silently stop matching
after a rebuild.

**Suggested fix** — either write the old machine string (add a `source_id` to the registry)
or state in the roxygen that `source` is now the display label, so the change is deliberate
and documented rather than incidental.

---

### 8. Unified-mitochondria merge is first-wins on colliding ids, where the old builder unioned their genes  [severity: low]

**R/gsdb-rebuild.R:76-85** — `merged[!duplicated(names(merged))]` drops the later of two
same-named sets. The old builder `rbind()`ed the T2G frames and `split()` by `gs_name` in
`deduplicate_genesets()`, which *unions* the two sets' genes.

**Evidence** — old-vs-new quote. Old:
```r
unified_T2G <- do.call(rbind, lapply(components, `[[`, "T2G"))
gs_list <- split(T2G$gene_symbol, T2G$gs_name)      # same name -> merged genes
```
New:
```r
merged <- c(unclass(out$mitopathways), unclass(out$mitoxplorer))
unified <- gs_db(merged[!duplicated(names(merged))], ...)
```

**Failure scenario** — unreachable with today's inputs, because the two components are
prefixed `MITOPATHWAYS_` and `MITOXPLORER_` and cannot collide. It becomes reachable the
moment a third component (or an unprefixed one) is added to the merge.

**Suggested fix** — leave the behaviour but say why it is safe in a comment ("component ids
are prefix-disjoint, so `!duplicated` never fires; revisit if a component is added"), or
union explicitly.

---

## Categories with nothing to report

- **Behavioural divergence from the old code (the priority-1 category)** — beyond items 7
  and 8 above, **nothing found**. I specifically checked the rebuild pipeline, which is
  where the ordering of operations changed: the old builder size-filtered *inside*
  `parse_gmx()` (on human symbols, before `convert_human_to_mouse()`) and again afterwards
  via `filter_by_size()`; the new `.gsdb_rebuild()` filters only after
  `.gsdb_human_to_mouse()`. I ran the real raw file in the container to see whether that
  changes the output:
  ```
  total sets: 149   sets with >500 human genes: 0   with <5: 30
  shipped mouse mitopathways sets: 119            # == 149 - 30
  ```
  Because no set exceeds `max_size` before mapping, and mapping only shrinks sets, the
  dropped-pre-filter is a no-op on the shipped data — the two pipelines agree. (It would
  diverge for a database containing sets of >500 genes; worth a one-line comment, not a
  finding.) `format_pathway_name()` is a faithful port of the old function plus an added
  `NA` guard; the `"^GO "` entry in `.gs_name_prefixes()` looked like a double-anchor bug
  but I verified it still strips correctly (`format_pathway_name("GO_APOPTOTIC_PROCESS")`
  → `"Apoptotic Process"`).
- **Namespace collisions / shadowing** — **nothing found**.
  `grep -rhoP '^`?[.A-Za-z][A-Za-z0-9._]*`?(?= <- function)' R/*.R | sort | uniq -d`
  returns empty: no duplicate top-level bindings anywhere in `R/`. The `.gs_filter_size()`
  rename holds, `test-namespace-hygiene.R` passes, and all provider internals correctly
  carry the leading dot.
- **Layer violations** — **nothing found**. No `gsdb_*` file references ggplot2, a
  `gs_result`, or a test/score function; `format_pathway_name()` is a pure string utility
  legitimately shared by the provider and renderer layers.
- **Dead code** — `.gsdb_rebuild()` has no non-test caller, but it is explicitly documented
  as a maintainer tool, so that is intentional, not dead.
