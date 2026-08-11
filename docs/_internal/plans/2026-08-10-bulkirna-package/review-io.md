# Review: IO / saving / shared utilities / gatom — bulkiRNA @ fc016ff

Scope: R/gs-io.R, R/gs-save.R, R/data-io.R, R/utils.R, R/gatom.R, R/bulkiRNA-package.R,
their tests, DESCRIPTION, .Rbuildignore, inst/extdata reachability.

All five scoped test files pass: test-utils (4/4), test-gs-io (5/5), test-gs-save (10/10),
test-data-io (11/11) in `scdock-r-dev:v0.5.11`; test-gatom 20/20 with **zero skips** in
`scbio-singleuser:v1.7.1` (the standard image skips the five real-pipeline tests).

**Namespace hygiene: nothing found.** A strict start-of-line scan of all 34 files in `R/`
found no duplicated top-level definition:
`grep -rnoE "^\`?[A-Za-z._%][A-Za-z0-9._%<-]*\`? <- function" R/ | ... | sort | uniq -d` → empty.
`R/utils.R` holds exactly three things (`ensure_dir`, `ensure_parent_dir`, `%||%`), the
`ensure_dir`/`ensure_parent_dir` split is documented as the fix for the old collision
(old code had `ensure_dir` twice: `scripts/utils_plotting.R:16` and
`scripts/GSEA/GSEA_processing/build_reference_databases.R:49`), and no layer file
redefines them. **inst/extdata: nothing found** — every access goes through
`.gsdb_extdata()` → `system.file("extdata", ..., package = "bulkiRNA")`
(R/gsdb-load.R:43); `grep -rn '"data/|"inst/' R/` returns nothing.
**.Rbuildignore: nothing found** — `docs`, `examples`, `data`, `tests/golden`,
`tests/fixtures`, `tests/output`, `tests/.*\.pdf`, all agent md files are covered.

---

### 1. Raw featureCounts output cannot be read at all — the `#` header line breaks it  [severity: high]
**R/data-io.R:12 (`.io_guess_sep`), :63 (`read_counts_matrix`)** — the docstring
advertises "featureCounts, raw -- `Geneid, Chr, Start, End, Strand, Length` then samples
from column 7 on", but `featureCounts` always writes a `# Program:featureCounts v...;
Command:...` line first. `.io_guess_sep()` sniffs **line 1 only**, so it sniffs the
comment; `utils::read.delim()` has `comment.char = ""`, so the comment becomes the header.
The old implementation used `data.table::fread()`, which auto-detects and skips it.

**Evidence** — ran in `scdock-r-dev:v0.5.11` on a two-sample raw featureCounts file whose
first line is the standard `# Program:featureCounts v2.0.1; Command:"featureCounts" "-a"
"g.gtf" ... "s1.bam" "s2.bam"`:
```
guess sep: ;
Error : Unsupported counts file format: `/tmp/fcx/counts.txt` has no `Geneid` or `gene_id` column.
```
Old vs new: `git show ff80de2^:scripts/General/io_helpers.R` → `dt <- fread(fp)`;
new → `utils::read.delim(path, sep = .io_guess_sep(path), ...)`.
`tests/testthat/test-data-io.R:44` ("raw featureCounts metadata columns are dropped")
passes only because its fixture is written by `write.table()` and has no comment line — it
tests the shape, not the file format.

**Failure scenario** — `read_counts_matrix("counts.txt")` on unmodified featureCounts
output → hard error claiming the file has no `Geneid` column, which it does.

**Suggested fix** — skip leading `#` lines: read the first non-`#` line in
`.io_guess_sep()`, and pass `comment.char = "#"` (or `skip = <n>`) to `read.delim()`. Add a
test whose fixture includes the real featureCounts comment line.

### 2. `ensure_dir()` reports success when directory creation failed  [severity: high]
**R/utils.R:22-24** — `dir.create(p, recursive = TRUE, showWarnings = FALSE)` suppresses
the only signal, the return value is discarded, and nothing checks `dir.exists(p)`
afterwards. Every writer in the package (`gs_write`, `gs_save`, `.gs_write_log`,
`write_session_provenance`, `download_gatom_references`, `gatom_save_html`) trusts it, so a
failed `mkdir` surfaces later as an unrelated error from the *writing* library.

**Evidence** — ran in `scdock-r-dev:v0.5.11` with a mode-0500 parent:
```
ensure_dir returns: /tmp/rotest/figures  dir.exists: FALSE
gs_save error: Error in `ggsave()`: ! Cannot find directory '/tmp/rotest/figures'.
               ℹ Please supply an existing directory or use `create.dir = TRUE`.
```
`ensure_dir()` returned its path as if it had worked.

**Failure scenario** — `gs_save(p, "/read/only/mount/fig/dot")` → the user is told to
"supply an existing directory", pointing at ggplot2 rather than at the permission problem;
with `formats` partly written, output is half-complete.

**Suggested fix** — capture `dir.create()`'s return and `stop()` with the offending path
and `file.access()` reason when `!dir.exists(p)` afterwards.

### 3. `gs_write()` never clears stale output, and `gs_read()`'s fallback silently row-binds it  [severity: medium]
**R/gs-io.R:60-74 (write), :137-146 (read)** — `gs_write()` only adds files; a re-run with
a different contrast set leaves the old `by_contrast/<contrast>/` directories in place. The
`gs_read()` fallback globs `by_contrast/*/<name>_*.tsv`, so it picks the leftovers up with
no warning. (With `overview = TRUE` the pooled file masks the problem, so this bites exactly
the `overview = FALSE` workflow `gs_write()` supports.)

**Evidence** — ran in `scdock-r-dev:v0.5.11`: wrote contrasts `A_vs_B, C_vs_D` to a dir, then
re-wrote only `A_vs_B` into the same dir with `overview = FALSE`:
```
files on disk: "A_vs_B/gsea_demo.tsv" "C_vs_D/gsea_demo.tsv"
gs_read rows: 2 contrasts: A_vs_B,C_vs_D
```
The second write's own `attr(,"files")` lists only `A_vs_B` — so the returned paths and what
`gs_read()` returns disagree.

**Failure scenario** — analyst drops a contrast, re-runs the stage, reads the tables back for
a figure → a contrast that is no longer in the analysis appears in the plot, dated from the
previous run.

**Suggested fix** — either add a `prune`/`clean` argument that `unlink()`s `by_contrast/`
before writing, or make `gs_read()` warn when the globbed file set does not match a manifest
written alongside it.

### 4. A truncated download is cached forever as "already present"  [severity: medium]
**R/gsdb-gatom.R:65-71** — the skip branch is `file.exists(dest_file) && !isTRUE(overwrite)`
→ message `[skip] ... (exists; ...)` and the path is added to `downloaded`. There is no size
or integrity check, and `download.file()` writes in place rather than to a temp file, so an
interrupted transfer leaves a partial `.rds` that all later runs report as fine. The
docstring's promise that "a partial run is visible rather than silent" holds only for the
current call.

**Evidence** — code quote above (no `file.info(dest_file)$size` guard anywhere in the
function; contrast with the success branch at :83 which *does* read `$size` purely for the
message). Downstream, `gatom_refs()` (R/gatom.R:138-140) calls `readRDS()` on it unguarded.

**Failure scenario** — network drops during the 24 MB `network.kegg.rds` fetch → every
subsequent `download_gatom_references()` prints `[skip] network.kegg.rds (exists)` and
`gatom_refs()` dies with `readRDS: error reading from connection`, which names neither the
file nor the fix.

**Suggested fix** — download to `paste0(dest_file, ".part")` and rename on success; treat a
zero-byte (or `readRDS`-unparseable) existing file as absent.

### 5. `gatom_refs(download = TRUE)` can silently ignore what it just downloaded  [severity: medium]
**R/gatom.R:106-127** — with `dir = NULL`, the download goes to
`dir %||% default_dir` = `"00_data/references/gatom"`, but the resolution order is
`unique(c(dir, "/opt/gatom-refs", default_dir))` — i.e. the staged container copy is
searched *before* the directory just written. `download = TRUE` then has no visible effect.

**Evidence** — old-vs-new is not applicable (new code); quoting the two statements:
`download_gatom_references(dest_dir = dir %||% default_dir, ...)` followed by
`search_dirs <- unique(c(dir, "/opt/gatom-refs", default_dir))`. Also note `test-gatom.R:270`
and :303 always pass `dir = "/opt/gatom-refs"` explicitly, so no test exercises the
`dir = NULL` + `download = TRUE` path.

**Failure scenario** — inside the container, `gatom_refs("Homo sapiens", download = TRUE)`
fetches fresh references into `00_data/references/gatom` and then loads the older
`/opt/gatom-refs` copies; module contents change with the image, not with the download, and
`print(refs)` is the only clue.

**Suggested fix** — when `download = TRUE`, put the download destination first in
`search_dirs` (and say in the docs that `/opt/gatom-refs` otherwise wins over the default
dir).

### 6. Three separate implementations of "require a Suggests package"  [severity: low]
**R/de-utils.R:71 (`.de_require`), R/gatom.R:45 (`.gatom_require`), plus raw
`requireNamespace()` blocks at R/gs-score.R:80, R/gs-test.R:308, R/gsdb-load.R:51 & :79,
R/gsdb-msigdb.R:33, R/gsdb-rebuild.R:233** — `.gatom_require()` is `.de_require()` with the
installer hardcoded to a `pkg == "gatom"` special case, which means the mwcsr/igraph
branches emit `install.packages(...)` from a switch rather than from the caller. This is the
utility-duplication class the brief calls out, one step short of a name collision.

**Evidence** — `grep -rn '_require <- function' R/` → the two definitions; the two bodies
differ only in how `install` is chosen (`.de_require` takes it as an argument, and all eight
call sites pass the right string, including `BiocManager::install("edgeR")`).

**Failure scenario** — none; maintainability only. The next person to add a Bioconductor
Suggests has three prior patterns to copy and no single place to fix wording.

**Suggested fix** — keep `.de_require()`, delete `.gatom_require()` (rename the survivor to
something layer-neutral, e.g. `.need_pkg()`), and route the five raw blocks through it.

### 7. `gs_write()` documents an invisible return but returns visibly  [severity: low]
**R/gs-io.R:32 (`@return The directory, invisibly, ...`) vs :94
(`structure(dir, files = files)`)** — `gs_save()` and `.gs_write_tsv()` do use
`invisible()`; `gs_write()` does not, so the path plus its whole `files` attribute prints at
the console on every call.

**Evidence** — from the run in finding 3, each `gs_write()` call auto-printed:
```
[1] "/tmp/gsw_stale"
attr(,"files")
[1] "/tmp/gsw_stale/by_contrast/A_vs_B/gsea_demo.tsv" ...
```

**Failure scenario** — none; noise in scripts and a docs/behaviour mismatch. Would also be
caught by an `expect_invisible()` in test-gs-io.R, which is absent (test-gatom.R:387 uses
one, so the convention exists).

**Suggested fix** — wrap in `invisible()` and add `expect_invisible(gs_write(...))`.

### 8. `pheatmap` is in Suggests but never used  [severity: low]
**DESCRIPTION:45** — `grep -rn "pheatmap" R/` hits exactly one line, a roxygen sentence in
R/gs-plot-heatmap.R:5 explaining that the *old* code was built on pheatmap. No code path
loads it.

**Evidence** — per-package usage count over `R/`: `pheatmap: 0` (every other Suggest is
non-zero; all 14 Imports are used, `methods` via R/deprecated-plot.R:34).

**Failure scenario** — none; maintainability. It invites `R CMD check` noise and misleads
readers into thinking a pheatmap path still exists.

**Suggested fix** — drop `pheatmap` from Suggests; reword the roxygen line to past tense.

---

## Not reviewed / caveats
- I did not test `gatom_refs()`/`gatom_module()` against a *second* reference version, so
  finding 5 is argued from source order rather than from an observed wrong module.
- Finding 4's truncation half is argued from source; I did not simulate a mid-transfer
  network failure (no outbound network in the container).
- The double `set.seed()` in `gatom_module()` (R/gatom.R:381-385) is correct and its comment
  explains *why* — `test-gatom.R` proves seed stability at both toy and realistic scale, and
  the BUM-fit test asserts the absence of the "parameters on the limit" warning. No finding.
