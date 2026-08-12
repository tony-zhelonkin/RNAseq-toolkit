# CoReSh: extraction plan, and the package-versus-skill decision

**Date:** 2026-08-11 · **Status:** plan of record, nothing implemented
**Parent:** [../2026-08-10-bulkirna-package/00_INDEX.md](../2026-08-10-bulkirna-package/00_INDEX.md) — this
document is Phase 4's CoReSh branch, and it also opens the Phase 5 (skills) question.

Evidence base: three read-only surveys on 2026-08-11 covering the DC-nexus `coresh-slice`,
the `14839-DM-cGAS` CoReSh scripts, the `coresh-signature-search` skill, and the
SciAgent-toolkit refactor. Every claim below carries a `file:line` so the next reader can
check it rather than trust it.

---

## 0. The answer to "does CoReSh work properly"

**The engine is correct and identical everywhere. The wrappers around it have drifted, and
one of them silently drops four columns.**

Three findings, in order of consequence:

1. **The 269-line CoReSh engine is already a de-facto package, vendored four times.**
   `lib/{symbols_to_entrez,coresh_batch,extract_gene_loadings}.R` are **md5-identical**
   across DC_hum_verse and `14839-DM-cGAS` (`coresh_batch.R` =
   `e34d1dc08deca13c3e98b9ca2c67e5b6` in four locations). `14839`'s `07_coresh_search.R:21-23`
   sources them out of `01_modules/.ref/dc-hum-verse-tooling/…/coresh-slice/lib` — one project
   vendors another project's scripts directory as a git ref. Zero project-specific lines.
   This is the highest-value extraction target in the repo and it needs packaging, not
   refactoring.

2. **The `neg_log_padj` defect is a column-projection bug, not an arithmetic one.**
   Verified directly in the live table, not inferred. See §1.

3. **Nothing CoReSh exists in `bulkiRNA` today.** No `R/coresh-*.R`, no export, no mention
   outside a single incidental row count in the parent plan. The folding has not happened,
   so there is nothing to un-fold — the question is entirely forward-looking.

---

## 1. The `neg_log_padj` defect, confirmed

`08_coresh_derived_gsea.R` contains no `neg_log_padj` expression at all. The column arrives
correct and then gets thrown away.

`normalize_gsea.R:214-215` computes it properly, so `nd` at `08_coresh_derived_gsea.R:157`
carries a valid value (`padj = 2.6e-33` → `32.58`). Then `08_coresh_derived_gsea.R:215-219`
projects to a ten-column allowlist that omits it:

```r
gsea_cols <- c("pathway_id", "pathway_name", "database", "nes", "pvalue", "padj",
               "set_size", "core_enrichment", "contrast", "direction")
append_coresh(file.path(mst_dir, "master_gsea_table.csv"),
              coresh_rows[, intersect(gsea_cols, names(coresh_rows)), drop = FALSE], ...)
```

`append_coresh:75` then calls `dplyr::bind_rows(existing, new_rows)`, which unions columns
and NA-fills the absent ones. The master table has fourteen columns, so every CoReSh row
reads `padj=2.6175484e-33, neg_log_padj=NA`. The same allowlist also nulls
`leading_edge_size` and `gene_ratio`; `master_unified.csv`'s separate 13-column `uni_cols`
(`:222-226`) recovers `genes_full_set` and `leading_edge_size` and still omits
`neg_log_padj` and `gene_ratio`.

**"Non-finite" means `NA` from a column union.** The fix is one line — or, better, deleting
the projection, which is what `1.12_coresh_derived_gsea.R:308-335` already does.

### The second-order drift, which matters more

One column, two saturation conventions:

| Where | Expression | Ceiling |
|---|---|---|
| `normalize_gsea.R:214-215` (CoReSh, custom DBs) | `-log10(padj)` then `[is.infinite] <- 16` | **16** |
| `05_gsea_msigdb_run.R:214` (migrated MSigDB) | `-log10(pmax(padj, .Machine$double.xmin))` | **≈307.65** |
| `12_gsea_viz.R:480` (the consumer) | `-log10(pmax(padj, .Machine$double.xmin))` | ≈307.65 |
| `08_coresh_derived_gsea.R` | absent | `NA` |

`12_gsea_viz.R:487` sizes points by this column, so once `08` is fixed the CoReSh rows still
render on a different scale from the MSigDB rows in the same figure. **Pick one convention
before extracting.** The recommendation is `pmax(padj, .Machine$double.xmin)`: it is
information-preserving, it already matches the migrated script and the renderer, and a cap
at 16 is a silent lie about anything below `1e-16`.

---

## 2. Duplication and drift across the consumers

`08_coresh_derived_gsea.R`'s own comments admit the copying — `:82` *"mirrors 1.12:105-128"*,
`:134` *"mirrors 1.12:228-280"*, `:211` *"mirrors 1.12:308-335"*.

| Unit | `1.12` (DC-nexus) | `08` (cGAS) | Drift |
|---|---|---|---|
| `parse_gmt()` | `:105-128` | `:82-93` | `08` strips blank lines and returns `NULL` on empty; `1.12` stops. `1.12` auto-prefixes `CORESH_`; `08` does not. |
| ranked vector | local `:167-209` | `de_gsea_helpers.R:31-46` | `1.12` also drops `is.infinite`; the helper drops only `NA`. `1.12` accepts four p-value column spellings, the helper one. |
| **GSEA call** | `:238-250` | `:147-150` | **`nPermSimple` 1000 vs 100000 (100×); `minGSSize` 10 vs 15.** `1.12` sets no seed; `08` does, and caches. |
| `rename(nes = NES)` | `:294-296` | `:162` **and** `append_coresh:73-74` | `08` does it twice. |
| master append | `:308-335`, all columns, `master_unified` only | `append_coresh:65-77`, allowlist, both masters | The defect in §1. `08` also bypasses the project's own `append_master_table()` (`de_gsea_helpers.R:117`) and its schema validation and numeric rounding. |
| summary | `n_total/n_sig_05/n_sig_01/n_up/n_down` | `n_sets_tested/n_sig_fdr/n_up_sig/n_down_sig/top_set/top_nes` | Different names for identical quantities. |
| rownames guard | absent | `:97-99`, `^ENSMUSG` | `08` only, hardcoded to mouse. |

The query-side pair (`1.11` vs `07_coresh_search.R`) agrees on every constant —
`TOP_N = 5`, `n_top = 50`, `min_size = 15`, `max_size = 500`, `jaccard = 0.8`,
`MIN_QUERY_SIZE = 3` — differing only in hardcoded-versus-`CONFIG`, human/`hsa` versus
mouse/`mmu`, and skip-versus-stop on missing chunks. **Agreed constants are safe package
defaults.** `nPermSimple` and `minGSSize` disagree and need a deliberate call.

Line accounting across the 1,404 surveyed lines: **~265 (19%) genuinely project-specific,
~1,139 (81%) generic machinery**, and the 269-line `lib/` is 0% project-specific.

### Latent defects worth fixing during extraction

- `1.11:453` — `max(which(grepl("^GSE", tokens)))` returns `-Inf` on no match, and the
  downstream `seq_len()` errors. Unreachable only because `build_coresh_gmt` always emits a
  GSE token. Both query scripts reconstruct provenance by re-parsing
  `CORESH_<query>_<GSE>` set names when they still hold `top_hits` — keep the table instead.
- `extract_gene_loadings.R:29-31` stops unless `sym2ent`/`ent2sym` are already in the global
  environment. Package namespacing removes the guard and the problem together.
- `extract_gene_loadings.R:35` calls `digest::digest` without declaring `digest`, and
  memoizes into `options()`. A package-local environment is the right cache.
- `sym2ent`'s `is.null(db)` branch is dead: the `::` errors first, so the friendly
  BiocManager message never prints.
- `coresh_batch.R:44` calls **`fgsea:::gesecaCpp`**, an unexported internal. This is the one
  real packaging hazard — `R CMD check` flags `:::` on another package, and an fgsea release
  can remove it without notice. See §5.
- The bridge documents a "drop sets >90% overlap with the query" rule that no code
  implements (`coresh-to-gsea-bridge.md:113`).
- `validate_coresh_install.R:25-29` requires the `coresh` package, which no script ever
  calls — the engine reimplements the kernel. Decide whether that dependency is real.

---

## 3. The decision: package owns the computation, skill owns the judgement

**Recommendation: extract the engine into `bulkiRNA` and keep `coresh-signature-search` as a
skill, rewritten to call the package.** Both, with a clean line between them.

The line: **anything that takes a chunk, a query, or a table and returns a table is a package
function. Anything that takes a result and returns a decision stays in the skill.**

The SciAgent-toolkit refactor states the same rule in its own words
(`docs/packaged-skills.md:22-31`): package it "when it ships executable logic that must run
identically every time"; when the content is "guidance, conventions, a decision tree, a thin
wrapper over a stable public API — a flat SKILL.md is correct." CoReSh is both things at
once, which is exactly why it splits cleanly.

Three facts settle it:

- **`coreshMatch()` is triplicated verbatim** — `SKILL.md:115-127`,
  `coresh-to-gsea-bridge.md:26-47`, `coresh_batch.R:32-51`. A numerical kernel pasted into
  prose in two places is a kernel with three futures.
- **The engine is already byte-identical in four checkouts**, re-vendored by hand with a
  `cp` instruction in a comment (`1.11:32-34`). That comment is the argument for packaging,
  written by the person who felt the cost.
- **The judgement content is genuinely irreducible.** Choosing whether a confounder class
  *is* the biology, pre-registering expected top-hit categories before running, the
  Expected/Adjacent/Novel/Cell-line/Noise taxonomy with its ≥5/10 rubric, the red-flag
  table, the Synapse token walkthrough with a human in the loop. None of it has a return
  value. An R function cannot hold it and should not try.

### Tradeoffs, honestly

**In favour of the package:** one tested copy; version pinning that makes a result citable;
`R CMD check` and `testthat` on the kernel; the `:::` hazard confined to one guarded place;
the source-order `stop()` guards disappear; `SKILL.md` shrinks to calls and its prose stops
drifting from its code.

**Against, and the answers:**

- *Dependency weight.* CoReSh needs `qs2`, `BiocParallel`, `data.table` and `digest` on top
  of what `bulkiRNA` already imports. Answer: **`Suggests`, not `Imports`** — `fgsea` is
  already an Import, and `AnnotationDbi`, `org.Hs.eg.db`, `org.Mm.eg.db`, `homologene` and
  `babelgene` are already in `Suggests`. The marginal cost is four suggested packages, and
  the `gatom_*` family already sets the precedent for a feature that installs on demand.
- *Scope creep — is CoReSh "bulk RNA-seq gene-set analysis"?* Its output is a GMT that feeds
  `gsdb_from_file()` and `gs_test()`. It is a **gene-set provider**, which is precisely the
  layer `bulkiRNA` already has, and `gsdb_coresh()` would sit beside `gsdb_msigdb()` without
  bending the architecture.
- *The 20 GB Synapse chunk tree.* A package cannot ship it and must not try. Answer: the
  same shape as `download_gatom_references()` — resolve from `CORESH_CHUNKS`, validate,
  and fail with instructions. Tests skip without the tree.
- *Two audiences.* An agent reading a skill and a human writing a script want different
  affordances. Answer: that is the split, not an objection to it.

**Against a package-only move (no skill):** the interpretation protocol is the part that
stops a CoReSh sweep from producing confident nonsense, and it has no home in `help()`.

**Against a skill-only move (status quo):** the kernel keeps being pasted, `neg_log_padj`-class
defects keep appearing per-project, and there is no version to cite in a methods section.

### Scope decision to make first

`gsdb_coresh()` inside `bulkiRNA`, or a sibling package `coreshtools`? **Recommendation:
inside `bulkiRNA`**, because the output is a `gs_db`, the split would duplicate the ID-mapping
and GMT machinery, and Phase 3 is about to install exactly one package into the image. Revisit
only if the chunk-I/O dependency set grows past `Suggests`.

---

## 4. Proposed surface

Names follow the house `verb_noun` / layer-prefix convention. The provider is the headline;
the rest are the primitives it composes.

**Provider layer — the entry point**

```r
gsdb_coresh(queries, chunk_dir = Sys.getenv("CORESH_CHUNKS"), species = "human",
            top_hits = 5L, n_top = 50L, min_size = 15L, max_size = 500L,
            jaccard_threshold = 0.8, max_query_overlap = 0.9, n_cores = 4L)  -> gs_db
```

Returns a `gs_db` with `database = "coresh"`, so `gs_test()`, `gs_plot_*()` and `gs_write()`
work unchanged, and the provenance table (`query_name`, `gse`, `chunk_path`,
`loading_cutoff`, `rank_in_coresh`) rides along as an attribute rather than being
reconstructed from set names.

**Search layer**

```r
coresh_search(queries, chunk_dir, n_cores = 4L, pvalues = FALSE)   -> tibble
coresh_match(obj, query, pvalues = FALSE, sample_size = 21L, seed = 1L, eps = 1e-300)
coresh_chunks(chunk_dir, cache = TRUE)                              -> tibble(gse, chunk)
coresh_validate(chunk_dir = Sys.getenv("CORESH_CHUNKS"), species = c("hsa", "mmu"))
coresh_convergence(ranking, top_n = 10L, min_queries = 2L)          -> tibble
```

**Set-building layer**

```r
coresh_loadings(chunk_path, gse_id, query, n_top = 50L)             -> tibble
coresh_sets(top_hits, queries, chunk_dir, ...)                      -> gs_db
```

**Gene-identifier layer (generally useful, beyond CoReSh)**

```r
gene_to_entrez(symbols, species = "human", multi_vals = "first")    -> integer
entrez_to_gene(entrez, species = "human")                           -> character
map_orthologs(symbols, from = "human", to = "mouse",
              method = c("homologene", "babelgene", "biomart"))
filter_confounder_genes(symbols,
    drop = c("ribosomal", "mito", "hemoglobin", "cell_cycle", "sex"))
```

The confounder patterns become a documented data object rather than a regex retyped in
`query-design.md:53-64` and `1.11:255`.

**Shared machinery the consumers keep re-implementing** (needed by `06`, `08` and `12` as
much as by CoReSh — today there are two `parse_gmt`s, two `build_ranked_vector`s and three
master-append implementations):

- `gsdb_from_file()` already covers `parse_gmt` and `gs_ranks()` already covers
  `build_ranked_vector` — the migration is to **delete** the local copies, not to port them.
- One derived-column contract in one place: `neg_log_padj`, `gene_ratio`,
  `leading_edge_size`, `genes_full_set`, one saturation convention, and **no column
  allowlist before `bind_rows`**. That is the whole §1 bug, and `gs_to_master()` is where it
  belongs once lifted into `de_gsea_helpers.R`.

`Suggests` grows by `qs2`, `BiocParallel`, `data.table`, `digest`.

---

## 5. The `fgsea:::gesecaCpp` hazard

`coresh_batch.R:44` reaches into fgsea's internals for the GESECA p-value. Options, best
first:

1. **Reproduce the p-value through public `fgsea::geseca()`** and prove equivalence on a
   fixture chunk, the same way `B2-fgsea-equivalence.R` discharged the clusterProfiler
   question. If the numbers match, the hazard is gone permanently.
2. **Keep `:::` behind a guard** — a version assertion plus `exists()`, erroring with a
   named workaround, in the shape of `.msigdbr_assert_ortholog_coverage()`. Correct but
   permanently fragile.
3. Vendor the C++. Rejected: it forks a numerical kernel.

`pvalues = FALSE` (the `pctVar` screen) touches none of this, so the provider stays useful
even if the p-value path degrades. Attempt 1 and fall back to 2.

---

## 6. What the skill becomes

`coresh-signature-search` stays, rewritten against the package and against the refactored
toolkit's conventions (canonical frontmatter only — `name`, `description`, `license`,
`allowed-tools`, `compatibility`; description ≤350 characters with no angle brackets; prose
`## See also` bullets; SKILL.md under ~400 lines).

It keeps: when-to-use and contraindications, the routing to `bulk-rnaseq-gsea` and RummaGEO
and gatom, the web-UI-first smoke test, query design and confounder judgement, pre-registering
expected categories, the top-20 taxonomy and its rubric, negative-control design, the
red-flag table, the Synapse walkthrough, and bridge failure diagnosis.

It loses: `scripts/` (four R files become package functions), the triplicated `coreshMatch`
body, and the manual `cp` re-vendoring instruction.

**Note for whoever edits it:** the copy under
`scbio-docker/toolkits/SciAgent-toolkit` is **ahead** — its `SKILL.md` carries the
post-refactor frontmatter, and the Meta-Aging checkout at `01_modules/` is an ancestor.
`references/` and `scripts/` are byte-identical in both. Edit the scbio-docker copy.

### And the `bulkiRNA` skill the user wants

The natural home is **a rewrite of the existing `bulk-rnaseq-gsea` skill**, not a new one.
It already has the right four-reference shape (219-line SKILL.md plus msigdb, custom-db,
master-tables, visualization), and it currently instructs
`source(file.path(DIR_TOOLKIT, "GSEA/GSEA_processing/run_gsea.R"))` (`SKILL.md:76-77`) and
points at `docs/GSEA-workflow/` — both of which the package retired. **`bulkiRNA` appears
nowhere in the toolkit today.** `annotate-bulk-rnaseq-data` carries the same coupling and
pins "RNAseq-toolkit v0.2.0"; the toolkit already tracks fifteen source-coupled skills as
known debt (commit `9abee81`). A link-integrity test now gates cross-skill references, so
`bulk-rnaseq-pathway-explorer:235`, `te-geneset-gsea`, `star-te-preprocessing` and
`gatom-metabolomic-predictions:99` need updating in the same pass.

---

## 7. Sequencing, folded into the parent plan

CoReSh extraction sits **after** Phase 3 and runs alongside Phase 4, because a package
function is only worth writing once the image can deliver it.

| Step | Work | Gate |
|---|---|---|
| **C0** | Settle the `neg_log_padj` convention and one derived-column contract; lift `gs_to_master()` into `de_gsea_helpers.R`. | Decision recorded here |
| **C1** | Fix `08_coresh_derived_gsea.R`'s allowlist as part of its Phase 4 migration — drop the projection. Re-run and diff. | 492 CoReSh rows carry finite `neg_log_padj` |
| **C2** | Discharge the `gesecaCpp` question (§5) on a fixture chunk. | Equivalence numbers recorded, or the guard written |
| **C3** | Port the 269-line engine into `R/coresh-*.R` with the §4 surface, plus tests that skip without the chunk tree. | `devtools::test()` green, `verify_golden.R` exit 0, `R CMD check` clean of `:::` |
| **C4** | Add `gsdb_coresh()` and prove it against DC-nexus's existing GMT — the same set names and memberships from the same queries. | Byte-level agreement on set contents |
| **C5** | Rewrite the skill against the package; delete `scripts/`. | `sciagent validate`, `tests/run-all.sh` |
| **C6** | Rewrite `bulk-rnaseq-gsea` (and `annotate-bulk-rnaseq-data`) onto `library(bulkiRNA)`; fix the four cross-links. | Link-integrity test passes |

C0 and C1 are cheap and unblock Phase 4. C3 needs the chunk tree mounted, so it needs a
container with `/data2/users/shared/refcache/coresh/current/preprocessed_chunks` available.

---

## 8. Open decisions for the integrator

1. **`neg_log_padj` convention** — `pmax(padj, .Machine$double.xmin)` (recommended, matches
   the migrated script and the renderer) or the cap at 16.
2. **`nPermSimple` and `minGSSize`** — `1.12` uses 1000/10, `08` uses 100000/15. One default.
   100000 is the safer statistic and the slower run.
3. **Scope** — `gsdb_coresh()` in `bulkiRNA` (recommended) or a sibling package.
4. **Whether the `coresh` GitHub package is a real dependency**, given no script calls it.
5. **Whether C5/C6 wait** on the wider SciAgent-toolkit refactor. That repo is a separate
   submodule, currently 26 commits ahead of `origin/dev` and dirty, so the skill rewrites
   need their own go-ahead.
