# CoReSh: extraction plan, and the package-versus-skill decision

**Date:** 2026-08-11 · **Status:** plan of record; C0-C3 done, C4 next, C5/C6 blocked. See §13 for C3's second half. §10/§11's p-value conclusion is **corrected by §12** — read §12 first.
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
`TOP_N = 5`, `top_n = 50`, `min_size = 15`, `max_size = 500`, `jaccard = 0.8`,
`MIN_QUERY_SIZE = 3` — differing only in hardcoded-versus-`CONFIG`, human/`hsa` versus
mouse/`mmu`, and skip-versus-stop on missing chunks. **Agreed constants are safe package
defaults.** `nPermSimple` and `minGSSize` disagreed; the call is now made — 100000 and 10
(§8.2).

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
- `validate_coresh_install.R:25-29` requires the `coresh` package. **Corrected twice.** This
  bullet first doubted the validator; §8.4 then defended it; §11 settles it. The validator
  requires a package that exports nothing at all — `alserglab/coresh` v0.1.0 has no `R/`
  directory — so neither `coreshMatch()` nor `queryGSE()` can be called from it, and the
  skill's documented R path cannot run as written. That is a C5 defect.

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
            top_hits = 5L, top_n = 50L, min_size = 15L, max_size = 500L,
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
coresh_loadings(chunk_path, gse_id, query, top_n = 50L)             -> tibble
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

> ✅ **Settled 2026-08-13 by measurement — see §10.** Neither option below is what happened:
> `pctVar` needs no fgsea entry point at all, and the p-value path delegates to
> `coresh::coreshMatch()`. The options are kept because the reasoning that ranked them is what
> the measurement tested.

`coresh_batch.R:44` reaches into fgsea's internals for the GESECA p-value. Options, best
first:

1. **Reproduce the p-value through public `fgsea::geseca()`** and prove equivalence on a
   fixture chunk, the same way `B2-fgsea-equivalence.R` discharged the clusterProfiler
   question. If the numbers match, the hazard is gone permanently. The skill's own SKILL.md
   describes CoReSh as "a thin R wrapper around `fgsea::geseca`", so the public entry point
   is very likely sufficient and the `:::` call is a shortcut rather than a necessity.
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
| **C0** | ~~lift `gs_to_master()` into `de_gsea_helpers.R`~~ → **superseded by [ADR-002](../../adr/ADR-002-master-table-schema.md)**: the versioned schema, `gs_to_master()` and `gs_validate_master()` go **into the package**. `neg_log_padj = -log10(pmax(padj, .Machine$double.xmin))`. | `devtools::test()` green; the validator rejects the ten-column projection |
| **C1** | ✅ **done 2026-08-12** — migrated `08_coresh_derived_gsea.R` onto `bulkiRNA` and dropped the projection. See §9. | ✅ 708 CoReSh rows, all carrying finite `neg_log_padj` |
| **C2** | ✅ **done 2026-08-13** — discharged on real chunk data. `pctVar` becomes package-owned arithmetic over the stored `totalVar`; p-values delegate to `coresh::coreshMatch()`, so no `:::` enters our namespace and §5's guard is dropped. See §10. | ✅ equivalence recorded: `pctVar` agrees to ≤1.03e-5 relative, cause fully accounted for |
| **C3** | 🟡 **half done 2026-08-13** — `R/coresh.R` (search layer) and `R/gene-ids.R` landed; the set-building layer, `coresh_loadings()` and `coresh_sets()`, is still to write. See §11. | ✅ 1061 tests pass, golden 20/20, no `:::` in the package; `pct_var` bit-identical to a hand-computed value on real chunks |
| **C4** | Add `gsdb_coresh()` and prove it against DC-nexus's existing GMT — the same set names and memberships from the same queries. | Byte-level agreement on set contents |
| **C5** | 🚫 **Blocked** on the SciAgent-toolkit refactor (§8.5). Rewrite the skill against the package; delete `scripts/`. | `sciagent validate`, `tests/run-all.sh` |
| **C6** | 🚫 **Blocked**, same reason. Rewrite `bulk-rnaseq-gsea` (and `annotate-bulk-rnaseq-data`) onto `library(bulkiRNA)`; fix the four cross-links. | Link-integrity test passes |

C0 and C1 are cheap and unblock Phase 4. C3 needs the chunk tree mounted, so it needs a
container with `/data2/users/shared/refcache/coresh/current/preprocessed_chunks` available.

---

## 8. Decisions — all settled as of 2026-08-12

1. **`neg_log_padj` convention** — ✅ **settled: `-log10(pmax(padj, .Machine$double.xmin))`.**
   Information-preserving, and already what two of three call sites do. The cap at 16 is
   retired. Now owned by [ADR-002](../../adr/ADR-002-master-table-schema.md).
2. **`nPermSimple` and `minGSSize`** — ✅ **settled: `nPermSimple = 100000`, `minGSSize = 10`.**
   The safer statistic on the permutation side, the more inclusive threshold on the set-size
   side, so neither script's value carries over wholesale: `1.12` keeps its `minGSSize`, `08`
   keeps its `nPermSimple`. Cost is runtime, which is the right thing to spend. Small sets stay
   testable, and `padj` is what protects against reading too much into them.
3. **Scope** — ✅ **settled: `gsdb_coresh()` lives in `bulkiRNA`.** Its output is a `gs_db`,
   so it enters at the provider layer and bends nothing.
4. **Whether the `coresh` GitHub package is a real dependency** — ⚠️ **settled: yes, then
   overturned by §11.** The package exists and is documented in the skill, but it exports
   nothing: v0.1.0 has no `R/` directory. It cannot be a dependency of anything, and the
   three consequences below are void. `bulkiRNA` neither imports nor suggests it, and no
   `Remotes:` line was added. The paragraph stands as written because the reasoning is what
   §11 tested.

   [`alserglab/coresh`](https://github.com/alserglab/coresh), used enough already that the
   question was mine, not the owner's. **My "no script calls it" was wrong**: it came from
   grepping analysis scripts only. `coresh-signature-search/SKILL.md:101` does
   `library(coresh)` and documents `coreshMatch()` as the R-package path for anything past one
   query, and `DC_hum_verse/docs/coresh-plan/.../consumer_quickstart.md:11` names
   `coresh::queryGSE()` as the call that silently returns nothing when misused. The skill
   inlines `coreshMatch()` only because the upstream vignette defines it inline.

   Three consequences:

   - **C3 shrinks.** Where `coresh` exports the primitive, `bulkiRNA` calls it rather than
     re-typing a 20-line body that has already been re-typed three times. Porting the engine
     means porting the orchestration — chunk iteration, Entrez conversion, result assembly —
     not the scoring kernel.
   - **It also weakens C2.** The `fgsea:::gesecaCpp` call is upstream's own, inside
     `coreshMatch()`. Calling `coresh::coreshMatch()` moves that `:::` out of our namespace and
     out of `R CMD check`'s reach. The equivalence question stays worth answering, but it stops
     blocking a clean check.
   - **`Remotes:` is needed after all.** §16 concluded otherwise because at that point nothing
     in `Suggests` was GitHub-only. `coresh` is, so `Remotes: alserglab/coresh` lands in
     `DESCRIPTION` with it, and the preflight row carries the `remotes::install_github()`
     command rather than a CRAN or Bioconductor one. It is not installed in `v0.5.13`; whether
     the image should carry it is a Phase 3 follow-up, not a blocker, since the guard degrades
     to a clear message.
5. **Whether C5/C6 wait** on the wider SciAgent-toolkit refactor — ✅ **settled: they wait.**
   That refactor is still in flight, so the skill rewrites are blocked and stay blocked until
   the owner says otherwise. Phase 3 has now discharged the *technical* prerequisite — nested
   submodules are forbidden, and a skill can now reach an installed `bulkiRNA` in
   `scdock-r-dev:v0.5.13` — so what remains is sequencing against the other track, not
   capability.

### Added by the 2026-08-12 ADR pass

6. **Chunk-tree resolution** — ✅ `Sys.getenv("CORESH_CHUNKS")` in §4 is **superseded** by the
   single `.ref_path("coresh", ...)` resolver: explicit argument →
   `$REFCACHE_ROOT/coresh/current` → an error naming the `refcache.sh` command. See
   [ADR-004](../../adr/ADR-004-reference-data-tiers.md). It must land before C4, and it means
   CoReSh needs **no new environment variable** and no host path anywhere in the package.
7. **Provenance** — ✅ whatever reads the chunk tree must record the **resolved snapshot tag**
   (e.g. `syn66227307_20260721`), not merely the fact that it followed `current`. That is how
   a re-run whose numbers move becomes explicable rather than mysterious.

---

## 9. C1 as executed (2026-08-12)

`08_coresh_derived_gsea.R` migrated onto `bulkiRNA` 0.4.0 in `scdock-r-dev:v0.5.13`:
`parse_gmt()` → `gsdb_from_file()`, `clusterProfiler::GSEA()` → `gs_ranks()` + `gs_test()`,
`normalize_gsea_results()` → `gs_to_master()`, and both allowlists deleted.

### The gate, and more than the gate

| | Before | After |
|---|---|---|
| CoReSh rows | 492 | **708** |
| finite `neg_log_padj` | **0** | **708** |
| finite `gene_ratio`, `leading_edge_size` | 0 | 708 |
| max `neg_log_padj` | — | 37.55 |
| columns (gsea / unified) | 14 / 13 | 14 / 15 |

The 37.55 is worth a second look: it is above the retired cap of 16, so it also demonstrates
the ADR-002 convention doing its job on real values rather than in a test.

### The row count moved, and the reason was a second defect

492 → 708 is not the allowlist. The old rows covered **41 of the GMT's 59 sets**, and nine of
those 41 named sets that **are not in the GMT at all**.

The cause is a stale `load_or_compute()` cache. `gsea_coresh_<contrast>.rds` is dated
**Jun 18 15:39**; `coresh_derived_sets.gmt` is dated **Jun 22 00:14**. The cached GSEA results
predate the gene sets they claim to describe by four days, and `load_or_compute()` invalidates
on nothing but the filename. Every CoReSh master row written after 22 June, the 11 August
rewrite included, was re-persisted from those objects. So the table simultaneously described
nine sets that no longer existed and omitted 27 that did.

**This generalises past CoReSh.** Any `load_or_compute()` cache whose input file changes
without its key changing has the same exposure. Worth a survey of the other cached stages
before trusting a number that came out of one.

### Two fixes the migration forced, neither of them planned

- **`gs_ranks()` does not drop empty gene names, and `build_ranked_vector()` did.** fgsea
  rejects `""` in `names(stats)`, so the first migrated run threw on all twelve contrasts. The
  script now cleans the DE table first, as `05_gsea_msigdb_run.R` already did.
- **A failure must not be cached, and a total failure must not look like a null result.** That
  first run went down the no-enrichment soft path: it purged all 492 CoReSh rows from both
  master tables and **exited 0**. `load_or_compute()` had also cached the `NULL` from each
  errored contrast, so the failure would have survived the fix that repaired it — the next run
  would read the files, skip every contrast, and report the same clean nothing. The script now
  deletes a failed contrast's cache entry and stops outright when every contrast errored,
  leaving the master tables untouched.

  This is the §0 lesson again, and the third time in this refactor: the empty result and the
  broken run were indistinguishable from inside the script, so the broken run took the empty
  result's path. Counting the errors is what separates them.

### What was left alone

The twelve pre-migration `gsea_coresh_<contrast>.rds` files are now unused but not deleted —
they are the owner's data and the evidence for the staleness finding above. The project's
`analysis_config.yaml` keeps `gsea_min_size: 15`; only the script's fallback moved to 10, to
match the package default. `14839-DM-cGAS` remains uncommitted, with fresh
`master_*.pre-coresh-c1.20260812-225911.csv` backups next to the tables.

---

## 10. C2 as executed (2026-08-13) — the `gesecaCpp` question, discharged

Run in `scdock-r-dev:v0.5.13` with the chunk tree mounted read-only. Resolved snapshot
`syn66227307_20260721`, 85 mouse chunks, 500 datasets in `chunk_1`. Measured on both synthetic
matrices with planted signal and real chunk data.

**The answer splits by quantity, and the split is better than either §5 option.**

### `pctVar` — the default screen, and the one that matters

Formula-identical to `fgsea::geseca()`'s own `pctVar`. On synthetic matrices where both sides
divide by the same total variance the two agree bit-for-bit. On real chunks they differ by at
most **1.03e-5 relative**, and the cause is fully accounted for: `coreshMatch()` divides by the
stored `obj$totalVar`, computed before the 1024-step quantization, while `geseca()` recomputes
`sum(E^2)` from the quantized matrix. Stored versus recomputed, on four datasets:

| GSE | stored `totalVar` | `sum(E^2)` | relative difference |
|---|---|---|---|
| GSE10000 | 22115.51873486 | 22115.74709606 | 1.033e-05 |
| GSE10000 | 80427.45593222 | 80427.41750622 | 4.778e-07 |
| GSE10001 | 11641.66123124 | 11641.77372837 | 9.663e-06 |
| GSE100012 | 37353.52744422 | 37353.56605911 | 1.034e-06 |

The stored value is the more faithful of the two, so the conclusion is not "call `geseca()`
instead" but **"compute `pctVar` in-package from the stored `totalVar`"** — three arithmetic
lines over `colSums()`, no fgsea entry point at all, public or internal. `pvalues = FALSE`
therefore needs neither `fgsea:::` nor the `coresh` package.

### `pval` — not reproducible as a value, and it does not need to be

Option 1 fails here, for three separate reasons:

- Both paths are Monte-Carlo estimates of the same tail and neither is deterministic. Across
  five seeds at `sampleSize = 21` on a null matrix, `gesecaCpp` returned 0.4431 0.4658 0.3976
  0.3976 0.3749 and `geseca()` returned 0.3337 0.3337 0.3576 0.3437 0.3397. With planted
  signal both collapse to the same order of magnitude — 1e-28 to 1e-30 at boost 0.4, 1e-44 to
  1e-53 at boost 1.0 — and neither series contains the other's values.
- **`sampleSize` does not propagate through the public path.** `geseca()` returned an identical
  five-seed series at `sampleSize = 21` and `101`, while `gesecaCpp`'s spread narrowed as
  expected. The public function is not a wrapper over the internal one with the same estimator
  controls.
- **`geseca()` rejects every real chunk.** `checkGesecaArgs()` errors with "Duplicate
  rownames(E) are not allowed"; chunk Entrez ids repeat. The kernel indexes by `match()` and
  does not care. A `make.unique()` workaround is required before the public function will run
  at all, which is a second reason not to route through it.

### Resolution

**`pctVar` in-package from the stored `totalVar`; p-values delegated to
`coresh::coreshMatch()`.** That keeps the `:::` inside upstream's namespace, out of ours and out
of `R CMD check`'s reach — which is what decision 4 already implied and this measurement now
confirms. §5's option 2, a guarded `:::` in `bulkiRNA`, is **not needed and is dropped.** The
`coresh` guard degrades to a clear message on the p-value path only, and the screen keeps
working without it.

### Two incidental findings

- **`qs2` is absent from `scdock-r-dev:v0.5.13`.** Any chunk reader needs it, so it joins the
  Phase 3 follow-up alongside `coresh` itself. Installed to a scratch library for this run;
  `qs2` 0.2.2 pulls `stringfish`.
- **Chunk Entrez rownames are not unique.** Whatever `bulkiRNA` does with a chunk must index by
  `match()`, never by name lookup.

---

## 11. C3, first half (2026-08-13) — search layer and gene identifiers

Two new files each, written by delegated agents in isolated clones and gated here.
**C3 is not finished:** the set-building layer is still to come. What landed:

| File | Exports |
|---|---|
| `R/coresh.R` | `coresh_chunks()`, `coresh_match()`, `coresh_search()`, `coresh_convergence()`, `coresh_validate()` |
| `R/gene-ids.R` | `gene_to_entrez()`, `entrez_to_gene()`, `filter_confounder_genes()` |

Exports go from 64 to 72. Gates: **1061 tests pass, 0 fail, 5 skip; golden 20/20.**
`Suggests` gains `qs2` and `BiocParallel`, and the dependency registry gains a
`"coresh"` feature so `bulkirna_check_deps("coresh")` names them.

### The p-value premise collapsed, and the agent caught it

The brief told the agent to call `coresh::coreshMatch()`, following §8.4 and §10. It
refused, reporting that the function does not exist. It was right, and verification
made it starker: `alserglab/coresh` at tag `v0.1.0` contains `DESCRIPTION`,
`NAMESPACE` (`exportPattern("^[[:alpha:]]+")`), `README.md`, two licence files and one
file `vignettes/coresh-local.Rmd`. **There is no `R/` directory, so the package
exports nothing.** The vignette defines `coreshMatch()` inline and calls
`fgsea:::gesecaCpp` itself. Three commits, all April 2025.

So §8.4's "C3 shrinks — call the primitive rather than re-typing it" is void: there is
no primitive. `coresh::queryGSE()`, which
`DC_hum_verse/.../consumer_quickstart.md:11` names, does not exist either — a finding
for C5, since the skill's documented R path cannot run as written and
`validate_coresh_install.R` requires a package that provides nothing.

That leaves the p-value with no public route at all. `fgsea::gesecaSimple()` is a plain
`nperm` permutation test whose p-value floors at roughly `1/nperm`, so it cannot
express the magnitudes CoReSh reports; `fgsea::geseca()` was already ruled out in §10.
**The only route is the unexported `fgsea:::gesecaCpp`, and whether this package takes
on that call is now an open decision for the owner** — §5 option 2, back from the dead
because option 1 is measurably impossible rather than merely inelegant.

Meanwhile `pvalues` stays in both signatures and errors immediately, before a chunk is
read, with a message that states all three dead ends and points at `pvalues = FALSE`.
`p_value` stays as an always-`NA_real_` column so wiring the path later changes no
result shape. Both real consumers (`07_coresh_search.R:142`, `1.11:347`) already pass
`pvalues = FALSE`, so nothing in flight is blocked.

### Verified against the real chunk tree, not just the test suite

Snapshot `syn66227307_20260721`, mouse, two chunk files, 1,000 datasets, indexed
through `REFCACHE_ROOT` and `.ref_path()`.

- `coresh_search()` returned 2,000 rows over two queries with contiguous per-query
  ranks and provenance carrying the snapshot tag.
- **`pct_var` is bit-identical to an independently hand-computed value in 12 of 12
  datasets**, maximum absolute difference exactly 0, including the zero-overlap rows.
  My first comparison disagreed and was wrong, not the code: `GSE10000` appears twice
  under different platforms, and the check matched on `gse` alone.
- `coresh_validate()` reported all eight checks, `coresh` among them, marked not
  installed and annotated with the fact that installing it would not help.
- `gene_to_entrez()` mapped mouse symbols to integers and warned that `Tmem173` does
  not map — it is the retired alias of `Sting1`, so the warning is doing its job.

### Two defects fixed in the ported rules

- **The documented hemoglobin regex removes a growth factor.**
  `^HB[ABDEG][0-9]*`, as written in `query-design.md`, matches `HBEGF` and `HBS1L`.
  Anchored to `^HB[ABDEGMQZ][0-9]*$|^HB[ABDEGMQZ]-`, it excludes both, still matches
  the mouse `Hba-a1` spelling, and gains mu, theta and zeta globin, which a filter
  named after hemoglobin should not have been missing. A regression test pins it.
- **The pseudogene rule was left unimplemented**, as briefed. `^[A-Z0-9]+P[0-9]+$`
  matches `DUSP1` and `RANBP1`.

`filter_confounder_genes()` has **no `species` argument**. Matching is
case-insensitive, so one rule set covers `MT-ND1` and `mt-Nd1` alike; an argument that
only validated its own value would imply the rules changed with it.

### Still to do in C3, and C4

The set-building layer — `coresh_loadings()` and `coresh_sets()`, from
`extract_gene_loadings.R` — is not written. `gsdb_coresh()` (C4) sits on top of it, so
both remain open. `.ref_path()` is already in place, so decision 6 is satisfied.

---

## 12. Correction (2026-08-13): the p-value path is public after all

**§10 and §11 reached the wrong conclusion, on a measurement I did not make.** The claim was
that `fgsea::geseca()` "does not honour `sampleSize`", and from that followed "the only route
to the published p-value is `fgsea:::gesecaCpp`, and adopting that internal call is a pending
owner decision".

The premise is false. `geseca()` honours `sampleSize`; its multilevel estimator's precision
tracks it exactly as documented:

| `sampleSize` | `log2err` across three seeds |
|---|---|
| 21 | 3.04, 3.12, 3.09 |
| 101 | 1.41, 1.39, 1.40 |
| 501 | 0.625, 0.629, 0.629 |

My earlier test used a **null matrix only**, where the p-value comes from the pre-permutation
screen and the multilevel estimator never escalates, so `sampleSize` cannot show an effect.
The identical five-seed series I reported at `sampleSize = 21` and `101` was real and was
evidence of nothing.

The owner's pointer is what prompted the recheck: GESECA is a first-class fgsea method —
`R/geseca-multilevel.R`, `R/geseca-simple.R`, `R/geseca-utils.R`, `R/geseca-plot.R`,
`src/geseca.cpp`, plus a tutorial vignette. `gesecaCpp` is one internal helper inside a fully
public method, and the vendored kernel reached past the front door for convenience.

**Agreement, measured properly.** On real chunks at `sampleSize = 21`, `fgsea:::gesecaCpp`
and `fgsea::geseca()` agreed within their summed `log2err` in **7 of 8 datasets**;
`|log2(ratio)|` reached 2.55 against a median summed bound of 1.53. Two Monte-Carlo
estimators of one quantity agree within their stated error, not to equality — which is the
comparison §10 should have made instead of expecting identical values.

**What stands from §10 and §11, unaffected:**

- `pctVar` uses the stored `obj$totalVar` and needs no fgsea call. Verified bit-identical to
  hand computation on 12 of 12 real datasets.
- `alserglab/coresh` v0.1.0 exports nothing; C5's validator and documented R path are broken.
- `gesecaSimple()` floors p near `1/nperm` and cannot express CoReSh's magnitudes. It is
  still useful as an independent cross-check at moderate p.
- Chunk Entrez rownames are not unique. `make.unique()` before calling `geseca()`; index by
  `match()` everywhere else.

**What changes.** §5's option 2 is dead rather than deferred; **no `:::` enters the package
and no owner decision is owed.** `coresh_match(pvalues = TRUE)` routes to `fgsea::geseca()`
with `center = FALSE`, `scale = FALSE`, uniquified rownames and `log2err` carried into the
result. The guard that currently stops that path, and the message naming three dead ends, are
both removed.

The work is scheduled as **G1** in
[../2026-08-13-analysis-api-roadmap/00_ROADMAP.md](../2026-08-13-analysis-api-roadmap/00_ROADMAP.md),
§3, which also opens `gs_coregulation()` as a first-class verb — GESECA takes a matrix and no
contrast, so it complements DE-driven GSEA rather than competing with it.


---

## 13. C3 completed (2026-08-18) — the set-building layer

`coresh_loadings()` and `coresh_sets()` landed as roadmap step G2. See
[../2026-08-13-analysis-api-roadmap/00_ROADMAP.md](../2026-08-13-analysis-api-roadmap/00_ROADMAP.md)
§10 for the full record. In brief:

- The ported math is **identical to the reference on 63 of 63 comparable real hits** — same top-50
  Entrez ids, same order.
- §2's latent defects are discharged: no `options()` memoization, no `digest` dependency, no
  global-environment guard, and provenance is a table rather than a re-parse of set names, so
  `1.11:453`'s `max(which(grepl("^GSE", ...)))` returning `-Inf` has no analogue here.
- The Jaccard rule is now explicit — keep the higher-ranked hit — where the reference dropped
  whichever set came later in the input.
- `gs_db()` gained `provenance` and `set_provenance`, both surviving subsetting.

**Two findings that change what C4 has to do.**

`gsdb_coresh()` cannot be gated against DC-nexus's stored GMT the way §C4 assumes. That GMT dates
from 2026-04-25 and the earliest surviving chunk snapshot is the `_migrated` rewrite of 2026-04-30;
28 of 58 sets still match exactly, but the 24 that differ are explained by chunk content that no
longer exists, verified four ways rather than assumed. Its `coresh_provenance.csv` cannot arbitrate
either: it holds 13 distinct `pctVar` values across 58 rows, having recorded each query's score
rather than each hit's. **C4's gate has to be the loading-level agreement above plus a fresh
end-to-end run, not byte agreement with a stored artefact.**

And **a GSE accession is not a unique key**: 1,635 of 42,465 human accessions appear more than once
on different platforms, sometimes in different chunk files, and the first-match lookup is therefore
order-dependent. Fix that before `gsdb_coresh()` builds on it.
