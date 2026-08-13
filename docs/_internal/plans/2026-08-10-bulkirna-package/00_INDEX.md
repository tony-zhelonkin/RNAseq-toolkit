# RNAseq-toolkit → `bulkiRNA` package — working plan and state

**Repo path:** `/data1/users/antonz/pipeline/bulkiRNA` (renamed from `RNAseq-toolkit`
during the 2026-08-10 session — a Phase 6 item, done early. The `hub` remote still points at
`/data1/users/antonz/git/RNAseq-toolkit.git`; `origin` is `git@github.com:tony-zhelonkin/bulkiRNA.git`.)
**Branch:** `feat/bulkirna-package` (off `dev`)
**Fork point / revert target:** `752481f` (`dev` tip at fork)
**Baseline at fork:** 9 legacy suites in `tests/`; 20 golden cases captured, 0 errors
**Current:** Steps 0, A, **B1–B6 and C all complete and merged.** 61 exports
(41 new API + 20 deprecation shims). `R CMD check --as-cran`: **0 errors, 0 warnings,
0 notes.** 700 tests (dev image) / 723 (gatom image), 0 failures. Golden **20/20** against
the package, with the perturbation self-test still failing it (19/1) so the gate is live.
**`scripts/` is deleted — the package is the only implementation.**

Step C, done: `R/deprecated-gs.R` + `R/deprecated-plot.R` shim all 20 frozen names that
needed one (4 of the 24 kept their own names); the golden harness now loads the package
instead of sourcing `scripts/`; ten goldens re-captured with a written reason each; the
freeze verified mechanically at 24/24 before `scripts/` went. Ledger and evidence in
`STEP-C-REQUIREMENTS.md` §§ 0b and 2.

**Done: the architecture review pass, and the medium findings after it.** Six Opus reviewers
over `R/` by layer produced 1421 lines of findings — 8 high, ~20 medium. Every finding acted
on was re-verified by the integrator in the container first, and every fix carries a
regression test **proven to fail without it**. See §9 for the highs and §10 for the mediums.
Tests **802** pass / 0 fail (was 700 at merge). Golden **20/20**. `R CMD check`
(`_R_CHECK_FORCE_SUGGESTS_=false`, no `--as-cran`) **OK — 0/0/0**. Under `--as-cran` there is
one NOTE (CRAN wants the Title as "RNA-Seq") and one ERROR that is environmental only
(`gatom`/`mwcsr` are Suggests and are absent from the dev image). An earlier claim of
"`--as-cran` 0/0/0" in this file was the non-`--as-cran` run; corrected here.

Note for anyone verifying against the old code: `scripts/` is deleted, but readable at
`ff80de2^` — `git show ff80de2^:scripts/<path>`.

**Then: Step 1b / Phase 6** — unfreeze the 24 names, rename where the design says to, update
consumers. Nothing blocks it. Note that "unfreezing" means retiring the shims, which breaks
the 64 legacy call sites, so it is coupled to Phase 4 consumer migration in `sciagent-rna`
and is a scope decision for the user, not a refactor decision.

**Pushed to `hub` only** (through `9b63f91`; four commits since are local). `origin` needs an ssh-agent with a passphrase key —
ask, do not work around it. Merged worktrees have been removed; the `wt/*` branches are kept
as provenance.

Resume by reading this header, then §1 for how the state was reached. **§3 and the two blocks
at the end of §1 are historical snapshots, explicitly labelled — do not act on them.**

---

## 0. The decision, in one paragraph

This script library becomes an installed R package named **`bulkiRNA`**. The GitHub repo is
already renamed (`tony-zhelonkin/bulkiRNA`); local paths and the hub still carry the old
name and move in Phase 6. Scope is deliberately narrow: **fgsea + GSVA + ORA only** — no
engine registry, no EnrichmentBrowser, and `clusterProfiler`/`enrichplot` leave entirely.
The public surface drops from 123 definitions to ~30 exports behind a four-layer namespace
(`gsdb_*` providers → `gs_test`/`gs_score` → `gs_result`/`gs_matrix` → `gs_plot_*`).

Design record lives **outside this repo** in `/data1/users/antonz/pipeline/sciagent-rna/docs/`
(hub: `/data1/users/antonz/git/sciagent-rna.git`):

| Doc | What |
|---|---|
| `00` | Overview: why, the three layers, the backdoor, phase status |
| `02` | API inventory + the **frozen 24-export list** + collision resolutions |
| `03` | Blast radius: 64 call sites, 24 files, 10 projects |
| `04` | Dependency ledger + reference-data layout |
| `05` | Namespace/rename plan; the near-duplicate families to collapse |
| `07` | **API design — the contract every agent codes against** |
| `08_refactor-execution-plan` | The multi-agent pass: Step 0/A/B1–B5/C, per-agent briefs |

**B6 (`gatom_*`) is not in those docs** — it was agreed in session on 2026-08-10 and is
specified in §5a of this file. `07` §7's "~30 exports" and `08`'s "five parallel agents"
both predate it.

Those are the record of *why* and *what shape*. **This file is the source of truth for
_what is done_ and _what is next_.**

---

## 1. State

### Step 0 — golden safety net ✅ done (`7ae41eb`, `1fef700`, `0c80549`)

`tests/fixtures/` — deterministic synthetic dataset (2000 genes × 8 samples, seed
`20260807`): planted up/down/**null** sets, cross-database overlap, 20 duplicated symbols,
80 scattered DE genes. Plus a **real-symbol companion** (4393-symbol MSigDB mouse universe,
signal in two real Hallmark sets) because `run_gsea()`/`load_reference_db()` resolve by real
symbol and would otherwise return empty.

`tests/golden/` — **20 cases captured, 0 errors.** `capture_golden.R` owns the case
definitions and runs in capture or `--verify` mode; `verify_golden.R` is the `--verify`
wrapper and exits non-zero on drift. Plots are stored as `ggplot_build(p)$data`, not pixels.

**The harness self-tests.** `GOLDEN_SELFTEST_PERTURB=1` wraps `format_pathway_name` in
`toupper()`; verify then reports 18 pass / 2 fail — catching the direct change *and* its
propagation into `normalize_gsea_results$pathway_name`. Clean exits 0, perturbed exits 1.
A gate nobody has watched fail is not a gate.

### Step A — skeleton and shared contracts ✅ done

`DESCRIPTION` (`bulkiRNA` 0.3.0, the `07` §8 ledger + `stats`), `LICENSE`,
`.Rbuildignore`, `NAMESPACE` (roxygen), `man/` (28 topics), `CONVENTIONS.md`.

- `R/bulkiRNA-package.R` — package doc; **the one `@import ggplot2`**.
- `R/utils.R` — `ensure_dir()` (directory semantics, exported), `ensure_parent_dir()`
  (file semantics, internal — the `02` §3.2 collision resolved by splitting, not renaming
  in place), `%||%`. Nothing else, ever, without a handback request.
- `R/gs-result.R` — `gs_result()` constructor / `new_gs_result()` / `validate_gs_result()`,
  the 12 core columns of `07` §4 with enforced types, `gs_stat_types()` (the only export
  besides `ensure_dir`), `gs_direction()`, `gs_stat_label()`, and S3 `print` (via
  `tbl_sum`), `summary`, `as_tibble`, `rbind`, `[`, `dplyr_reconstruct`.
- `R/gs-matrix.R` — `gs_matrix()` + validator, five metadata attributes, accessors, and
  S3 `print`, `summary`, `[`, `as_tibble` (long form joined to `sample_data`).
- `inst/extdata/` — the 152 KB processed tree + `METADATA.yaml` + three `CITATIONS.bib`.
- `tests/testthat/` — 67 tests across `test-gs-result.R`, `test-gs-matrix.R`,
  `test-utils.R`, plus `helper-gs.R` builders. All pass.

`R CMD check --no-manual`: **0 errors, 0 warnings, 1 NOTE** — unused `Imports` (`fgsea`,
`msigdbr`, `ggrepel`, `scales`, `stringr`, `tidyr`, `rlang`), which B1–B5 consume. Golden
verify still 20/20.

Constructors are **internal** (`@keywords internal`, no `@export`): `07` §7 lists no
`gs_result()` in the ~30, and B agents reach them package-internally. Both result classes
downgrade to a plain tibble when a verb drops a core column, so a broken object cannot
survive a `select()`.

### Step B1 — databases (`gsdb_*`) ✅ done, merged

Commits `bfe34e9` → `3f8eb8d` → `17481af` → `ce82974`, merged to base. **7 exports:**
`gsdb_list`, `gsdb_load`, `gsdb_register`, `gsdb_from_file`, `gsdb_msigdb`, `gsdb_info`,
`download_gatom_references`. Files: `R/gs-db.R`, `R/gsdb-load.R`, `R/gsdb-file.R`,
`R/gsdb-msigdb.R`, `R/gsdb-rebuild.R`, `R/gsdb-gatom.R`. 143 tests.

`.resolve_toolkit_dir()` and every `toolkit_dir`/`helper_root` argument are gone; paths
resolve via `system.file("extdata", ...)`. `build_reference_databases.R`'s source-time
executable body became the internal `.gsdb_rebuild()`, which errors "source checkout
required". `gsdb_from_file()` sniffs GMT vs GMX. Deleted with evidence: `parse_geneset_file`,
`convert_geneset_ids` (both verified absent from the frozen 24).

### Step B2 — compute (`gs_*`) ✅ done, merged

Commit `28dc231`, merged to base. **9 exports:** `gs_ranks`, `gs_test`, `gs_score`,
`gs_leading_edge`, `gs_filter`, `gs_top`, `gs_split`, `gs_write`, `gs_read`. Files:
`R/gs-test.R`, `R/gs-score.R`, `R/gs-ranks.R`, `R/gs-ops.R`, `R/gs-io.R`. 101 tests.
`gs_test()` is S3-dispatched on `numeric`/`integer`/`character`/`gs_matrix`/`default`.
**clusterProfiler and enrichplot appear nowhere in the package.**

**B2 never delivered a handback report** — it went idle three times. Its gates were verified
by the integrator instead, so treat its *reasoning* as unrecorded: what it deleted and why
is not written down anywhere. See §2 finding 8.

### ✅ THE EQUIVALENCE RESULT — the one that justifies the whole migration

`docs/_internal/plans/2026-08-10-bulkirna-package/B2-fgsea-equivalence.R`, run on the
real-symbol fixture, old `run_gsea()` (clusterProfiler) vs new `gs_test()` (direct fgsea):

```
OLD run_gsea(): 50 pathways    NEW gs_test(): 50 pathways
Pathways only in OLD: 0        Pathways only in NEW: 0
max |dNES|  = 0.000e+00
max |dp|    = 0.000e+00
max |dpadj| = 0.000e+00
```

**Bit-identical**, down to `p = 6.220419e-131`. Re-run it after any change to `.gs_fgsea()`.
It is a diagnostic, not a package test; it sources `scripts/`, so it dies with `scripts/` in
Step C — capture its output in the Step C commit message before deleting it.

### Cross-module integration — verified live by the integrator

B1 and B2 never saw each other's code. Confirmed working on the merged branch:
`gsdb_load("mitopathways")` → `gs_test()` = 39 rows, `database = "mitopathways"`,
`stat_type = "NES"`; `gs_test(ranks, list(hallmark=, mito=))` = 89 rows over two databases;
`rbind()` pooling preserves class; character-vector input → ORA with
`stat_type = "log2_fold_enrichment"`; both `gs_filter(res, padj = 0.05)` and
`dplyr::filter(res, padj < 0.05)` work and keep the class. Zero-row `gs_result` constructs
(12 cols) and `rbind`s — so `empty_gsea_tibble()` is shimmable.

### Gate on the merged branch (`633dddc`)

**311 tests, 0 failures** · **golden 20/20, exit 0** · `R CMD check --no-manual`
**0 errors, 0 warnings, 1 NOTE** (unused Imports `ggrepel`, `rlang`, `scales`, `stringr`,
`tidyr` — all claimed by B3/B5) · **18 exports** so far.

### Steps B3, B4, B5, B6, C — ✅ all done and merged

**The two blocks above are a historical snapshot, kept deliberately.** They record the state
at `633dddc` — 311 tests, 18 exports, `scripts/` still present — because the sentence
"`scripts/` is untouched, which is why golden still passes" is the single most important
thing this document says. It was true, and it meant the golden gate was **green and
structurally blind**: it exercised only the old code. Migrating it onto the package is what
found the five defects in §8. Read those two blocks as evidence, not as current state.

Current state is in the header: 61 exports, 700/723 tests, golden 20/20 against the
**package**, `scripts/` deleted.

---

## 2. Findings that change the work (discovered during Step 0)

1. **`build_reference_databases.R` executes at source time** and resolves paths from its own
   location — it cannot be sourced as a definition file. `filter_by_size` had to be reached
   by parsing function definitions out of the AST. Concrete case for splitting it → **B1**.
2. **`.resolve_toolkit_dir()` auto-detects via `sys.frame()$ofile`**, which only works
   *during* sourcing and fails when the function is called later. Callers must pass
   `toolkit_dir` explicitly today. `system.file("extdata", ...)` deletes this → **B1**.
3. **msigdbr 26.1.0 downloads MSigDB at runtime** to `$HOME/.cache/R/msigdbr` — no longer
   bundled. Needs `-e HOME=<writable>` + a mounted cache, and network on a cold cache →
   affects **B2** and every agent touching `run_gsea`.
4. **`clusterProfiler::GSEA` is already `by = "fgsea"`** (`run_gsea.R:181`). Going direct to
   `fgsea::fgseaMultilevel()` is the same computation. **The golden diff on `run_gsea` must
   be empty**; a non-empty diff is a stop signal, not a tolerance to widen.
5. **`fgsea::plotEnrichmentData()` returns `curve`/`ticks`/`stats`** — exactly the three
   panels the running-sum plot needs. **B4 does not hand-roll the cumulative sum.**
6. **`fora()` returns `foldEnrichment`**, not an odds ratio → `stat_type =
   "log2_fold_enrichment"`.
7. **The image gap is closed**: `scdock-r-dev:v0.5.11` has all 10 Imports and all 14
   Suggests, including `org.*.eg.db`, `FactoMineR`, `factoextra`, `homologene`. The pending
   fleet action is a `provc` bump from `v0.5.10` → `v0.5.11`, not a rebuild.

### Discovered during Step A / B1 / B2 (2026-08-10)

8. **Two design docs gave contradictory orders, and every B brief cited only one.**
   `02_api-inventory.md` §5 freezes 24 exports; `07_api-design.md` §7 lists the ~30 new
   ones; they are different sets. B1 read §7 plus "deletion is the job" and deleted
   `download_gatom_references()`, a frozen export. **Ruling, now in `CONVENTIONS.md` §12 and
   `08` §5.1:** a frozen name may be renamed or internalised (Step C shims it) but never
   deleted without leaving a body the shim can call. The test is not "does anything call it"
   but "can Step C still make the old call work". **B3–B5 face the identical trap.**
9. **`data/used-functions.tsv` is not in this repo.** It is at
   `/data1/users/antonz/pipeline/sciagent-rna/docs/data/used-functions.tsv`. Earlier drafts
   pointed agents at a path that does not exist; that caused finding 8.
10. **`08` §4's `{t2g, t2n}` provider sentence was stale** — the clusterProfiler input
    format. Corrected in `21e0655` (sciagent-rna). `CONVENTIONS.md` §11a governs.
11. **`gs_result$database` is the stable snake_case registry key**, not a display label —
    it is a join and filter key. `database_label` carries the display string.
    **B3 must read `database_label` for facets and legends**, never `database`.
12. **Process defect (integrator's).** `CONVENTIONS.md` and `DESCRIPTION` were edited on the
    base branch *after* the worktrees existed, so both agents worked from a stale contract
    and B1 needed a merge mid-flight. **For B3–B5: freeze the shared contract for the batch,
    or message all agents simultaneously with a merge instruction.**
13. **Never measure a live worktree.** Test counts and three "failures" observed mid-edit
    were artifacts of sampling a worktree while its agent was still writing. Verify only
    after the agent commits and reports idle.
14. **`parse_transportdb()` had a dead `org_db` argument** and `METADATA.yaml` disagrees with
    the parser about the TransportDB raw format (METADATA: headerless 7-column CSV needing
    AnnotationDbi; parser: `read.csv(header = TRUE)` hunting `Symbol`/`Family`). The shipped
    RDS contains mouse symbols, so it was not built by this path. **Pre-existing latent
    rebuild bug on a path with no golden baseline.** Recorded, not chased.
15. **Rebuild size-filter ordering changed**: the old code filtered 5–500 inside `parse_gmx()`
    *before* human→mouse conversion and again after; `.gsdb_rebuild()` filters only after.
    Shipped RDS untouched so golden is unaffected, but a future rebuild could differ.
16. **`mito_unified`'s RDS carries an unused `merge_map`** element. Dropped, no contract slot.
    Recoverable: `readRDS(system.file("extdata",
    "mitochondria_unified/processed/Mus_musculus/unified_mito_pathways.rds",
    package = "bulkiRNA"))$merge_map`.
17. **`readxl`, `org.Mm.eg.db`, `org.Hs.eg.db` are in `Suggests` but used by nothing** so far.
    If B5 does not claim them, remove them in Step C.
18. **`msigdbr` prints a once-per-session ortholog notice** that is not suppressible via its
    API. Documented exception to `CONVENTIONS.md` §5. **Do not write a test asserting
    MSigDB providers are silent.**

---

## 3. HISTORICAL — the dispatch brief for B3, B4, B5, B6

**All four shipped and merged; this section is retained as the record of what they were told,
not as an instruction.** It is worth keeping because §8's assessment of the fleet only makes
sense against the briefs they actually received. Do not act on it. For current next steps see
the header.

Briefs are `08_refactor-execution-plan.md` §4; the rules are §5 + **§5.1 (the freeze rule)**,
and the in-repo contract is **`CONVENTIONS.md`** — agents read that first.

### Setup

```bash
cd /data1/users/antonz/pipeline/bulkiRNA
git worktree add -b wt/b3-renderers ../.bulkirna-wt/b3 feat/bulkirna-package
git worktree add -b wt/b4-running   ../.bulkirna-wt/b4 feat/bulkirna-package
git worktree add -b wt/b5-de-io     ../.bulkirna-wt/b5 feat/bulkirna-package
git worktree add -b wt/b6-gatom     ../.bulkirna-wt/b6 feat/bulkirna-package
```

`../.bulkirna-wt/b1` and `b2` still exist and are fully merged — remove them with
`git worktree remove` once you are confident, or leave them as reference.

**Freeze `CONVENTIONS.md` and `DESCRIPTION` for the duration of the batch** (finding 12). If
one genuinely must change, message all three agents simultaneously with a merge instruction.

### What each agent must be told beyond its §4 brief

All three:
- **The freeze rule (§5.1 / `CONVENTIONS.md` §12).** State it explicitly — B1 got it wrong.
  Give each agent the frozen names *in its own module* and require a per-name statement of
  how Step C reaches it.
- `@import ggplot2` already lives in `R/bulkiRNA-package.R`. **Do not add a second one**, and
  no `ggplot2::` prefixes are needed.
- `@export` on an S3 method whose generic belongs to another package emits `export()`, not
  `S3method()`. Use `@exportS3Method pkg::generic`. **Check `NAMESPACE` after `document()`.**
- ggplot2 is **4.0.3**: `colour = "transparent"`, never `colour = NA`.
- Never edit `DESCRIPTION`, `NAMESPACE`, `CONVENTIONS.md`, `R/utils.R`, `R/gs-result.R`,
  `R/gs-matrix.R`, `R/bulkiRNA-package.R`, `R/gs-db.R`, or anything B1/B2 own. Report instead.
- **Do not modify `scripts/`.** Step C deletes it; golden depends on it.
- Deliver a handback report. B2 did not, and its reasoning is lost.

**B3 (renderers + theme)** — frozen names: `gsea_dotplot`, `gsea_dotplot_facet`,
`gsea_barplot`, `format_pathway_name`, `custom_minimal_theme_with_grid`, `save_gsea_log`,
`plot_all_gsea_results`. Exports `gs_plot_dot`, `gs_plot_bar`, `gs_plot_heatmap`, `gs_save`,
`theme_bulki`, `format_pathway_name`. **Read `database_label`, not `database`** (finding 11).
Axis labels come from `gs_stat_label(res)` — never a literal `"NES"`. Claims `ggrepel`,
`scales`, `stringr` from the unused-Imports NOTE.

**B4 (running-sum rewrite)** — frozen name: `gsea_running_sum_plot`. Exports
`gs_plot_running`. **`gs_leading_edge()` and `gs_ranks()` now exist for real** — build on
them rather than the contract. `fgsea::plotEnrichmentData()` returns `curve`/`ticks`/`stats`,
the exact three panels; **do not hand-roll the cumulative sum.** Colours keyed by pathway id
via explicit `scale_colour_manual(values = named_vector)`, never positional. Read the old
file's `@note` block (`scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R:95-101`) — it
documents the trap not to reproduce. This is the one deliberate golden change.

**B5 (DE + IO)** — frozen names: `create_standard_volcano`, `create_MD_plot`, `build_dge`,
`ensure_dir` (already exported by Step A — do not redefine). Exports `de_volcano`,
`de_md_plot`, `de_bfc_plot`, `de_pca`, `de_pca_3d`, `de_volcano_grid`, `build_dge`,
`annotate_genes`, `read_counts_matrix`, `read_metadata`, `write_session_provenance`.
**Preserve decision-by-FDR volcano semantics exactly** — dashed line at the p-value
corresponding to the FDR boundary, plus `fixed_p_boundary`. Golden-tested, load-bearing.
Delete the byte-identical duplicate `aggregate_duplicate_ids` (`io_helpers.R:231` vs `:177`),
`source_if_present()` and the `here` dependency; rename `save_plot` → internal
`famd_save_plot`. Claims `tidyr`; guard `plotly`, `FactoMineR`, `factoextra`.
**`FactoMineR`/`factoextra` are not in DESCRIPTION** — request them if FAMD stays.

### Then Step C

1. Merge b3/b4/b5. **`NAMESPACE` will conflict — resolve by regenerating**, never by hand
   (this worked cleanly for b1+b2).
2. Apply collected `DESCRIPTION`/`R/utils.R` requests.
3. `devtools::document()`; confirm exports match `07` §7 plus B6's five (~35 total;
   currently 18). **`07` §7's "~30" predates the B6 decision — update that doc in Step C.**
4. **`R/deprecated.R`** — the 24 frozen names → new names with `.Deprecated()`. Each B
   agent's handback states how its names are reached. Known shim details:
   - `list_reference_dbs()`: `gsdb_list()` then `x$species[!x$bundled] <- "(not bundled)"` and
     reorder to `database, name, bundled, description, species`. Accept and ignore
     `toolkit_dir`.
   - `load_reference_db()`: `gsdb_load()` → `.gsdb_as_t2g()`, plus re-add the old `source` and
     `created` list elements.
   - `filter_by_size(result, min_size = 5, max_size = 500)`: `.gsdb_from_t2g()` → internal
     `filter_by_size()` → `.gsdb_as_t2g()`. **The 5/500 defaults live in the shim**; the
     internal version defaults to `NULL`/`NULL`.
   - `empty_gsea_tibble()`: a zero-row `gs_result` — verified constructible.
5. `R CMD check --as-cran`.
6. **`verify_golden.R` — the real gate.** Then migrate the golden harness off `scripts/` onto
   the new API, and only then delete `scripts/` (the user's explicit sequencing).
7. Capture the equivalence-script output in the commit message before deleting it.

---

## 4. Environment — no R on the host

```bash
DK() { docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
  -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
  -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 "$@"; }

DK Rscript tests/golden/verify_golden.R          # the Step C gate
DK Rscript tests/fixtures/make_fixture.R         # regenerate fixture (deterministic)
DK Rscript -e 'devtools::load_all(); testthat::test_local()'
```

Both `--user` and `HOME` are mandatory; without them `saveRDS` and msigdbr's cache
respectively fail. R 4.5.3 · ggplot2 **4.0.3** (so `colour = "transparent"`, never
`colour = NA`) · fgsea 1.36.2 · msigdbr 26.1.0 · GSVA 2.4.9 (param-object API).

---

## 5. Expected golden changes

Record differences here; do not loosen the comparison.

| Case | Expected | Why |
|---|---|---|
| ORA | will differ | `enrichGO`/`enrichKEGG` → `fgsea::fora`; GO/KEGG now from MSigDB |
| `gsea_running_sum_plot` | complete redesign | enrichplot dropped; B4 is from scratch |
| all | names change | `gs_*`/`de_*`; `R/deprecated.R` keeps old names callable |
| `run_gsea` | **must not change** | already fgsea underneath — a diff here means stop |

---

## 5a. GATOM — IN SCOPE, as agent B6 in the B3–B5 batch

**Corrected 2026-08-10.** An earlier revision of this section called the `gatom_*` module
"deferred by the user's decision". **That was wrong on both counts** — it was the
integrator's recommendation, not the user's instruction, and the user's actual position is
that GATOM is a routine part of their workflow and belongs here. It is **not deferred.**

The deferral rested on a false premise: that `gatom` would need a new dependency baked into
an image before anything could be verified. It does not. **`scbio-singleuser:v1.7.1` already
has `gatom` 1.8.4, `mwcsr` 0.1.11, `igraph` 2.3.1, `devtools`, `testthat`, the full bulkiRNA
Imports at identical versions, and the staged references at `/opt/gatom-refs/`**
(`network.kegg.rds`, `met.kegg.db.rds`, `org.Hs.eg.gatom.anno.rds`). So B6 is testable today
in a second container; `scdock-r-dev:v0.5.11` lacks `gatom`/`mwcsr`, and package tests
`skip_if_not_installed("gatom")` there.

```bash
# B6's test image — note HOME=/tmp, and no msigdb cache mount needed
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PWD":/pkg -w /pkg scbio-singleuser:v1.7.1 <cmd>
```

`download_gatom_references()` is **already done and exported** (`3f8eb8d`) — that part was
never deferred and is not B6's to redo. B6 absorbs it into the wider surface.

### Where the presets actually live

`/data1/users/antonz/pipeline/scbio-instruct/04_provision/seed_content/06_gatom/06_network_handoff.qmd`
(291 lines) is the working pipeline and the real spec. Also read
`scbio-instruct/docs/ai-research/gatom-validation.md`. The SciAgent skill
(`SciAgent-toolkit/skills/gatom-metabolomic-predictions/SKILL.md`) is **prose only, no code**
— it re-teaches the raw API. B6's job is to turn its four documented traps into enforced
invariants, after which SKILL.md shrinks to a pointer:

1. `pval` must be **raw**, never `padj` — BUM scoring breaks silently otherwise.
2. `baseMean` must be **linear** scale (`2^AveExpr`), not log.
3. In `topology = "atoms"` graphs, **genes live on EDGES**, not vertices —
   `igraph::as_data_frame(m, "edges")$Symbol`.
4. `met.db` is **required even when `met.de = NULL`**.

Plus two operational traps from the qmd: `saveModuleToHtml()` needs pandoc on PATH
(`Sys.setenv(RSTUDIO_PANDOC = ...)`), and `k.gene` is the module-size dial (smaller k →
larger module; 50 is the gatom default, 25/75 are the sensitivity branches).

### B6 — the agreed surface (user-chosen, 2026-08-10): 5 exports

```r
refs <- gatom_refs(species = "Homo sapiens")
de   <- gatom_de(tt, id = gene, pval = p.value, log2FC = abs_zscore, baseMean = 1)
m    <- gatom_module(de, refs, k_gene = 50, seed = 42)
gatom_genes(m)
gatom_save_html(m, "kyn_module.html", name = "Kynurenine")
```

| Export | Does | Encodes which trap |
|---|---|---|
| `gatom_refs(species, dir = NULL, download = FALSE)` | Loads + validates the three reference files as one object. Searches `/opt/gatom-refs`, then `download_gatom_references()`'s `dest_dir`, then `dir`. `download = TRUE` delegates to `download_gatom_references()` — does not reimplement it. Errors naming the missing file and how to fetch it. | 4 — carries `met.db` so it can never be omitted |
| `gatom_de(x, id, pval, log2FC, baseMean)` | Builds `gene.de`: `arrange(pval)`, `distinct(ID)`. **Validates:** errors if the `pval` column name matches `adj\|fdr\|padj\|q.?val`; errors if any `pval` is outside `[0,1]`; **warns** if `baseMean` has negatives or `max(baseMean) < 30` (log-scale smell); errors on all-`NA`. | 1 and 2 |
| `gatom_module(de, refs, k_gene = 50, k_met = NULL, met_de = NULL, seed = 42, solver = "rnc")` | `makeMetabolicGraph(topology = "atoms", keepReactionsWithoutEnzymes = FALSE)` → `scoreGraph()` → `solve_mwcsp()`. Returns the module igraph with `k_gene`/`seed`/`solver`/node+edge counts as attributes for provenance. | 4 |
| `gatom_genes(m)` | `unique(igraph::as_data_frame(m, "edges")$Symbol)` | **3 — the big one** |
| `gatom_save_html(m, path, name)` | Wraps `saveModuleToHtml()`, sets `RSTUDIO_PANDOC` if unset, `ensure_dir()`s the parent. The only B6 function that writes. | pandoc-on-PATH |

Constraints for B6:
- `gatom`, `mwcsr`, `igraph` are **`Suggests`** (added to DESCRIPTION `2026-08-10`, pre-dispatch).
  Every entry point guards with `requireNamespace()` and an actionable install message.
- Tests: `skip_if_not_installed("gatom")` so the suite passes in **both** images. Validation
  tests for `gatom_de()` need no gatom at all — **those must run unskipped everywhere**, since
  they encode the traps.
- **No golden baseline exists for GATOM** (`download_gatom_references` is in the harness's
  skip list). B6 is new surface, so nothing to preserve — but it must not perturb the 20
  existing cases.
- `set.seed()` before `solve_mwcsp()`: the solver is a heuristic. Pin it and test that the
  same seed gives the same module size.
- Only human refs are staged (`org.Hs.eg.gatom.anno.rds`). Mouse must fail with a clear
  "run download_gatom_references(species = \"Mus_musculus\")" message, not a missing-file error.
- Do **not** touch `R/gsdb-gatom.R`'s frozen `download_gatom_references()` formals.

**Superseded text removed 2026-08-10.** This section previously ended with a "Why not now"
paragraph and an earlier four-function sketch (`gatom_graph()`, `gatom_plot()`). Both are
dead: the surface is the five functions in the table above, and the deferral reasoning was
withdrawn (the `Suggests`-ledger objection was resolved pre-dispatch; the no-golden-baseline
point survives only as the constraint already listed above). Ignore any older copy.

`download_gatom_references()` stays **exported with its exact current name and
formals** — it is on the frozen 24-export list, and `02_api-inventory.md` §5 says to export
it regardless of having no in-repo caller. B1 deleted it on a "no consumer" reading; that
test does not override the freeze. Restored.

---

## 6. Revert

```bash
git checkout dev && git branch -D feat/bulkirna-package
git worktree remove ../.bulkirna-wt/b1 && git branch -D wt/b1-databases
git worktree remove ../.bulkirna-wt/b2 && git branch -D wt/b2-compute
```

`scripts/` is still unmodified, so abandoning this branch costs the package skeleton, the
`gsdb_*`/`gs_*` layers, the fixture and the golden harness — the last two are worth keeping
regardless. **Nothing has been pushed**, so revert is purely local.

---

## 7. Kickstart prompt for a fresh session

```
Resume the bulkiRNA packaging work. cd /data1/users/antonz/pipeline/bulkiRNA

Read, in order:
  1. docs/_internal/plans/2026-08-10-bulkirna-package/00_INDEX.md   (state + next action)
  2. CONVENTIONS.md                                                (the in-repo contract)
  3. /data1/users/antonz/pipeline/sciagent-rna/docs/07_api-design.md §7 (export list)
  4. /data1/users/antonz/pipeline/sciagent-rna/docs/08_refactor-execution-plan.md §4 (B3/B4/B5
     briefs), §5 and §5.1 (rules + the freeze rule)

State: branch `feat/bulkirna-package` at `fc016ff`, pushed to `hub` only. Steps 0, A,
B1–B6 and C are ALL done and merged — 61 exports (41 new API + 20 deprecation shims),
700/723 tests passing, golden 20/20 against the PACKAGE, R CMD check 0E/0W/0N.
`scripts/` is DELETED; the package is the only implementation. Old code is still
readable at `ff80de2^` (`git show ff80de2^:scripts/<path>`) for verification.
`v0.3.0` is deliberately untagged pending the architecture review pass.

Task: see the header of 00_INDEX.md for what is actually in flight. As of the last
session that was the six-reviewer architecture pass (read-only reviewers, integrator
applies fixes). After it: Step 1b / Phase 6, which is coupled to consumer migration
and is the user's scope call.

Hard constraints:
- No R on this host. Everything runs in a throwaway container:
    docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
      -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
      -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 <cmd>
  Both --user and HOME are mandatory (saveRDS permissions; msigdbr's runtime cache).
- NAMESPACE is roxygen-generated, never hand-edited. On merge conflict, REGENERATE.
- The freeze rule (CONVENTIONS.md §12): a name on the frozen 24 may be renamed or made
  internal, but never deleted without leaving a body Step C's shim can call. B1 got
  this wrong once; tell every agent explicitly.
- FREEZE CONVENTIONS.md and DESCRIPTION while agents are running. Editing the base
  branch under a live worktree cost real time this session.
- Never measure a live worktree — verify only after the agent commits AND reports idle.
- `Rscript tests/golden/verify_golden.R` must exit 0. Confirm at the start too, so a
  later failure is attributable. The gate now loads the PACKAGE; re-capture only with
  `--cases=<name>` so blessing one baseline cannot silently bless unrelated drift.
- The golden gate is BLIND to label text and theme element sizes (ggplot_build carries
  only numeric layer data). Assertions about those belong in tests/testthat/.
- tests/golden/ and tests/fixtures/ are integrator-owned; agents may read, never write.
- Push to `hub` and `origin` only when the user says so (origin needs an ssh-agent with
  a passphrase key — ask rather than working around it).

Verify agents' claims yourself; do not take a handback at face value. B2 never reported
at all and its gates had to be checked independently. Two Opus reviewers with tight
briefs returned nothing across four idle cycles — do not wait on a silent reviewer;
run the mechanical check yourself. See §8.
```

To resume the **SciAgent** track instead, point a session at
`scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/plans/2026-08-07-demolish-role-layer/00_INDEX.md`
(branch `refactor/demolish-role-layer`) and read §1 then §3.

---

## 8. Step C outcome, and the one lesson that generalises

**The golden gate was green and blind for the entire refactor.** `capture_golden.R` sourced
`scripts/`, so 20/20 PASS meant "the old code still works" -- it could not observe a single
thing the new code did. Migrating it onto the package (`9c7e1a5`) found five defects in one
sitting, none of which `document()`, `test()`, `check --as-cran` or the previous green golden
could see:

1. **The dotplot plotted the wrong variable.** Old x is GeneRatio (`count / setSize`); the
   shims left `gs_plot_dot()`'s `aes_x = "stat"` default, so the toolkit's most-used GSEA
   figure silently changed what its x axis *meant*. The capability already existed --
   `aes_x = "gene_ratio"` -- and the shim simply never asked. After the fix, x matches on
   14/15 pathways and all 15 dot sizes match exactly; the 15th is an exact `padj` tie
   (`HALLMARK_MYOGENESIS` vs `HALLMARK_TGF_BETA_SIGNALING`, both 0.2913391).
2. **`.de_theme()` discarded `base_size` entirely** (`base_size = 20` rendered at 14) and let
   `theme_bulki()`'s 14pt floor scale every volcano and MD plot by 14/12.
3. **A cross-agent seam**: `run_gsea()` returns a `gs_result`; the four `gsea_*` renderer
   shims accepted only S4 `gseaResult`. `obj <- run_gsea(...); gsea_dotplot(obj)` -- the
   commonest legacy call in the toolkit -- was broken end to end, and passed both agents'
   suites because each tested only its own half.
4. **`filter_by_size` was an unreachable export**, shadowed by collation order. Now guarded by
   `test-namespace-hygiene.R`, verified to fail on a planted duplicate.
5. **The golden baseline itself held a bug**: `ensure_dir`'s stored value came from the loser
   of a documented name collision, because the harness evaluated that definition last.

Three volcano cases and `create_MD_plot` pass **byte-identically**, which is the first hard
evidence `de_volcano()`/`de_md_plot()` reproduce their frozen predecessors rather than an
agent asserting they do.

### On the fleet

Seven agents. Three (C1, C2, C3) delivered clean work and, more usefully, pushed back: C1
found the `filter_by_size` collision and correctly refused to fix a file it did not own; C2
found the `gs_source` row-order divergence that is invisible in the figure but wrong in the
table `gs_save()` writes; C3 flagged two doc contradictions outside its scope rather than
silently widening its diff. **Two Opus reviewers returned nothing at all**, going idle without
reporting even when asked directly. Their checks -- the 24-name freeze diff, the namespace
scan, the golden triage -- were done by the integrator instead and found more than the reviews
did.

The generalisable lesson, consistent with findings 21 and 24: **running the gate found real
defects; asking an agent to reason about the code did not.** Prefer a check that executes.

### Two silent breaks for downstream code

- `direction` is `"up"`/`"down"`, was `"Up"`/`"Down"`. Filtering `direction == "Up"` yields
  **zero rows, not an error**.
- `empty_gsea_tibble()` and `run_gsea()` return `gs_result`'s 12 core columns; `$NES` on the
  result is now `NULL` rather than a vector.

---

## 9. The architecture review pass (2026-08-10)

Six Opus reviewers, read-only, one per layer: gsdb providers, compute/objects, `gs_plot_*`
renderers, DE, IO/utils/gatom, plus one cross-cutting architecture-and-shims remit. Brief:
namespace separation, organisation, readability, architecture, and behaviour against the old
implementations (still readable at `ff80de2^`).

### The eight high-severity findings, all confirmed and fixed

| # | Defect | Consequence |
|---|---|---|
| 1 | `.gs_empty_plot()` built on `ggplot()` with no data, so `$data` was a `waiver` and `nrow()` returned `NULL` | `gsea_barplot()` **errored outright** whenever nothing passed `padj_cutoff` — the normal outcome at the default 0.05 |
| 2 | `padj == 0` → `-log10` → `Inf`, unmappable by a size scale | The **strongest hit** silently absent from the dotplot, no warning |
| 3 | Two pathways formatting to the same label | Collapsed onto **one axis row**, bars stacked, figure shows n−1 of n |
| 4 | `de_volcano_grid()` read the caption via `ggplot_build(p)$layout$plot`, removed in ggplot2 4.0.3 | Threshold caption **lost in both modes** — the only place the realised raw-p boundary is stated |
| 5 | `de_bfc_plot()` appended highlights with `rbind()`, which uniquifies rownames | Printed **invented gene names** (`Gene11`, `Gene51`); the dedup that followed could never fire |
| 6 | `.io_guess_sep()` sniffed line 1, `read.delim(comment.char = "")` | Raw `featureCounts` output **unreadable**, erroring that a present `Geneid` column was absent |
| 7 | `ensure_dir()` discarded `dir.create()`'s result and never checked `dir.exists()` | **Reported success on failure**; every writer trusted it |
| 8 | Six shims named a private successor, one named a function that no longer existed | Deprecation advice **a user cannot act on** |

### What made this pass work, and what nearly broke it

- **Requiring an on-disk artefact.** All six reviewers went idle *without a final message* —
  the exact failure mode that produced nothing in the previous session. The briefs required
  findings written to `/tmp/review-*.md`, so all 1421 lines survived. **Never rely on a
  subagent's closing message as the deliverable.**
- **Re-verifying every high finding before fixing it.** Cheap, and it is the step that earns
  the right to change code on a reviewer's word.
- **Proving each regression test fails without its fix.** Done by reverting only the fixed
  hunk. One test (`de_volcano_grid`) initially appeared to fail for an unrelated reason — a
  stash had reverted a file to a version calling a helper deleted elsewhere — which is exactly
  the false signal that makes an unverified "it failed before" worthless.
- **Not accepting a fix that traded one bug for another.** The renderer agent fixed the
  raw-id-in-legend bug by formatting *every* label, and updated a pre-existing assertion to
  match `"Beta response"` → `"beta Response"`. `format_pathway_name()` is built for
  ALL_CAPS_SNAKE ids and is not idempotent on prose. Re-implemented to format only labels
  still equal to their own id. **An agent changing an existing assertion to match new output
  is a signal to look harder, not a completed task.**
- **The worktree base fault.** Every agent worktree was created at `25d9a6c`, predating the
  restructure, so `R/` did not exist. Two agents correctly refused to work around it and
  escalated. Fixing it repo-side was right; *telling the agents a constraint had been "worded
  too broadly" right after one was blocked by it* was not, and was correctly refused as
  permission laundering. If a brief's own wording blocks legitimate work, the fix is a new
  brief from the human, not a relaxation announced mid-flight.

### Declined, deliberately

- Barplot fill limits stay per-call rather than fixed at ±3.5 — cosmetic consistency, would
  move many baselines. Documented on `@param limits` instead, and `.gs_plot_all()` now passes
  one shared `limits` across the databases it renders, which was the case that actually
  misled a reader (see §10).
- `gs_source` row order ≠ draw order: no correctness impact, and reordering would fight
  `gsea_barplot()`'s deliberate re-sort for parity with the old function.
- `gsdb_list()`'s yaml-less fallback still omits the `gatom` row — a registry-design question,
  not a mechanical bug.

### One behaviour change for external callers

`gsdb_from_file(database = "x")` no longer sets `database_label` to `"x"`; the label defaults
to `basename(path)` and is set independently. A display label and a machine key are different
concepts. Every internal caller passes an explicit label, and golden stayed 20/20.

### Medium findings — all actioned, see §10

The six reviewer files are preserved in this directory as `review-*.md` (commit `9b63f91`),
so they no longer depend on `/tmp`.

---

## 10. The medium-findings pass (2026-08-10/11)

Same discipline as §9: verify, fix, and prove the regression test fails against the reverted
hunk (done in a throwaway copy under `/tmp`, never by stashing in the repo — see §9 for why).
Four commits, one per layer.

| Layer | Finding | What it did |
|---|---|---|
| compute | `set.seed()` never restored | Any random draw after a `gs_test()` call came from the seed-123 stream, not the script's own — once per database |
| compute | `names(db) %||% …` on a partly named list | Half the rows carried `database = ""` and grouped under a blank label |
| compute | `method` hard-coded `"gsva"` | ssgsea/zscore/plage runs recorded — and exported — as GSVA |
| compute | `[` and dplyr verbs checked column *presence* only | `mutate(res, padj = "oops")` still claimed to be a `gs_result`; `gs_filter(direction = "up")` then returned zero rows silently |
| compute | `[.gs_matrix` funnelled through the full constructor | A pathway filter matching nothing **aborted the script** instead of yielding an empty matrix |
| compute | Empty results dropped the method's optional columns | An empty contrast turned `gs_leading_edge()` into a usage error |
| DE | `xlim()`/`ylim()` in `de_pca()` | `xlim_abs` **deleted** out-of-range samples instead of zooming; only signal a warning at print time |
| DE | `x_breaks` in `orientation = "vertical"` | Retuned the *p* axis, not the fold-change axis it documents |
| DE | `[sig_logic]` kept `NA` positions | An all-`NA` `adj.P.Val` captioned the figure `p ≤ NaN` and skipped the documented "no genes pass" path |
| DE | `-log10(0)` | A `P.Value` of 0 made the y limit `Inf`, squashing every point at the bottom of a blank panel |
| IO | `gs_write()` only ever added files | A contrast dropped from the analysis reappeared in the figure, dated from the previous run |
| IO | No integrity check on a cached download | An interrupted transfer was reported as `[skip] (exists)` forever; `gatom_refs()` then died in `readRDS()` |
| IO | `/opt/gatom-refs` searched before the download destination | `gatom_refs(download = TRUE)` had no visible effect |
| plot | Per-figure fill limits | Same NES pale in one panel of a figure, saturated in the next |
| tests | Guards saw only the exact bug found | Deprecation targets that resolve to nothing, and S3 methods for absent generics, passed silently |
| tests | Every shim test picked a cutoff that keeps rows | The nothing-significant path was untested through the shims — which is why the §9 barplot bug shipped |

### Judgement calls worth keeping

- **`prune` is opt-in.** The stale-output fix could have been "`unlink()` the tree on every
  write". A function whose job is to write does not get to delete by default; the default
  writes a manifest and `gs_read()` warns, naming the files and the fix.
- **`gs_matrix` persistence was declined, not forgotten.** Giving it a table format is an API
  decision that belongs with the consumer migration, not a review side effect. The asymmetry
  is now stated in `gs_write()`'s docs, pointing at `saveRDS()`.
- **The two new structural guards were checked against planted faults** — a bogus
  `.Deprecated()` target and an S3 method for a generic nothing provides. A guard that has
  never been seen to fail is a guess about what it covers.
- **`download_gatom_references()` has frozen formals**, so there was no URL seam to inject; its
  two new tests mock `utils::download.file` and exercise the skip/rename logic, not the network.
- **`.gatom_search_dirs()` was extracted** so the resolution order is testable without a
  `/opt/gatom-refs` present. The old path had no test at all.

---

## 11. Release and dev-loop verification (2026-08-11)

**Origin was already current.** An earlier report in this plan said
`origin/feat/bulkirna-package` was 64 commits behind; that was a stale remote ref read
before a fetch. `git rev-list --count origin/feat/bulkirna-package..HEAD` is `0` at
`a484acd`.

**`v0.3.0` was force-moved from `163cc4f` to `a484acd`** and pushed to both `origin` and
`hub`. The old tag predated all eight high-severity and roughly twenty medium fixes and
had never been published, so nothing pinned it; the alternative considered was a fresh
`v0.3.1`, rejected because the restructure plan and the `scdock-r-dev` bump both name
`v0.3.0` by version. Its message now carries the corrected gate claim: `R CMD check` is
OK under `_R_CHECK_FORCE_SUGGESTS_=false`, and under `--as-cran` there is one Title-case
NOTE plus one environmental ERROR from `gatom`/`mwcsr` being absent from the dev image.

Gates re-run on `a484acd` immediately before tagging: **802 tests / 0 failures**, 5
expected gatom skips, **golden 20/20, exit 0**.

### Phase 2 — what is proven, and what is not

| Check | Result |
|---|---|
| `remotes::install_github("tony-zhelonkin/bulkiRNA@v0.3.0")` into a clean lib | installs from source, loads from both temporary and final location |
| Installed package exports | 61, as designed |
| Real analysis against the **installed** copy: `gs_test()` -> `gs_write()` -> `gs_read()` -> `gs_plot_dot()` | `leading_edge` round-trips identically; renderer returns a ggplot |
| `devtools::load_all()` over the dev checkout while the installed copy is on `.libPaths()` | the dev namespace shadows the installed one |
| A real source edit, then reload | takes effect (sentinel observed) |

The edit half was done in a `/tmp` copy of the checkout, not in the repo, so no dev-loop
proof left a working-tree change; the repo was confirmed clean afterwards.

**Phase 2 is discharged against a real consumer.** My earlier note in this section said the
consumer repos were absent from this host; they are at `/scratch/current/antonz/projects`
(`14839-DM-cGAS`, `STING-JR`, `DC-nexus`) -- my search had only covered `/data1` and
`/home`. `14839-DM-cGAS` carries `01_modules/RNAseq-toolkit` as a submodule pinned to the
old `scripts/` layout, on branch `dev`.

### Equivalence on real production data

`load_all()` over the dev checkout, then `gs_ranks()` + `gsdb_msigdb()` + `gs_test()` on
`14839-DM-cGAS`'s real DE table for contrast `STING_KO_vs_WT_t4` (18,558 symbols),
Hallmark, `seed = 123`, `n_perm_simple = 100000`, size bound 15/500 -- compared against
that project's *cached production* `gseaResult` from the old `run_gsea()`:

| Check | Result |
|---|---|
| Pathways | 50 old, 50 new, 50 matched |
| `setSize` vs `n_genes_tested` | identical |
| max abs dNES (`NES` vs `stat`) | `0` |
| max abs dp, max abs dpadj | `0`, `0` |
| Sign agreement, significant at 0.05 | 1.00, 23 old / 23 new |
| Top 10 by padj, in order | identical |

So the compute layer is bit-identical to what this project has already published, on its
own data -- not just on the fixture.

### What the shims do *not* carry (the Phase 4 work list)

Probed with the same real data, through `run_gsea()` and `normalize_gsea_results()`:

| Legacy pattern in `05_gsea_msigdb_run.R` | Behaviour now |
|---|---|
| `nrow(res@result)`, `res@result <- res@result[keep, ]` | **error** -- no `@` method for `gs_result`; the shim returns a tibble, not an S4 `gseaResult`. This is the blocking break, and it is loud. |
| `res@result$setSize` | error, same cause; the column is `n_genes_tested` |
| `if ("NES" %in% colnames(nt)) rename(nt, nes = NES)` | condition is `FALSE`, so the rename never fires and `nt$nes` is never created |
| `filter(sig, nes > 0)` | **error**, `object 'nes' not found` -- consequence of the line above |
| `filter(nt, direction == "Up")` | **silently 0 rows**; the values are `"up"`/`"down"`. The only genuinely silent break of the set. |
| `normalize_gsea_results(atlas_universe = ...)` -> `genes_full_set` | column absent. The shim warns that `atlas_universe` is ignored, so this is loud, but it is a real feature loss: this project's pathway-explorer input spec depends on `genes_full_set`, and the replacement is `gs_leading_edge()` against the original `gs_db`. |

The `ranks` and `gene_sets` attributes the `run_gsea()` shim promises are confirmed
present on the real result, so `gs_plot_running()`'s fallback path is intact.

Nothing in `14839-DM-cGAS` was modified: its working tree was already dirty from
unrelated devcontainer edits, the project was mounted read-only, and every output went to
`/tmp`.

## 12. The 677-vs-1050 Reactome gap: an upstream msigdbr bug (2026-08-11)

Phase 4's migration of `05_gsea_msigdb_run.R` reproduced the published Hallmark rows
exactly but appeared to test 1050 Reactome pathways where the published table has 677.
That gap is now explained, and it is **not** in this refactor.

### Cause

`msigdbr::msigdbr()` (confirmed 26.1.0) builds its ortholog table from *the collection
currently being queried*:

```r
species_genes <- babelgene::orthologs(genes = unique(mdb$db_ensembl_gene), species = species)
pkg_env[[paste0("orthologs", species_id)]] <- species_genes
```

The cache key names only the taxon, not the collection. The next call for a *different*
collection finds the key, reuses the previous collection's ortholog table, and
`inner_join`s against it -- dropping every gene absent from the first collection. It only
triggers when ortholog mapping is active (`db_species = "HS"` with a non-human `species`),
which is exactly this project's configuration.

Measured, `species = "Mus musculus"`, `db_species = "HS"`, deterministic and reproduced on
a pristine cache with live network (so it is not a cache artifact):

| call order | Reactome rows | sets | unique symbols |
|---|---|---|---|
| Reactome alone | 105806 | 1839 | 10762 |
| Hallmark, then Reactome | 44635 | 1817 | **3688** |

GO:BP is worse: 604644 rows / 15988 symbols alone, 248104 / 4313 after Hallmark.

### What this means for the published tables

`05_gsea_msigdb_run.R` loops the collections with Hallmark first, so **every collection
after Hallmark was tested against Hallmark's ~4.4k-gene space instead of the full ~10.8k**.
Consequences for `master_gsea_table.csv`: under-tested gene sets (677 of 1050 in-bounds
Reactome sets, all sizes agreeing exactly for the 677 that survived), truncated set
memberships, and -- because the BH family is the tested family -- wrong `padj` for every
non-Hallmark collection. Hallmark itself is correct; it was first.

Confirmation that the extras were never tested rather than merely filtered: across *all*
contrasts the published table contains exactly the same 677 Reactome ids, and 0 of the 387
extras appear under any contrast. Per-contrast size filtering could not produce that.

This is a data-correctness problem in already-published results, and it is a scientific
call for the user, not a refactor decision. Nothing in `14839-DM-cGAS` was regenerated.

### Fix in this package

`.msigdbr_drop_ortholog_cache()` (`R/gsdb-msigdb.R`, internal, `@noRd`) clears the
`^orthologs` keys from msigdbr's `pkg_env` before every `gsdb_msigdb()` query, so each call
is independent -- which is what `gsdb_msigdb()`'s docs already promise. It reaches into
another package's internals, so every step is guarded and failure is silent: if upstream
changes its caching the workaround degrades to a no-op, never to an error in a provider
call. Cost is one `babelgene` ortholog recompute per call, seconds against a warm cache.

Two false starts are recorded because both *passed* and neither tested anything:

1. First diagnostic compared `run_one_db()` against a bare `gs_test()` in **separate
   sessions** and got 1050 both times -- which is why the earlier "677 matched,
   max|dpadj| 0.0381" number could not be reproduced and was retracted. The variable was
   call order within a session, not the code path.
2. First regression test called `gsdb_msigdb()` for both arms. Earlier tests in the file
   had already poisoned the cache, so both arms were equally truncated and the test passed
   with the fix disabled. It now poisons through `msigdbr` directly, and takes its baseline
   from a cleared cache. Poisoning also only bites from a *cleared* cache -- msigdbr keeps
   whatever table is already there, so a small collection must go first.

The test skips (rather than fails) if upstream ever stops leaking, and a negative control
confirms it fails with the workaround disabled. Gates after the change: **806 pass / 0
fail** (was 802; 4 new assertions), golden **20/20**, `verify_golden.R` exit 0.

## 13. Hardening the guard, the documentation sweep, and the re-run (2026-08-11)

### The guard grew a second, independent layer

`.msigdbr_drop_ortholog_cache()` alone was brittle in the worst way: it reaches into
msigdbr's `pkg_env`, so an upstream rename turns it into a silent no-op -- and a silent
no-op restores exactly the silent truncation it exists to prevent. A fix whose failure mode
is invisible is not a fix.

`.msigdbr_assert_ortholog_coverage()` now checks the *result* and ignores how the mapping
was produced. It fetches the same collection with `species = "Homo sapiens"` (a path msigdbr
skips the ortholog branch for entirely, so it can neither be poisoned nor poison, and costs
one cached read) and compares unique gene counts. Below a 0.45 floor it errors.

The floor is measured. Legitimate coverage, `Mus musculus` / `db_species = "HS"`: H 1.002,
CP:KEGG_LEGACY 0.959, CP:REACTOME 0.940, GO:MF 0.904, GO:CC 0.897, GO:BP 0.878, C7 0.819,
C3 TFT:GTRD 0.597. Truncated: Reactome 0.32, GO:BP 0.24. 0.45 sits in the gap. A threshold
tight enough to grade ortholog *quality* would fire on legitimate collections; this one only
catches catastrophic loss.

Error rather than warning, because the bug's danger is that it yields plausible numbers -- a
warning in a pipeline log is indistinguishable from success.

Verified end-to-end: with the cache clear disabled, `gsdb_msigdb()` errors ("dropped 7,793
of 11,448 genes ... 32% retained") instead of returning 1817 truncated sets.

### Documentation sweep

Everything user-facing described the pre-package world. Rewritten and verified:

| File | Before | After |
|---|---|---|
| `README.md` | 534 lines of `source()`-era toolkit | 184 |
| `MIGRATION.md` | obsolete submodule branch strategy | the real script -> package guide |
| `docs/API_REFERENCE.md` | 947 lines duplicating roxygen | 92-line index |
| `docs/WORKFLOWS.md` | 396 lines of `source()` recipes | 108 |
| `docs/GSEA-workflow/` | 3739 lines, clusterProfiler-era, project-specific | archived to `docs/_internal/legacy-pre-package/` |

The archive is deliberate: those documents are the written record of the pipeline the golden
baselines were captured against, so they explain several design choices. Their README states
plainly that their examples describe the superseded API.

**Every documented code path was executed** against the installed package and real DE data:
21/21. That caught two errors that would have shipped -- `de_md_plot()` takes a limma `fit`
plus `coef` rather than a DE table, and `gsdb_info()` takes a database *name*. Documentation
that has not been run is a guess.

Released as **v0.3.1** (tagged, pushed to origin and hub) so the documented install line
carries the guard; `babelgene` joined `Suggests`, which the regression tests skip on.

One blemish left in history: a backtick in the commit message for `449e120` hit shell
substitution, so it reads "the silent lowercase  break". Fixing it needs a force-push to
rewrite pushed history, which is not worth it for a typo.

### The re-run

Authorized explicitly after the auto-mode classifier -- correctly -- refused a bare
"proceed" for it. Prepared: `bulkiRNA` 0.3.1 installed into a container library
(`.rlib-bulkirna`), both master tables backed up as
`master_{gsea_table,unified}.pre-msigdbr-fix.20260811-193204.csv`.

Baseline row counts in `master_gsea_table.csv` (12 contrasts x 6 collections):
GO_BP 29412, Reactome 8124, GO_MF 4548, GO_CC 3756, KEGG 1548, Hallmark 600 = 47,988 MSigDB
rows. The 1,488 rows from other scripts (MitoPathways 636, CoReSh 492, MitoXplorer 264,
TransportDB 96) survive, because `append_master_table()` replaces per database.

Expectations to check against: Hallmark holds at 50 per contrast (it ran first, so it was
never truncated), Reactome rises from 677 to ~1050, and GO_BP/GO_CC/GO_MF/KEGG all rise.
`padj` moves for every collection except Hallmark.

## 14. CoReSh: extraction plan and the package-versus-skill decision (2026-08-11)

Surveyed the three independent CoReSh implementations (DC-nexus `coresh-slice`,
`14839-DM-cGAS` scripts `07`/`08`, and the `coresh-signature-search` skill) against the
SciAgent-toolkit refactor. The plan of record lives at
[../2026-08-11-coresh-extraction/00_PLAN.md](../2026-08-11-coresh-extraction/00_PLAN.md).

Headlines:

- **The 269-line CoReSh engine is md5-identical in four checkouts** and is 0%
  project-specific. Across the 1,404 surveyed lines, ~81% is generic machinery.
- **The CoReSh `neg_log_padj` defect flagged in §13 is a column-projection bug**, not
  arithmetic. `08_coresh_derived_gsea.R:215` projects to a ten-column allowlist that omits
  `neg_log_padj`, `gene_ratio` and `leading_edge_size`; `bind_rows()` then NA-fills them.
  Confirmed in the live table: `padj=2.6175484e-33, neg_log_padj=NA`.
- **One column carries two saturation conventions** — `normalize_gsea.R` caps at 16, the
  migrated `05_gsea_msigdb_run.R` and the renderer `12_gsea_viz.R` reach ~307.65. Settle this
  before `06` and `08` migrate.
- **Recommendation: both.** `bulkiRNA` gains `gsdb_coresh()` plus the search, set-building
  and gene-identifier primitives; the skill keeps the judgement — query design, the top-20
  categorisation rubric, negative controls, the Synapse walkthrough. Boundary rule: a table
  in and a table out belongs to the package, a result in and a decision out belongs to the
  skill. `Suggests` grows by `qs2`, `BiocParallel`, `data.table`, `digest`.
- **`coresh_batch.R:44` calls `fgsea:::gesecaCpp`**, an unexported internal. Discharge it
  against public `fgsea::geseca()` the way `B2-fgsea-equivalence.R` discharged
  clusterProfiler; guard it only if that fails.
- **Phase 5 gains a concrete target.** `bulkiRNA` appears nowhere in SciAgent-toolkit; the
  existing `bulk-rnaseq-gsea` skill still instructs
  `source(.../GSEA/GSEA_processing/run_gsea.R)`, and `annotate-bulk-rnaseq-data` pins
  RNAseq-toolkit v0.2.0. The high-level bulkiRNA skill is a rewrite of those two, not a new
  skill. Edit the `scbio-docker/toolkits/` copy — it is ahead of the Meta-Aging checkout.

---

## 15. The four architecture decisions, settled (2026-08-12)

The seams surfaced by the CoReSh survey and the Phase 3 re-orientation were put to the
integrator as four questions. All four are now decided and have a durable home:
**[../../adr/](../../adr/README.md)** — a new directory, because `plans/` records what we
were about to do and these record what we chose and why.

| ADR | Decision | Implementation |
|---|---|---|
| **001** | The **package version** is the unit of reproducibility. An image tag only promises which tool versions it contains; analysis repos carry their own lock. | provenance must emit the version closure; `Version` may not move without a tag |
| **002** | **`bulkiRNA` owns the master-table schema and validator**, not the file. `neg_log_padj = -log10(pmax(padj, .Machine$double.xmin))`. | T1 — `gs_to_master()`, `gs_validate_master()`, versioned schema in `inst/extdata/` |
| **003** | **`Suggests` stays**; no packages move. Add a bulk preflight and `Remotes`. | T2 — `bulkirna_check_deps()` |
| **004** | **Two reference-data tiers only**: bundled (no fetchable source) or refcache (has one). One resolver; provenance records the resolved snapshot tag. | T3 — `.ref_path()` |

### What these change in the existing plan

- **C0 is revised.** `gs_to_master()` goes **into the package**, not into
  `02_analysis/helpers/de_gsea_helpers.R`. One tested implementation retires all three
  drifted ones, and the saturation convention stops being a rule that three files each
  remember differently.
- **`Sys.getenv("CORESH_CHUNKS")` in the CoReSh plan §4 is superseded** by `.ref_path()`.
  T3 must land before `gsdb_coresh()`.
- **Phase 3 is no longer the only thing that can move.** T1–T3 need no rebuild, and T2 is
  what makes the eventual rebuild *verifiable* — `bulkirna_check_deps("all")` as a
  post-build smoke test, instead of grepping for the five missing packages by hand.
- **`download_gatom_references()` is on notice.** Its replacement is `sources/gatom.sh` in
  the refcache repo, which is a separate repo and needs its own go-ahead.

### Delegation

T1–T3 are delegated to **codex (`gpt-5.6-sol`)** via `codex exec`, run sequentially — they
all touch `NAMESPACE`, and parallel runs would collide on it. **Codex never runs R**: there
is no R on this host, and the container gates are integrator-owned. Each handback is gated
by `devtools::document()`, `devtools::test()` and `Rscript tests/golden/verify_golden.R`
before it is committed.

---

## 16. v0.4.0, and Phase 3 staged (2026-08-12)

### The release

`v0.4.0`, tagged and pushed to `origin` and `hub`. Three new exports, so minor rather than
patch. Verified to install from GitHub into a scratch library: 64 exports, and — the check
`devtools` cannot make — the schema CSV resolves through `system.file()` from an *installed*
package, `gs_to_master()` returns its 14 columns, and the validator reports zero problems.

`pheatmap` left `Suggests`; `babelgene` stayed with a comment. The correction worth recording:
**`babelgene` was never dead.** `msigdbr` calls it internally for ortholog mapping and
`test-gsdb-msigdb.R:53,115` skip on it, so a search for direct `babelgene::` calls finds
nothing and concludes wrongly. Those are the tests guarding the §12 cache bug. The comment
exists so the next reader does not repeat the finding.

### Phase 3 — staged, not built

Branch `feat/bulkirna-v0.4.0` in `scbio-docker`, commit `0672cb1`. **Not pushed, not built.**

| Change | Detail |
|---|---|
| `docker/base/R/install_core.R:148` | `"tony-zhelonkin/bulkiRNA@v0.4.0"` in `github_packages` |
| `install_gh_pkg()` | now splits `owner/repo@ref`; **no `GH_PKG_NAME` entry needed** |
| preflight | writes `/opt/settings/bulkirna_optional_deps.csv` after the GitHub installs |
| `VERSION` | `v0.5.12` → `v0.5.13`, plus a `docs/changelog.md` entry |

Three findings from doing it:

- **The plan's prescription was incomplete.** A `GH_PKG_NAME` entry fixes the derived package
  name but not `install_gh_pkg()`'s `install_git` fallback, which builds its URL with
  `sprintf("https://github.com/%s", slug)` — for a ref-suffixed slug that is not a URL, so the
  fallback would fail for a reason unrelated to the network. Teaching the function about refs
  fixes both paths and makes the `GH_PKG_NAME` entry unnecessary, since with the ref stripped
  the repo and package names match.
- **The five-missing-`Suggests` sub-decision is moot.** `install_core.R:117,122,125` already
  install `homologene`, `babelgene`, `FactoMineR`, `factoextra`, `org.Hs.eg.db` and
  `org.Mm.eg.db`. Only `gatom` and `mwcsr` are genuinely absent, and they are optional.
  `pheatmap` is installed image-wide at `:74`, so dropping it from `Suggests` costs consumers
  nothing.
- **The preflight must not write to `install_failures.csv`.** That file means "requested and
  failed", `AGENTS.md:82` tells the reader to check it after every build, and its "no failures"
  line is a real signal — earned after seven absent packages hid for two releases. Optional
  packages that are *expected* to be absent would put permanent rows in it and retire the
  signal. Separate artifact; a `bulkiRNA` that will not load is still a genuine failure and
  still recorded there.

### What is left in Phase 3

The rebuild, which runs long on the owner's shared infrastructure, and pushing the branch.
Both are the owner's to trigger. Open question surfaced by the preflight: whether to add
`gatom` and `mwcsr` to the image so the network-module features work out of the box, or leave
them to on-demand install.

---

## 17. Phase 3 complete (2026-08-12)

`scdock-r-dev:v0.5.13` built in 150m 2s. Zero `ERROR` lines, zero `FAILED:` lines,
`/data1` steady at 471G free.

The gatom question from §16 was answered by installing both: `mwcsr` first, since `gatom`
calls its solvers. `mwcsr` compiled its native code cleanly at `-std=gnu++17`, which works
because of the `Makevars.site` fix that shipped in v0.5.12.

### Verified in the built image

| Check | Result |
|---|---|
| `/opt/settings/install_failures.csv` | header only — the signal `AGENTS.md:82` relies on is intact |
| `/opt/settings/bulkirna_optional_deps.csv` | 16 rows, all `installed = TRUE` |
| `library(bulkiRNA)` | 0.4.0, 64 exports |
| `system.file("extdata", "master-schema-v1.csv")` | resolves from the installed package, 15 rows |
| `bulkirna_check_deps("all")` | 16/16 present |
| `gatom` / `mwcsr` | 1.8.4 / 0.1.12 |

The preflight replaced a manual grep on its first real run, and it caught nothing because
there was nothing to catch. That is the outcome to want from a check that costs seconds.

### What Phase 3 unblocks

Phase 4's remaining consumer migrations, `10_gatom_modules.R` among them, which the image can
now actually run. And C5/C6, which needed an *installed* `bulkiRNA` because nested submodules
are forbidden — though those wait on the SciAgent-toolkit refactor for separate reasons.

---

## 18. Phase 4 — two scripts now, in `14839-DM-cGAS` (2026-08-12)

### `08_coresh_derived_gsea.R` — migrated, and C1 discharged

Full account in [../2026-08-11-coresh-extraction/00_PLAN.md §9](../2026-08-11-coresh-extraction/00_PLAN.md).
The short version: both column allowlists deleted, 708 CoReSh rows all carrying finite
`neg_log_padj` against 492 rows carrying none, and a stale `load_or_compute()` cache found
underneath — cached results four days older than the GMT they claimed to describe, which is
why the table listed nine sets that no longer existed and omitted 27 that did.

### `05_gsea_msigdb_run.R` — off its local `gs_to_master()`

The script carried its own copy, written before the package had one, with a note to lift it
when 06 and 08 migrated. It now calls `bulkiRNA::gs_to_master()`.

Proven equivalent before the swap, not after: both implementations were run over all twelve
cached results and compared column by column. **72,408 rows, identical on every value of every
column.** Re-running the script then left both master tables byte-identical in row count,
per-database breakdown and key uniqueness.

Two things the package function leaves to its caller, and the script now does explicitly:
`pathway_name` formatting, and `database`/`contrast`, which come from the `gs_result`'s own
columns rather than from arguments.

Also fixed a message that read `nrow(gt)` and called it "MSigDB rows". The file also holds
CoReSh, MitoPathways, MitoXplorer and TransportDB rows, so it overstated this script's output
by 1,716. It now reports 72,408 — the same number the equivalence check covered, which is a
free cross-check that the two agree on what the script produced.

### Phase 4 remaining

`06_gsea_custom_run.R`, `12_gsea_viz.R`, `10_gatom_modules.R` in this project — the last of
which `v0.5.13` can run now that `gatom` and `mwcsr` are installed. Then STING-JR (42 refs),
then DC-nexus.

**One thing to decide first, given §9's finding:** whether to survey the other `load_or_compute()`
stages for the same staleness. The cache key is the filename, so any stage whose input changed
without its name changing carries the same exposure. `05`'s caches are newer than their inputs,
so it is clean; nothing else has been checked.

`14839-DM-cGAS` stays uncommitted, now with `master_*.pre-coresh-c1.20260812-225911.csv`
backups beside the tables. Committing there remains the owner's call.

---

## 19. Phase 4 complete for `14839-DM-cGAS` (2026-08-13)

The remaining three scripts migrated, delegated to `codex gpt-5.6-sol` as three parallel
fan-outs, each in an isolated scratch copy with the package source as read-only reference. I
gated every one by running it in `scdock-r-dev:v0.5.13`. Running codex with write access and no
approval gate against the live project was refused, correctly — the scratch-copy-then-apply shape
is the one to keep.

### `06_gsea_custom_run.R`

`gsdb_register()` around the existing phase-04 checkpoints rather than `gsdb_load()`, so the
bundled snapshot cannot silently replace the ids that key existing master rows. Both allowlists
gone. **1,008 custom rows, exactly the previous count**, all with finite `neg_log_padj`,
`gene_ratio` and `leading_edge_size`; MSigDB and CoReSh rows untouched.

**One defect I caught in review.** The draft formatted `pathway_id` into `pathway_name`, copying
05. That is right for MSigDB and wrong here: the custom collections ship a real `TERM2NAME`
description, which `gs_test()` already carries through from the `gs_db`. `TRANSPORTDB_AAAP` would
have gone from "The Amino Acid/Auxin Permease (AAAP) Family" to "Aaap" — 84 sets relabelled with
their own prefixes. Formatting `pathway_name` instead fixes it.

**A divergence that remains, and is worth a decision.** 33 of the 84 custom names still changed,
and every one is pure capitalisation — case-insensitively identical, nothing substantive. The
package's `format_pathway_name()` is deliberately smarter than the toolkit's ("Metabolism of
Lipids", "Fe-S Cluster" rather than "Metabolism Of Lipids", "Fe-s Cluster"), with one real bug:
a sentence-initial small word stays lowercase, so 9 names now begin "the" instead of "The".
Narrow fix, wide blast radius — it would move display names in every project and may move golden
baselines, so it is the owner's call rather than something to slip in here.

### `12_gsea_viz.R`

Off all seven toolkit `source()` calls and `enrichplot`. A documented local
`master_to_gs_result()` adapter feeds `gs_plot_dot()` and `gs_plot_bar()` from the master table.
`get_db_plot_params()` and `smart_wrap()` were sourced but never called — verified against the
pre-migration file before deleting the dependency.

**The first draft dropped the running-sum figures**, on the reasoning that `gs_plot_running()`
needs ranks and set membership that master rows cannot supply. Both halves true, conclusion
wrong: `.grs_ranks()` and `.grs_sets()` take explicit arguments, so ranks rebuild from the DE
table and the `gs_db` rebuilds per database. That omission would have destroyed **120 figures**,
and `CONFIG$figures$running_sum_ylim` was added on 2 July precisely to make those curves
comparable. Sent back as a correction. After it: **484 PDFs, byte-for-byte the same inventory as
before — 0 missing, 0 new** — with the 120 running-sum figures intact.

### `10_gatom_modules.R`, and the finding that matters most

The first attempt **stopped without editing** and reported two package gaps. Both were real; I
verified each against the source, closed them in `8a442eb`, and relaunched. It then found a
third and again stopped rather than working around it. That is the behaviour the briefs asked
for and it held three times out of three.

The gaps: `gatom_refs()` accepted only `"kegg"` while the config asks for `["kegg", "combined"]`;
`gatom_module()` computed the solver's objective weight and discarded it; and it cannot pass
`gene2reaction.extra`, which the combined network needs — so combined stays a direct
`gatom::` call, and threading that argument through is the next package task.

**The first gate was hollow and I nearly accepted it.** The script exited 0 having printed
`cache hit: gatom_results.rds`, so the migrated code never ran. Moving the cache aside and
recomputing exposed the real result: **every `kegg` cell changed, and every `combined` cell was
bit-identical.** Since combined is the unmigrated path, the cause is `gatom_module()` itself.

It is the seeding. `gatom_module()` seeds before `scoreGraph()`; the script only seeded before
solving. Two runs of the script's own sequence, same input, same session:

```
scoreGraph unseeded (the script): 165.8478, then 169.3227   -> not reproducible
both seeded (the package)       : 169.3227, then 169.3227   -> reproducible
```

**`gatom::scoreGraph()` is stochastic, and this project's GATOM modules were never reproducible.**
The cached 165.85 was one draw. The migration does not introduce variance, it removes it, and the
new numbers are the ones that repeat. Pre-migration caches preserved as
`gatom_*.pre-migration-20260812.rds`. 24/24 modules non-empty, 24 rows in `master_unified.csv`.

That makes four instances in this refactor of a number nobody could have known was wrong by
reading the code, and the second where the evidence only appeared because a cache was forced to
recompute.

### Staleness survey, discharged

§18 left open whether other `load_or_compute()` stages carried CoReSh's staleness. Checked:
they do not. `analysis_config.yaml` is newer than every GSEA cache, but its `thresholds.gsea_*`
values and its whole `databases:` block are **unchanged** since those caches were built, so the
newer mtime is drift in the TE and integration sections. The TE GMTs are newer than the TE
caches, but they are that stage's own **outputs**, not its inputs.

So CoReSh was singular for a specific reason: its GMT is a genuine cross-stage input, written by
phase 07 and read by phase 08, and phase 07 was re-run afterwards. The pattern to watch is a
cache whose input is another stage's generated artifact — not one whose input is the config.

### Phase 4 state

`14839-DM-cGAS` is done: all five GSEA-path scripts on the package. Next are STING-JR (42 refs)
and DC-nexus. The project tree remains uncommitted, with backups beside every table and object
this pass rewrote. Committing there is the owner's call.

## 20. The two open nicks, closed (2026-08-13)

Both delegated to `codex gpt-5.6-sol` in scratch clones, both gated here. `483754c`.

### `format_pathway_name()` capitalises sentence-initial function words

Only the nine function words are affected. A Greek letter, a chemical prefix, a Roman numeral
and a deliberately lowercase symbol like `mTOR` are all correct at the start of a label and stay
as they are; mid-label the function words still lowercase, so the documented
`TNF-alpha Signaling via NF-kappaB` is untouched. This is one of the frozen 24, and the freeze
covers signatures, so a behaviour change is permitted — but it changes display names everywhere,
which is why it waited for an explicit decision.

**I broke it once while tidying.** The word list was duplicated between the exceptions map and
the new rule, so I gave it one source in `.gs_function_words()` — and spliced a character vector
into a `list()`, which buries all nine words in a single unnamed element. Every lookup for them
missed, and mid-label `via` silently became `Via`. `as.list()` fixes it. It was caught by
diffing against the previous implementation on real inputs, not by reading the change, which is
worth remembering: the diff looked right.

The golden `format_pathway_name` baseline did not move, as codex predicted from its four inputs.

### `gatom_module(gene2reaction_extra = )`

Optional, `NULL` by default, last in the formals, so no existing call changes — and codex added a
test asserting the first eight formals keep their names and positions, which is a better freeze
guard than the prose rule. It checked upstream's own default in `gatom`'s source rather than
assuming. Verified against the real reference tree: the mapping takes the combined graph from
7,608 to 8,179 edges, the wrapper reproduces the direct call's graph exactly, and a bad value is
rejected with a message naming the fix. The combined network no longer needs to bypass the
wrapper, though `10_gatom_modules.R` has not been rewritten onto it yet.

Declined with reasoning: having `gatom_refs()` resolve and return the parsed table would change
the shape of an object that already has consumers.

### The consumer, settled

`06` re-run and `12` re-run against the fixed package. **Zero custom labels now begin lowercase.**
33 of 84 still differ from the pre-migration table, all capitalisation-only, and the package's
form is the better one throughout — `TCA Cycle` for `Tca Cycle`, `ROS and Glutathione Metabolism`
for `Ros And Glutathione Metabolism`, `Branched-Chain` for `Branched-chain`. 1,008 custom rows,
all with finite `neg_log_padj`. Figures regenerated: **484 PDFs, 120 running-sum, 0 missing,
0 new.**

Package state: 942 tests passing, golden 20/20, 64 exports. Every gate green.

---

## 21. Phases 8-11: the analysis-API roadmap (2026-08-13)

This plan covers Phases 0-7 — getting the GSEA layer out of the scripts and into a versioned
package. The owner has now named the destination beyond that: a thin uniform API that carries
his methodological preferences as defaults, so an agent implementing an analysis inherits them
instead of reading four previous projects and guessing.

Four documents own that work. **The sequencing was revised the same day**: stabilize the
current surface first, complete one CoReSh/GESECA run at full fidelity, and park the rest.

- [../2026-08-13-analysis-api-roadmap/00_ROADMAP.md](../2026-08-13-analysis-api-roadmap/00_ROADMAP.md)
  — **Phase 8 stabilize the surface** (stability contract, name and argument audits, deprecation
  clock, vignettes) and **Phase 9 one standard CoReSh and GESECA run**. Includes what each ADR
  now has to carry and the tradeoffs of the alternatives.
- [../2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md](../2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md)
  — what CoReSh's score is, traced to `fgsea:::calcGesecaScores`, why the PCA compression is
  lossless for it, and exactly what running the authors' intended analysis still needs.
- [../2026-08-13-analysis-api-roadmap/03_DEFERRED.md](../2026-08-13-analysis-api-roadmap/03_DEFERRED.md)
  — TF activity, PROGENy and WGCNA, parked by decision with the survey findings kept intact.
- [../2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md](../2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md)
  — the canonical reference implementation for each capability, with paths, conventions,
  parameter choices and their reasoning, and the defects not to carry forward. Read it before
  designing an API.

**Position after this session:** Phase 3 done; Phase 4 done for `14839-DM-cGAS` with STING-JR
and DC-nexus remaining; CoReSh C0-C2 done and C3 half done; Phases 5-6 blocked on the
SciAgent-toolkit refactor; Phase 8 (stabilize) and Phase 9 (one CoReSh/GESECA run) planned and
unstarted, with the activity and network layers parked. `v0.5.0` tagged. Package state:
**1060 tests passing, golden 20/20, 72 exports.**

`R CMD check` on `v0.5.0` is **0 errors, 0 warnings, 0 notes**, so Phase 8 is a contract job
rather than a repair job: the surface is technically sound and rhetorically unfinished. 72
exports, 20 of them deprecated shims with no removal date, and nothing anywhere saying which
names are safe to build on.

**Start at S1, the stability contract, then G1.** G1 removes a guard, redirects one argument to
`fgsea::geseca()`, and closes a decision that should never have been opened — see §12 of the
CoReSh plan for the measurement that corrected it.
