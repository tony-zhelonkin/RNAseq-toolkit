# RNAseq-toolkit → `bulkiRNA` package — working plan and state

**Branch:** `feat/bulkirna-package` (off `dev`)
**Fork point / revert target:** `752481f` (`dev` tip at fork)
**Baseline at fork:** 9 legacy suites in `tests/`; 20 golden cases captured, 0 errors
**Current:** Steps 0 and A complete. Nothing in `scripts/` modified yet. Next: dispatch B1–B5.

Resume by reading §1 (state), then §3 (next action). Everything needed to continue is here.

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

### Steps B1–B5, C — ⬜ not started

`scripts/` is untouched.

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

---

## 3. Next action — dispatch B1–B5 in worktrees

Briefs are `08_refactor-execution-plan.md` §4; the rules every agent follows are §5, now
also written down in-repo as **`CONVENTIONS.md`** (read that first — it is the contract).

Ownership rule: `DESCRIPTION`, `NAMESPACE`, `R/utils.R`, `R/gs-result.R`, `R/gs-matrix.R`,
`R/bulkiRNA-package.R` belong to Step A and the integrator; B agents **report** needed
changes rather than making them. `tests/golden/` and `tests/fixtures/` are off limits.

Two Step-A facts the briefs did not anticipate:

- The `@import ggplot2` lives in `R/bulkiRNA-package.R`, so B3/B4/B5 need **no** `ggplot2::`
  prefixes and must not add a second `@import ggplot2`.
- `@export` on an S3 method whose generic is owned by another package emits
  `export()`, not `S3method()`. Use `@exportS3Method pkg::generic` (e.g.
  `@exportS3Method tibble::as_tibble`). Verify `NAMESPACE` after `document()`.

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

## 5a. Deferred: invert the GATOM ownership (post-Step C)

**Decided in principle, deliberately not in this pass.**

Today: `SciAgent-toolkit/skills/gatom-metabolomic-predictions/SKILL.md` is **prose only** —
no code — re-teaching GATOM's raw API, while this repo vendors only the reference-file
downloader. The skill carries four traps in prose that a function signature should enforce:
raw p-values not `padj`; `baseMean` on **linear** scale; genes on graph **edges**, not
vertices; `met.db` required even when `met.de = NULL`.

Target: a `gatom_*` layer here — `gatom_refs()` (absorbing `download_gatom_references()` and
its cache), `gatom_graph()` (validating the `gene.de` contract instead of documenting it),
`gatom_module()`, `gatom_plot()` — with `SKILL.md` reduced to a pointer at bulkiRNA.

Why not now: `07_api-design.md` scopes this package to fgsea + GSVA + ORA, so this is a **new
feature module, not a refactor**; there is no golden baseline for GATOM behaviour, making it
the one part of the pass without a safety net; and `gatom` is a new heavy Bioconductor
`Suggests`, i.e. a ledger change mid-flight.

Meanwhile `download_gatom_references()` stays **exported with its exact current name and
formals** — it is on the frozen 24-export list, and `02_api-inventory.md` §5 says to export
it regardless of having no in-repo caller. B1 deleted it on a "no consumer" reading; that
test does not override the freeze. Restored.

---

## 6. Revert

```bash
git checkout dev && git branch -D feat/bulkirna-package   # nothing else was touched
```

`scripts/` is unmodified as of Step 0, so abandoning this branch costs only the fixture and
golden harness — both of which are worth keeping regardless.

---

## 7. Kickstart prompt for a fresh session

```
Resume the bulkiRNA packaging work. Read, in order:

  1. docs/_internal/plans/2026-08-10-bulkirna-package/00_INDEX.md   (this repo — state + next action)
  2. /data1/users/antonz/pipeline/sciagent-rna/docs/07_api-design.md (the API contract)
  3. /data1/users/antonz/pipeline/sciagent-rna/docs/08_refactor-execution-plan.md §3 (Step A spec)

You are on branch `feat/bulkirna-package`. Step 0 is done: 20 golden cases in
tests/golden/, fixture in tests/fixtures/, `scripts/` untouched.

Task: execute **Step A** — the skeleton and shared contracts. DESCRIPTION,
inst/extdata/, R/gs-result.R, R/gs-matrix.R, R/utils.R, CONVENTIONS.md,
tests/testthat/ scaffold, .Rbuildignore. Do NOT start B1–B5 and do NOT
modify anything under `scripts/` yet.

Hard constraints:
- No R on this host. Everything runs in a throwaway container:
    docker run --rm --user "$(id -u):$(id -g)" -e HOME=/cache \
      -v /data1/users/antonz/pipeline/.msigdb-cache:/cache \
      -v "$PWD":/pkg -w /pkg scdock-r-dev:v0.5.11 <cmd>
  Both --user and HOME are mandatory (saveRDS permissions; msigdbr's runtime cache).
- NAMESPACE is roxygen-generated, never hand-edited.
- The 24-export list in 02_api-inventory.md §5 is FROZEN: no renames, no argument
  changes, until Step 1b.
- `Rscript tests/golden/verify_golden.R` must exit 0 before you hand back. Confirm it
  still passes at the start too, so a later failure is attributable.
- Commit on the branch; push to `hub` and `origin` (origin needs an ssh-agent with a
  passphrase key — ask rather than working around it).

Report: what you created, the R/utils.R and DESCRIPTION contents, golden-verify result,
and anything in the plan you found wrong.
```

To resume the **SciAgent** track instead, point a session at
`scbio-docker/toolkits/SciAgent-toolkit/docs/_internal/plans/2026-08-07-demolish-role-layer/00_INDEX.md`
(branch `refactor/demolish-role-layer`) and read §1 then §3.
