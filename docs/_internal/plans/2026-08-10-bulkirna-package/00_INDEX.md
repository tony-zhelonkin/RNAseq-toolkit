# RNAseq-toolkit → `bulkiRNA` package — working plan and state

**Branch:** `feat/bulkirna-package` (off `dev`)
**Fork point / revert target:** `752481f` (`dev` tip at fork)
**Baseline at fork:** 9 legacy suites in `tests/`; 20 golden cases captured, 0 errors
**Current:** Step 0 complete. Nothing in `scripts/` modified yet.

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

### Steps A, B1–B5, C — ⬜ not started

`scripts/` is untouched. No `DESCRIPTION`, no `R/`, no `NAMESPACE`.

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

## 3. Next action — Step A (serial, blocks all five B agents)

Create the skeleton and the shared contracts. Per `08_refactor-execution-plan.md` §3:

- `DESCRIPTION` — `Package: bulkiRNA`, `Version: 0.3.0`, the Imports/Suggests ledger from
  `07_api-design.md` §8 verbatim.
- Layout: `R/`, `inst/extdata/` (the 128 KB `data/references/*/processed/` tree +
  `METADATA.yaml`), `tests/testthat/`, `.Rbuildignore` excluding `data/references/*/raw/`.
- `R/gs-result.R` — constructor, validator, S3 `print`/`summary`/`as_tibble`/`rbind`.
  Core columns exactly as `07` §4.
- `R/gs-matrix.R` — constructor + methods.
- `R/utils.R` — the **single** surviving `ensure_dir()` (directory semantics, per `02` §3.2),
  `%||%`, nothing else. The one place shared helpers live.
- `CONVENTIONS.md` (repo-internal) — roxygen style, `@keywords internal` for non-exports,
  `@import ggplot2` in exactly one file, the `stat_type` vocabulary.
- NAMESPACE is **roxygen-generated, never hand-edited** — this is what stops five agents
  conflicting on one file.

Then dispatch B1–B5 in worktrees. Ownership rule: `DESCRIPTION`, `NAMESPACE`,
`R/gs-result.R`, `R/gs-matrix.R`, `R/utils.R` belong to Step A and the integrator; B agents
**report** needed changes rather than making them.

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

## 6. Revert

```bash
git checkout dev && git branch -D feat/bulkirna-package   # nothing else was touched
```

`scripts/` is unmodified as of Step 0, so abandoning this branch costs only the fixture and
golden harness — both of which are worth keeping regardless.
