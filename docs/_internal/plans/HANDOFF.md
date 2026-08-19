# bulkiRNA — handoff and plan of record

**Updated:** 2026-08-18 · **Branch:** `feat/bulkirna-package` · **Version:** `v0.5.0` tagged
**State:** review fixes in flight; gates require a fresh run · 78 exports

This is the entry point. Everything else is reachable from here.

| Document | What it owns |
|---|---|
| [2026-08-10-bulkirna-package/00_INDEX.md](2026-08-10-bulkirna-package/00_INDEX.md) | Phases 0–7 and their execution record |
| [2026-08-13-analysis-api-roadmap/00_ROADMAP.md](2026-08-13-analysis-api-roadmap/00_ROADMAP.md) | Phases 8–9, and §7–§9 record what was executed |
| [2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md](2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md) | **The canonical previous implementation of each analysis, with paths.** Read before designing anything |
| [2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md](2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md) | What CoReSh's score is, traced to its implementation |
| [2026-08-13-analysis-api-roadmap/03_DEFERRED.md](2026-08-13-analysis-api-roadmap/03_DEFERRED.md) | TF activity, PROGENy, WGCNA — parked, with findings intact |
| [2026-08-11-coresh-extraction/00_PLAN.md](2026-08-11-coresh-extraction/00_PLAN.md) | CoReSh extraction, C0–C6 |
| [../adr/](../adr/) | The four architectural decisions |
| [../../../AGENTS.md](../../../AGENTS.md) | The rules an agent must follow in this repo |

---

## 0. Why this project exists, in one paragraph

`RNAseq-toolkit` was a folder of scripts consumed as a git submodule and pulled in with
`source()`. 24 copies, each pinned to a different commit, so a fix in one reached nobody.
Nothing was tested. Nothing carried a version, so no result could name the code that produced
it. **The goal is not tidiness. It is that a number in a figure can be traced to a version of
the code that produced it, and that a defect fixed once is fixed everywhere.**

The second goal, stated by the owner and now the reason Phases 8–11 exist:

> *a thin clean API that just does the job the way I like it, without having to re-state to
> agents implementing analysis for me to go search for my previous experience reference repos —
> a zero-token architecture for my personal preferences on how bulk and pseudo-bulk should be
> done.*

That reframes the package as **the executable record of a set of methodological preferences**.
Two consequences set the acceptance bar:

1. **A function needing a paragraph of instructions to use correctly has failed.** The defaults
   *are* the preferences, documented once.
2. **A choice made the same way three times is a default, not a decision.** Three projects
   agreeing on `top = 500` for PROGENy is no longer a judgement call.

---

## 1. The pain points

The first four are original; the rest were found during the work.

**1. Silent wrongness that a green run cannot distinguish from a null result.** Ten confirmed
instances, every one measured rather than inferred:

| # | Defect | How it hid |
|---|---|---|
| 1 | `msigdbr` keys its ortholog cache on species alone | Second collection in a session truncated to the intersection; filed as issue 62 |
| 2 | A ten-column allowlist NA-filled four columns | 492 rows with `padj = 2.6e-33` beside `neg_log_padj = NA` |
| 3 | `annotate_genes(use_biomart = TRUE)` returned all-NA | Silently, after being explicitly asked for biomaRt |
| 4 | `safe_install()` logged two rows per failure | One falsely "returned without error" |
| 5 | `load_or_compute()` served results older than their input | Table listed 9 sets that no longer existed, omitted 27 that did |
| 6 | `gatom::scoreGraph()` is stochastic and was never seeded | Same input scored 165.85 then 169.32 in one session |
| 7 | A cached `NULL` made a failure outlive its own fix | Next run reads the file, skips the contrast, reports the same clean nothing |
| 8 | A green gate that never ran the code | Exit 0 printing `cache hit`; forcing a recompute showed every cell had moved |
| 9 | **`bplapply()` switches `RNGkind()`, so seeding without pinning changes answers** | `gs_test()` returned different p-values inside a parallel worker, no warning |
| 10 | **`NA` Entrez ids in 0.7% of real CoReSh datasets** | Passed the entire test suite; found by a compendium sweep |

Numbers 5–10 share one shape: **an empty result and a broken run were indistinguishable from
inside the code, so the broken run took the empty result's path.** In one case that purged 492
rows from two master tables and exited 0.

**2. Machinery retyped per project, then drifted.** ~81% of the surveyed CoReSh lines were
generic; the engine was md5-identical in four checkouts, re-vendored by a `cp` in a comment.

**3. The same drift reproduced *inside* the package.** Seeding was implemented three times,
three ways. TF activity is next: two estimators and three network provenances across three
projects.

**4. Reference data fetched live, so a result depends on the day.** `msigdbr`'s cache, the GATOM
downloads, and a documented OmniPath breakage.

**5. Provenance that can be silently meaningless.** `.ref_path()` recorded
`basename(Sys.readlink(current))`; a `current` pointing outside its source, or being a real
directory, recorded a snapshot identifying nothing while the call succeeded.

**6. The agent-facing documentation described the dead predecessor.** `AGENTS.md` — the one file
every assistant reads, via `CLAUDE.md` and `GEMINI.md` — was 225 lines about a script library,
including "there is no `library(RNAseqToolkit)`" and submodule-pinning instructions. Rewritten
2026-08-14.

---

## 2. The ADRs, with premise and rejected alternatives

**ADR-001 — the package version is the unit of reproducibility.** *Premise:* four version
authorities floated free and no result could name its code. *Decisive argument:* under
image-as-unit the `msigdbr` bug would have been **invisible** — tag unchanged, numbers changed.
*Rejected:* image-as-unit, which hides dependency drift and couples a scientific claim to an
artifact rebuilt for unrelated reasons; a lock file inside `bulkiRNA`, which a library cannot
impose on its consumers. *What it grew:* the stability contract is its missing half. A version
says *that* something changed; a lifecycle tier says *whether you were entitled to rely on it*.

**ADR-002 — the package owns the master-table schema and validator, not the file.**
`neg_log_padj = -log10(pmax(padj, .Machine$double.xmin))`. *Premise:* whoever computes a column
owns its definition, or the definition lives somewhere untestable. *Vindicated:* `37.55` in the
live table, above the retired cap of 16, on real values. *Rejected:* toolkit ownership — a skill
cannot be `library()`d, so the schema reverts to prose, which is how it drifted three ways;
per-project schema files, which drift by default. *Held under pressure:* when `gs_test()` gained
`log2err`, `gs_to_master()` still returned exactly 14 columns, verified rather than assumed.

**ADR-003 — `Suggests` stays.** *Premise:* it is the machine-readable list `R CMD check`,
`install_deps` and the image build read. *Decisive argument:* 14 light `Imports` are why
`install_github` works on bare R in seconds. *Rejected:* collapsing into `Imports`, making every
optional feature a hard install blocker; core-plus-feature packages, cleanest isolation but a
second release train for one maintainer. *Amended twice:* `Remotes:` was declared unnecessary,
then necessary once `coresh` looked GitHub-only, then unnecessary again once `coresh` turned out
to export nothing. *Watch item:* WGCNA's dependency tree is the first that makes the split
arguable; a second like it should reopen the decision.

**ADR-004 — two reference-data tiers, bundled or refcache.** *Premise:* three mechanisms had
grown and a fourth environment variable was about to appear. Bundling MitoCarta stands: a source
with no API cannot be re-fetched reproducibly, so vendoring **is** the archive, and 168 KB makes
the size objection moot. *The insight:* nobody needs to own pinning, someone needs to own
**recording**. Freshness stays the default; the drift stops being silent. *Earning its keep
twice:* CoReSh chunks resolve through it, and the OmniPath networks are exactly its case —
fetchable in principle, unfetchable in practice at current package versions.

**ADR-005 was proposed and withdrawn.** It was to settle ULM versus MLM for TF activity. With
that layer parked there is no decision to record, and an ADR for unscheduled work rots. It
survives as an open question in `03_DEFERRED.md`.

---

## 3. Where we are

| Phase | Status |
|---|---|
| 0 Inventory & freeze | ✅ |
| 1 Skeleton, internal `source()` removed | ✅ |
| 2 Prove the dev loop | ✅ against real data |
| 3 Install into the image | ✅ `scdock-r-dev:v0.5.13`, 16/16 optional deps |
| 4 Migrate heavy consumers | 🟡 `14839-DM-cGAS` done (5/5 scripts); **STING-JR and DC-nexus remain** |
| 4c CoReSh extraction | 🟡 C0–C3 ✅ · C4 open · C5/C6 🚫 blocked |
| 5 Bind the skills | 🚫 blocked on the SciAgent-toolkit refactor |
| 6 Retire the submodule | ⬜ |
| 7 Distribution | ⬜ |
| **8 Stabilize the surface** | 🟡 S1 ✅ S2 ✅ S3 ✅ S4 ✅ · **S5, S6 open** |
| **9 One CoReSh/GESECA run** | 🟡 G1 ✅ G2 ✅ G4 ✅ G5 half · **G3 open** |
| 10–11 Activity layer, WGCNA, uniform surface | 🚫 parked by decision |

### What Phase 8–9 delivered so far

- **`bulkirna_api()`** — 78 exports with lifecycle and signature-freeze as *independent* axes:
  46 stable, 11 experimental, and 21 deprecated, with 24 signature-frozen names. All deprecated
  names are removed in `1.0.0`. Its load-bearing test compares the registry against `NAMESPACE`
  in both directions at test time.
- **`bulkirna_stochastic()`** — the five functions that consume randomness, each seed argument,
  default and source of randomness. `write_session_provenance()` records `RNGkind()`.
- **`R/rng.R`** — the only place allowed to mutate RNG state, with a parse-tree test that fails
  if a fourth site appears, verified in both directions.
- **CoReSh** — `coresh_chunks/match/search/convergence/validate/loadings/sets`, p-values through
  `fgsea::geseca()`, `pct_var` bit-identical to hand computation on 12 of 12 real datasets, and
  the loading projection identical to the reference implementation on 63 of 63 comparable hits.
  Dataset identity is the `(gse, gpl)` pair: 1,635 of 42,465 human accessions appear on more than
  one platform, and on the one I measured the two platforms shared 19 of 50 top-loading genes.
- **`gs_coregulation()`** — GESECA on any expression matrix, `stat_type = "pct_var"`, `direction`
  `NA` because the statistic is unsigned.
- **`gs_test()`** keeps fgsea's `log2err` instead of discarding it.
- **A real-data fixture** carrying the three structures that have broken this package.

### The strongest evidence the port is correct

`HALLMARK_HYPOXIA` against the whole human compendium: **44,253 datasets in 87 seconds**, and
**7 of the top 10 are explicitly hypoxia or HIF experiments**, with VHL — the canonical HIF
degradation pathway — at rank 4. Validated against biology, not against itself.

---

## 4. Next immediate steps

1. **G5's other half** — the same query through <https://alserglab.wustl.edu/coresh>, compared
   accession by accession. Needs a browser, so it is the owner's step. It is the only remaining
   independent check on the CoReSh port.
2. **G3** — `gsdb_coresh()` on `.ref_path("coresh")`, recording the snapshot tag in the
   `gs_db` provenance that G2 added. Its gate cannot be byte agreement with DC-nexus's stored
   GMT: that artefact's input chunks no longer exist, and its companion provenance file recorded
   each query's score rather than each hit's. Gate it on the loading-level agreement plus a
   fresh end-to-end run.
4. **Finish Phase 4 — STING-JR first.** It is also the TF reference project, so migrating it
   puts its conventions in front of us before the activity layer is designed. Two unmigrated
   consumers are what keep the deprecated shims alive, so this is on the critical path to
   `v1.0.0`.
5. **S5, S6** — return-type audit and one vignette per layer.

**Then** reconsider Phases 10–11 with `03_DEFERRED.md` in hand.

---

## 5. Open items that block nothing

- The image still pins `bulkiRNA@v0.4.0` (`scbio-docker/docker/base/R/install_core.R:174`), so
  nobody using `scdock-r-dev:v0.5.13` has the GATOM reproducibility fix. A 150-minute rebuild
  for it is a judgement call.
- `qs2` is absent from the image; every CoReSh chunk reader needs it, and `R CMD check` now
  needs it too, since it is a declared `Suggests`. Until the image carries it, the check has to
  run with a scratch library mounted at `R_LIBS`.
- The consumer GMT that G2 was gated against is **not reproducible from any snapshot on disk**:
  it predates the earliest surviving chunk rewrite. Its companion `coresh_provenance.csv` is
  also wrong — 13 distinct `pctVar` values across 58 rows, one per query rather than per hit.
- `gatom_download_refs()` still defaults to the repo-relative `00_data/references/gatom`. It is
  the historical default and the frozen shim must keep it, but a downloader whose default
  destination is relative to the working directory is how the test suite ended up writing 16 MB
  into `tests/testthat/`.
- `10_gatom_modules.R`'s combined-network path still calls `gatom::` directly.
- `14839-DM-cGAS` is uncommitted, with the moved GATOM numbers in it. Committing in a tree
  holding live research data is the owner's call.
- The `msigdbr` re-run's 47,988 → 72,408 MSigDB rows are still unread — the largest unexamined
  result in the project.
- Two `scbio-docker` commit messages say "the integrator's"; fixing means force-pushing
  published history.
- `README.md` and `CONVENTIONS.md` still carry small pre-package traces.

---

## 6. How this work has been run, and what it cost

Implementation is delegated to `codex gpt-5.6-sol` in isolated clones; review to independent
Opus agents; gating and merging stay with the integrator. Two things about that are worth
carrying forward.

**Agents stopping on a false premise is the highest-value behaviour in the loop.** It has
happened six times and been right six times: an off-by-one export count, five shims with no
machine-readable successor, three formals assumed to exist, an undecided RNG policy, an
inconsistency between two briefs, and a clone made from a commit predating the work it was
meant to extend. Briefs say so explicitly, and should keep saying so.

**A review finding is not automatically a defect.** The most severe finding of the last review
argued from R's C sources that a bare `RNGkind()` initialises `.Random.seed`, making a branch
dead. One command showed it false in R 4.5.3. Its second half — that the branch had no test —
was right and more useful. **Check before acting; the check is usually cheap.**

**Real data finds what fixtures cannot.** The `NA`-Entrez defect and the `R CMD check` fixture
failure were both invisible to a green `devtools::test()`. Run the real thing early.
