# bulkiRNA — handoff and plan of record

**Updated:** 2026-08-19 · **Branch:** `feat/bulkirna-package` · **Last release:** `v0.5.0` tagged
**Gates, freshly run:** 1,697 tests passing · golden 20/20 · `R CMD check` 0/0/0 (vignettes rebuild)
**Surface:** 79 exports — 46 stable, 12 experimental, 21 deprecated · 24 signature-frozen · 7 stochastic

This is the entry point. Everything else is reachable from here.

| Document | What it owns |
|---|---|
| [2026-08-10-bulkirna-package/00_INDEX.md](2026-08-10-bulkirna-package/00_INDEX.md) | Phases 0–7 and their execution record |
| [2026-08-13-analysis-api-roadmap/00_ROADMAP.md](2026-08-13-analysis-api-roadmap/00_ROADMAP.md) | Phases 8–9. **§7–§16 are the execution record**, including every gate that was superseded and why |
| [2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md](2026-08-13-analysis-api-roadmap/01_REFERENCE_PROJECTS.md) | **The canonical previous implementation of each analysis, with paths, parameters and the reasoning.** Read this instead of asking where prior work lives |
| [2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md](2026-08-13-analysis-api-roadmap/02_CORESH_METHOD.md) | What CoReSh's score is, traced to its implementation |
| [2026-08-13-analysis-api-roadmap/03_DEFERRED.md](2026-08-13-analysis-api-roadmap/03_DEFERRED.md) | TF activity, PROGENy, WGCNA — parked, with findings intact |
| [2026-08-13-analysis-api-roadmap/04_NAME_AUDIT.md](2026-08-13-analysis-api-roadmap/04_NAME_AUDIT.md) | Why each public name and argument is spelled as it is |
| [2026-08-13-analysis-api-roadmap/05_RETURN_AUDIT.md](2026-08-13-analysis-api-roadmap/05_RETURN_AUDIT.md) | What every export returns, with the reasoned exceptions |
| [2026-08-11-coresh-extraction/00_PLAN.md](2026-08-11-coresh-extraction/00_PLAN.md) | CoReSh extraction, C0–C6 |
| [../adr/](../adr/) | The four architectural decisions |
| [../../../AGENTS.md](../../../AGENTS.md) | The rules an agent must follow in this repo |

---

## 0. Why this project exists

`RNAseq-toolkit` was a folder of scripts consumed as a git submodule and pulled in with
`source()`. 24 copies, each pinned to a different commit, so a fix in one reached nobody. Nothing
was tested. Nothing carried a version, so no result could name the code that produced it. **The
goal is not tidiness. It is that a number in a figure can be traced to a version of the code that
produced it, and that a defect fixed once is fixed everywhere.**

The second goal, stated by the owner, is why Phases 8–11 exist at all:

> *a thin clean API that just does the job the way I like it, without having to re-state to agents
> implementing analysis for me to go search for my previous experience reference repos — a
> zero-token architecture for my personal preferences on how bulk and pseudo-bulk should be done.*

That reframes the package as **the executable record of a set of methodological preferences**, and
it sets the acceptance bar in three ways:

1. **A function needing a paragraph of instructions to use correctly has failed.** The defaults
   *are* the preferences, documented once.
2. **A choice made the same way three times is a default, not a decision.** Three projects agreeing
   on `top = 500` for PROGENy is no longer a judgement call.
3. **Prose is not enough.** `bulkirna_api()` and `bulkirna_stochastic()` are enforced in code with
   tests; the pipeline conventions in `01_REFERENCE_PROJECTS.md` are prose and can drift the way the
   originals did. Where a preference matters it should become a default with a test, and the
   remaining prose should say plainly that it is prose.

### For an agent asked to run an analysis "the way it is usually done here"

Read in this order, and do not re-derive any of it:

1. `AGENTS.md` — the rules, the container invocation, the gates.
2. `vignettes/` — the order of operations per layer: `gene-sets`, `differential-expression`,
   `gatom`, `coresh`.
3. `01_REFERENCE_PROJECTS.md` §0 and §0b — the conventions that hold across every stage, and the
   house `1.x` ingest / `2.x` compute / `3.x` render pipeline shape with its config keys and output
   layout.
4. `01_REFERENCE_PROJECTS.md` §1–§4 — the canonical prior implementation of the specific analysis,
   with the project path, the parameters used, and where those projects disagreed with each other.
5. `bulkirna_api()` for what is safe to build on; `bulkirna_stochastic()` for what is random.

---

## 1. The pain points

The first four are original; the rest were found during the work, each measured rather than inferred.

**1. Silent wrongness that a green run cannot distinguish from a null result.** Ten confirmed:

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
| 9 | `bplapply()` switches `RNGkind()`, so seeding without pinning changes answers | `gs_test()` returned different p-values inside a parallel worker, no warning |
| 10 | `NA` Entrez ids in 0.7% of real CoReSh datasets | Passed the entire test suite; found by a compendium sweep |

Numbers 5–10 share one shape: **an empty result and a broken run were indistinguishable from inside
the code, so the broken run took the empty result's path.** In one case that purged 492 rows from two
master tables and exited 0.

**And the shape recurred inside the package, after being named.** `coresh_sets()` counted extractions
rather than sets, so a run where every set fell outside the size bounds returned a clean empty `gs_db`
with no message — the reference implementation's `skip:` hole moved one step downstream rather than
closed. **This is the failure mode to assume is present.**

**2. Machinery retyped per project, then drifted.** ~81% of the surveyed CoReSh lines were generic;
the engine was md5-identical in four checkouts, re-vendored by a `cp` in a comment.

**3. The same drift reproduced *inside* the package.** Seeding was implemented three times, three
ways. Then five species handlers disagreed about what a species is —
`gene_to_entrez(species = "Homo sapiens")` worked while `coresh_search(species = "Homo sapiens")`
errored, same session, same species. TF activity is next in line: two estimators and three network
provenances across three projects.

**4. Reference data fetched live, so a result depends on the day.** `msigdbr`'s cache, the GATOM
downloads, and a documented OmniPath breakage.

**5. Provenance that can be silently meaningless.** `.ref_path()` recorded
`basename(Sys.readlink(current))`; a `current` pointing outside its source, or being a real
directory, recorded a snapshot identifying nothing while the call succeeded. **And the opposite
failure:** `gsdb_coresh()`'s first version `dput()`d the queries and the whole hits table into the
provenance to satisfy a scalar validator — kilobytes of escaped R code that looked like completeness
and was unreadable. An honest gap a reader can see beats a blob that looks complete.

**6. The agent-facing documentation described the dead predecessor.** `AGENTS.md` — the one file
every assistant reads, via `CLAUDE.md` and `GEMINI.md` — was 225 lines about a script library,
including "there is no `library(RNAseqToolkit)`". Rewritten 2026-08-14.

**7. A document can record a method that provably has a gap.** `04_NAME_AUDIT.md` said its table was
"measured from the formals of all 57 live exports" — the technique that cannot see an argument
declared on an S3 method behind a generic's `...`, which is how `group` and `samples` went unaudited.
The counts were right; the counting was not. **Those are separate claims and only the first was
true.**

**8. A false premise attached to an enforcement-test exception does active work.** Two exceptions in
the argument audit claimed their spellings were "fixed by frozen signatures". None of the seventeen
exports involved was frozen and none had a consumer call site. A wrong row in a plan document is
merely wrong; **a wrong reason in a test holds open the gap the test exists to close**, under a
sentence that reads as settled.

**9. An enforcement test's success path is untested until the day it succeeds.** Removing the last
multi-spelling exception turned the suite red: the comparison had only ever run with a non-empty
allowlist, and at zero the two sides had different shapes. The same inversion is scheduled for
`v1.0.0`, when the 21 shims go and every shim allowlist becomes permanently empty. Checked early and
fixed; the helpers now carry "do not delete as dead code" and the reason.

---

## 2. The ADRs, with premise, rejected alternatives, and what each has since had to survive

**ADR-001 — the package version is the unit of reproducibility.** *Premise:* four version
authorities floated free and no result could name its code. *Decisive argument:* under
image-as-unit the `msigdbr` bug would have been **invisible** — tag unchanged, numbers changed.
*Rejected:* image-as-unit, which hides dependency drift and couples a scientific claim to an
artifact rebuilt for unrelated reasons; a lock file inside `bulkiRNA`, which a library cannot impose
on its consumers. *What it grew:* the stability contract is its missing half — a version says *that*
something changed, a lifecycle tier says *whether you were entitled to rely on it*. `v1.0.0` is now
a defined event: the release where the 21 shims go and `stable` means stable.

**ADR-002 — the package owns the master-table schema and validator, not the file.**
`neg_log_padj = -log10(pmax(padj, .Machine$double.xmin))`. *Premise:* whoever computes a column owns
its definition, or the definition lives somewhere untestable. *Vindicated:* `37.55` in the live
table, above the retired cap of 16, on real values. *Rejected:* toolkit ownership — a skill cannot be
`library()`d, so the schema reverts to prose, which is how it drifted three ways; per-project schema
files, which drift by default. *Held under pressure twice:* when `gs_test()` gained `log2err`,
`gs_to_master()` still returned exactly 14 columns; and when `gs_coregulation()` introduced an
unsigned statistic, the NES guard was left **untouched** — a `pct_var` result reaches the master table
only with `stat_as_nes = TRUE`, and a test pins that the guard fires by default. A master table whose
`nes` column silently holds three different statistics is the defect this project began with.

**ADR-003 — `Suggests` stays.** *Premise:* it is the machine-readable list `R CMD check`,
`install_deps` and the image build read. *Decisive argument:* 14 light `Imports` are why
`install_github` works on bare R in seconds. *Rejected:* collapsing into `Imports`, making every
optional feature a hard install blocker; core-plus-feature packages, cleanest isolation but a second
release train for one maintainer. *Amended twice on `Remotes:`* — unnecessary, then necessary when
`coresh` looked GitHub-only, then unnecessary again when `coresh` turned out to export nothing.
*Its seam proved itself:* adding `knitr` and `rmarkdown` for S6 immediately failed the
registry-versus-`Suggests` drift test, which is what that test is for. They are recorded as
**dev-only beside `testthat`**, because `bulkirna_check_deps()` answers "what must I install to use
this feature" and that answer never includes knitr. *Watch item:* WGCNA's dependency tree is the
first that makes the core/feature split arguable.

**ADR-004 — two reference-data tiers, bundled or refcache.** *Premise:* three mechanisms had grown
and a fourth environment variable was about to appear. Bundling MitoCarta stands: a source with no
API cannot be re-fetched reproducibly, so vendoring **is** the archive, and 168 KB makes the size
objection moot. *The insight:* nobody needs to own pinning, someone needs to own **recording**.
*Earning its keep three times:* CoReSh chunks resolve through it; the OmniPath networks are exactly
its case, fetchable in principle and unfetchable in practice at current versions; and `gsdb_coresh()`
now carries `snapshot = syn66227307_20260721` in the returned object's provenance rather than in a
log line, which is what the ADR asked for and had not yet been made to do.

**ADR-005 was proposed and withdrawn.** It was to settle ULM versus MLM for TF activity. With that
layer parked there is no decision to record, and an ADR for unscheduled work rots. It survives as an
open question in `03_DEFERRED.md`.

---

## 3. Where we are

| Phase | Status |
|---|---|
| 0 Inventory & freeze | ✅ |
| 1 Skeleton, internal `source()` removed | ✅ |
| 2 Prove the dev loop | ✅ against real data |
| 3 Install into the image | ✅ `scdock-r-dev:v0.5.13`, 16/16 optional deps |
| **4 Migrate heavy consumers** | 🟡 **1 of 3 — `14839-DM-cGAS` done; STING-JR and DC-nexus remain** |
| 4c CoReSh extraction | ✅ C0–C4 (C4 landed as G3) · C5/C6 🚫 blocked on the skills refactor |
| 5 Bind the skills | 🚫 blocked on the SciAgent-toolkit refactor |
| 6 Retire the submodule | ⬜ waits on Phase 4 |
| 7 Distribution | ⬜ |
| **8 Stabilize the surface** | ✅ **S1–S6 complete** |
| **9 One CoReSh/GESECA run** | ✅ **G1–G4** · 🚫 **G5's web-UI half needs a browser** |
| 10–11 Activity layer, WGCNA, uniform surface | 🚫 parked by decision, findings intact |

### What Phases 8–9 delivered

- **`bulkirna_api()`** — 79 exports with lifecycle and signature-freeze as *independent* axes,
  because 20 of the 24 frozen names are also deprecated and one column would have to lie about one or
  the other. Its load-bearing test compares the registry against `NAMESPACE` in both directions.
- **`bulkirna_stochastic()`** — the 7 functions that consume randomness, each seed argument, default
  and source. **The three default seeds deliberately disagree** — `coresh_*` uses upstream's literal
  `1L`, `gatom_module()` uses `42`, the fgsea adapter uses `123L` — because each matches numbers
  already published. Documented deviations, not untidiness.
- **`R/rng.R`** — the only place allowed to mutate RNG state, with a parse-tree test that fails if a
  fourth site appears, verified in both directions.
- **One species resolver** where five disagreed, widening every alias and narrowing none except
  `annotate_genes(species = NULL)`, which used to mean mouse silently.
- **One spelling per concept**: `top`/`top_n`/`n_top` → `top_n`, `color_palette`/`colours`/`palette`
  → `palette`, decided by the owner once the cost was measured at zero.
- **CoReSh, complete** — `coresh_chunks/match/search/convergence/validate/loadings/sets` plus
  `gsdb_coresh()`. `pct_var` bit-identical to hand computation on 12 of 12 real datasets; the loading
  projection identical to the reference implementation on **63 of 63** comparable hits; the provider
  identical to its own hand composition, provenance included. Dataset identity is the `(gse, gpl)`
  pair: **1,635 of 42,465** human accessions appear on more than one platform, and on the one
  measured the two platforms shared **19 of 50** top-loading genes.
- **`gs_coregulation()`** — GESECA on any expression matrix, `stat_type = "pct_var"`, `direction`
  `NA` because the statistic is unsigned.
- **Four vignettes**, two genuinely evaluated, and `R CMD check` at 0/0/0 with vignettes rebuilding.
- **A real-data fixture** carrying the three structures that have broken this package, with a
  generating script that reproduces it byte-for-byte.

### The strongest evidence the port is correct

`HALLMARK_HYPOXIA` against the whole human compendium: **44,253 datasets in 87 seconds**, and **7 of
the top 10 are explicitly hypoxia or HIF experiments**, with VHL — the canonical HIF degradation
pathway — at rank 4. Validated against biology, not against itself.

### Gates not met as written, recorded as superseded rather than quietly reinterpreted

| Gate | Why not met | What replaced it |
|---|---|---|
| G2: memberships match DC-nexus's stored GMT | 28 of 58 match; the chunks that produced that GMT no longer exist, established four ways | Loading projection identical to the reference on 63 of 63 hits |
| G4: a golden baseline added | That harness exists to prove **frozen legacy names** did not drift; a new function would also pin fgsea's estimator version into it | Agreement with the hand-computed GESECA formula, plus planted-signal and seed-relation tests |
| S6: vignettes runnable on the shipped fixture | `tests/fixtures/` and `data/` are both `.Rbuildignore`d | Gene sets on real `inst/extdata` sets; DE on declared synthetic data; GATOM and CoReSh honest and unevaluated |
| §5: run G5 before G2–G4 | G5's remaining half needs a browser | Gate each step against the reference implementation instead; G5 stays the only external falsifier |

---

## 4. Next immediate steps

1. **Phase 4 — STING-JR, then DC-nexus.** The critical path to `v1.0.0`: the 21 shims exist only
   because those two consumers still call the old names. **STING-JR first**, because it is also the TF
   reference project, so migrating it puts its conventions in front of us before the activity layer is
   designed. Migrating a consumer means committing in a live research tree, so it needs the owner's
   go-ahead per project.
2. **G5's other half** — the same query through <https://alserglab.wustl.edu/coresh>, compared
   accession by accession. Needs a browser, so it is the owner's step, and it is the only remaining
   independent check on the CoReSh port.
3. **Phase 6** — retire the submodule, once no consumer sources it.
4. **Then reconsider Phases 10–11** with `03_DEFERRED.md` in hand: TF activity, PROGENy, WGCNA. The
   sequencing argument for parking them was that a default in a package is load-bearing in a way a
   default in a script is not, and those three are where the methodological drift is worst. That
   argument is unchanged; what has changed is that the surface they would sit on has stopped moving.

---

## 5. Open items that block nothing

- The image still pins `bulkiRNA@v0.4.0` (`scbio-docker/docker/base/R/install_core.R:174`), so nobody
  using `scdock-r-dev:v0.5.13` has the GATOM reproducibility fix. A 150-minute rebuild is the owner's
  judgement call.
- `qs2` is absent from the image; every CoReSh chunk reader needs it, and since it is a declared
  `Suggests`, `R CMD check` needs a scratch library mounted at `R_LIBS` until the image carries it.
- The consumer GMT that G2 was gated against is **not reproducible from any snapshot on disk**. Its
  companion `coresh_provenance.csv` is also wrong — 13 distinct `pctVar` values across 58 rows, one
  per query rather than per hit.
- `gatom_download_refs()` still defaults to the repo-relative `00_data/references/gatom`. The frozen
  shim must keep it, but a downloader whose default destination is relative to the working directory
  is how the test suite once wrote 16 MB into `tests/testthat/`.
- **Both source-tree enforcement tests skip under `R CMD check`** — the RNG-ownership parse-tree test
  and the name audit — because an installed package has no `.R` files. A green check does not cover
  them; only `devtools::test()` does, and the skips say so.
- The argument-vocabulary test is **157 of 178 formals in two generic buckets**. It cannot rot, and it
  enforces snake_case with reasoned exceptions, but a new spelling can still be waved through as a
  singleton concept. Stated plainly rather than dressed up.
- `10_gatom_modules.R`'s combined-network path still calls `gatom::` directly.
- `14839-DM-cGAS` is uncommitted, with the moved GATOM numbers in it. Committing in a tree holding
  live research data is the owner's call.
- The `msigdbr` re-run's 47,988 → 72,408 MSigDB rows are still unread — the largest unexamined result
  in the project.
- Two `scbio-docker` commit messages say "the integrator's"; fixing means force-pushing published
  history.

---

## 6. How this work has been run, and what it cost

Implementation is delegated to `codex gpt-5.6-sol` in isolated clones; review to independent Opus
agents; **gating and merging stay with the integrator**, because the delegated agents have no Docker
socket and therefore cannot run a single gate. Every "tests pass" in this document was run here.

Four things are worth carrying forward.

**Agents stopping on a false premise is the highest-value behaviour in the loop.** Twelve times,
right twelve times: an off-by-one export count; five shims with no machine-readable successor; three
formals assumed to exist; an undecided RNG policy; an inconsistency between two briefs; a clone made
from a commit predating the work; a `gs_db` provenance premise that was false in three ways; a
`stat_type` vocabulary that would have had to be invented; a golden-case registration my own
read-only rule forbade; a fifth species handler I had missed; and vignette data that does not exist
in a built package. Briefs say "if a premise is false, stop and tell me", and should keep saying it.

**A review finding is not automatically a defect, and neither is a reviewer's method.** One review's
most severe finding argued from R's C sources that a bare `RNGkind()` initialises `.Random.seed`; one
command showed it false in R 4.5.3. Its second half — that the branch had no test — was right and
more useful. Later, checking a reviewer's own verification method found two public formals that
nobody's audit had ever seen. **Check before acting; the check is usually cheap.**

**Measure the claim, including your own prose.** Four findings across two reviews surfaced only
because a sentence was executed instead of read: the 88% singleton figure against a stated threshold,
seventeen non-frozen exports behind a "frozen signatures" exception, two unaudited S3 method formals,
and an empty-versus-empty comparison that had never run. The findings that took longest to surface
were the ones where a sentence sounded true.

**`devtools::test()` is not the gate.** It missed a semantic merge conflict between two branches that
each passed alone, stale `man/` twice, a fixture-backed test that must skip in a built package, and
an error-message rewrite that broke a matcher. **Run `R CMD check` on the built package after every
merge**, not only after an edit.
