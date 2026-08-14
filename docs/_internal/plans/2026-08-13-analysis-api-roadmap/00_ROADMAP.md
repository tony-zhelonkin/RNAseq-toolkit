# bulkiRNA: stabilize the surface, then one honest CoReSh/GESECA run

**Date:** 2026-08-13 · **Status:** plan of record for the work after Phase 7
**Supersedes:** the first draft of this file, which scheduled TF activity, PROGENy and WGCNA
as Phases 9–11. Those are **parked** — see [03_DEFERRED.md](03_DEFERRED.md).

**Companions**
- [01_REFERENCE_PROJECTS.md](01_REFERENCE_PROJECTS.md) — canonical reference implementations,
  their conventions, their defects.
- [02_CORESH_METHOD.md](02_CORESH_METHOD.md) — what CoReSh's score actually is, where it is
  implemented, and what running "their" analysis requires.
- [03_DEFERRED.md](03_DEFERRED.md) — the activity and network layers, deliberately not
  scheduled.

---

## 0. The decision that shapes this plan

The first draft of this roadmap proposed three new subsystems. The owner's call:

> *before we hop onto that more complexity we need to stabilize the surface of the API in its
> current form, maybe with standard coresh and geseca run, but without all the rest for now.*

That is the right call and it is worth saying why, because the argument is not only about
appetite.

**The package has 72 exports, 20 of which are deprecated shims, and no stability contract at
all.** Nothing says which names are safe to build on. Meanwhile 24 exports are frozen by an
explicit promise, 5 `coresh_*` functions are two days old and unproven against a consumer, and
13 live exports sit outside every layer prefix. Adding `tf_activity()` on top of that would
mean a consumer cannot tell a settled name from a provisional one — which is exactly the
failure the refactor was meant to end, reintroduced one level up.

**And the second reason is sharper.** The activity and network layers are where the
methodological drift is worst, which makes them the most valuable to unify and the most
expensive to get wrong. A default in a package is load-bearing in a way a default in a script
is not. They deserve to be built on a surface that has stopped moving.

So: **two phases, both small, both finishable.** Phase 8 makes the current surface something
you can depend on. Phase 9 completes one analysis end to end at full fidelity, which is what
proves the surface is actually usable rather than merely tidy.

---

## 1. Where we stand, measured

| Fact | Value |
|---|---|
| Version | `v0.5.0`, tagged and pushed |
| Exports | **72** — 52 live, **20 deprecated shims** |
| Live exports under a layer prefix | 39 (`gs_` 17, `gsdb_` 6, `de_` 6, `gatom_` 5, `coresh_` 5) |
| Live exports outside every prefix | **13** |
| Tests | 1060 passing, 0 failing, 5 skipped |
| Golden baselines | 20/20 |
| **`R CMD check`** | **0 errors, 0 warnings, 1 note** — and the note was a `NEWS.md` heading, now fixed |
| Consumers migrated | 1 of 3 (`14839-DM-cGAS`) |

The check result is better than I expected and it changes the character of Phase 8: this is
not a repair job. **The surface is technically sound and rhetorically unfinished.** What is
missing is a contract, not a fix.

The 13 unprefixed live exports, which is the concrete audit target:

| Group | Exports |
|---|---|
| I/O | `read_counts_matrix`, `read_metadata`, `ensure_dir`, `write_session_provenance` |
| DE construction | `build_dge`, `annotate_genes` |
| Gene identifiers | `gene_to_entrez`, `entrez_to_gene`, `filter_confounder_genes` |
| Presentation | `format_pathway_name`, `theme_bulki` |
| Package meta | `bulkirna_check_deps` |
| GATOM, misnamed | `download_gatom_references` |

---

## 2. Phase 8 — stabilize the surface

The goal is a package where **a name tells you what you are allowed to rely on**. No new
capability. Six steps, in order.

| Step | Work | Gate |
|---|---|---|
| **S1** | ✅ **done 2026-08-13.** `bulkirna_api()` ships the contract as a table. **Write the stability contract.** Three tiers: `stable` (the 24 frozen plus anything a migrated consumer depends on), `experimental` (the 5 `coresh_*`, the 3 gene-id functions — may change without a major bump), `deprecated` (the 20 shims, with a removal version). Publish it as a vignette and a machine-readable table, and make `bulkirna_api()` return it. | A test asserts every export carries exactly one tier, so a new export cannot be added without classifying it |
| **S2** | **Name audit on the 13 unprefixed exports.** Decide per name: keep as a deliberate top-level verb, or move under a prefix with a shim. `download_gatom_references` → `gatom_download_refs` is the clear one, because it already belongs to a family whose other five members are prefixed. `ensure_dir` is a utility that probably should not be exported at all. | Every decision recorded with a reason in the contract vignette; no silent renames |
| **S3** | **Argument-name consistency audit.** One spelling per concept across all 52 live exports: `species`, `db`, `contrast`, `seed`, `quiet`, `verbose`, `path`, `min_size`/`max_size`. Today `gs_test()` and `coresh_search()` do not agree on everything, and the species aliases accepted by `gatom_refs()`, `gsdb_msigdb()` and `gene_to_entrez()` were each written separately. | A test enumerates formals across exports and fails on a known-bad spelling; one shared `.species()` resolver |
| **S4** | ✅ **done with S1** — `removed_in` is `1.0.0` for all 20 shims, stated in `MIGRATION.md`. **Set the deprecation clock.** The 20 shims exist because 64 call sites could not migrate at once. Two of three consumers are still unmigrated, so they stay — but with a named removal version (`v1.0.0`) and a warning that says it. | `MIGRATION.md` states the removal version; a test asserts every shim warns |
| **S5** | **Return-type and error-message consistency.** Every compute function returns a tibble; every renderer returns a ggplot; every validation error names the fix. Mostly true already — this step is the audit that proves it and the tests that keep it true. | Tests assert the class of every export's return on the shipped fixture |
| **S6** | **`R CMD check` stays at zero notes, and vignettes build.** One vignette per live layer: gene sets, DE, GATOM, CoReSh. Each runnable on the shipped fixture with no network and no refcache. | `rcmdcheck` in the image, 0/0/0, vignettes built |

**What Phase 8 deliberately does not do:** rename anything frozen, remove any shim, or add any
function other than `bulkirna_api()`.

---

## 3. Phase 9 — one standard CoReSh and GESECA run

Small, and it is the acceptance test for Phase 8. See
[02_CORESH_METHOD.md](02_CORESH_METHOD.md) for the method and the fidelity gaps.

| Step | Work | Gate |
|---|---|---|
| **G1** | ✅ **done 2026-08-13.** Route `coresh_match(pvalues = TRUE)` to `fgsea::geseca()` — `center = FALSE`, `scale = FALSE`, `make.unique()`d rownames, stored `totalVar` kept for `pct_var`, `log2err` carried into the result. Restore p-value ordering in `coresh_search()`. Delete the guard. | p-values on real chunks; agreement with `gesecaCpp` within summed `log2err` on ≥7 of 8 datasets; both rankings reproduce the vignette's shape |
| **G2** | `coresh_loadings()` and `coresh_sets()` — the rest of C3. | Set names and memberships match DC_hum_verse's existing GMT |
| **G3** | `gsdb_coresh()` (C4), on `.ref_path("coresh")`, recording the snapshot tag. | Byte-level agreement on set contents |
| **G4** | `gs_coregulation()` — GESECA as a first-class verb on any expression matrix, returning a `gs_result` with `stat_type = "pct_var"`. **Its own verb, not a `gs_test()` method**: it takes a matrix where `gs_test()` takes ranks, and a shared signature would lie about the input. | Existing `gs_plot_*` renderers work unchanged; a golden baseline added |
| **G5** | 🟡 **half done** — the `pct_var` sweep ran and validated against GEO; the web-UI comparison is outstanding. The end-to-end reference run: `HALLMARK_HYPOXIA` against the mouse and human compendia, compared against the web UI at <https://alserglab.wustl.edu/coresh>. | Top accessions agree with the web UI. Two independent implementations agreeing beats any unit test |

**G5 is the real gate.** Everything else is internally consistent by construction; only G5 can
tell us the port is correct rather than merely self-consistent.

---

## 4. The ADRs, and what these two phases ask of them

**ADR-001 — the package version is the unit of reproducibility.** *Premise:* four version
authorities floated free and no result could name its code. *Decisive argument:* under
image-as-unit the msigdbr bug would have been invisible — tag unchanged, numbers changed.
*What Phase 8 adds:* a stability contract is the missing half of this ADR. A version number
tells you *that* something changed; a tier tells you *whether you were entitled to rely on it*.
`v1.0.0` becomes meaningful — the release where the shims go and `stable` means stable.
*Rejected:* image-as-unit; a lock file inside a library.

**ADR-002 — the package owns the master-table schema.** *Premise:* whoever computes a column
owns its definition, or the definition lives somewhere untestable. *What these phases add:*
nothing. Untouched, and that is the point — no new entity classes until the deferred work
starts. *Rejected:* a schema per entity class.

**ADR-003 — `Suggests` stays.** *Premise:* it is the machine-readable list `R CMD check`,
`install_deps` and the image build read; 14 light `Imports` are why `install_github` works on
bare R in seconds. *What these phases add:* `qs2` and `BiocParallel` already landed, and
**that is the whole cost of Phases 8 and 9.** Parking the activity and network layers parks
`decoupleR`, `OmnipathR` and the entire WGCNA dependency tree with them, which is a large part
of why the owner's sequencing is correct. *Rejected alternative, still rejected:* core plus
feature packages.

**ADR-004 — two reference-data tiers.** *Premise:* three mechanisms had grown and a fourth
environment variable was about to appear; nobody needs to own pinning, someone needs to own
recording. *What these phases add:* G3 is its second real exercise, and the snapshot tag
`syn66227307_20260721` has to reach the `gs_result`, not just the log.

**Proposed ADR-005 is withdrawn for now.** It was to settle ULM versus MLM. With the activity
layer parked there is no decision to record, and writing an ADR for work that is not scheduled
is how plans rot. It moves to the deferred document as an open question.

---

## 5. What to do first

1. **S1**, the stability contract. Everything else in Phase 8 is an audit whose findings need
   somewhere to be recorded, and this is that place.
2. **G1**, because it deletes a guard built on a premise I got wrong, and because the p-value
   ranking is the one upstream calls *more specific*.
3. **S2 and S3** together — both are name decisions and both touch the same 13 exports.
4. **G5 early, not last.** It needs no new code beyond G1 and it is the only check that can
   falsify the port. Run it as soon as G1 lands, before G2–G4 build on top.
5. Then finish **Phase 4** — STING-JR, then DC-nexus. Two unmigrated consumers are what keep
   the 20 shims alive, so this is on the critical path to `v1.0.0`.

**Open, none of it blocking:** the image pin still says `bulkiRNA@v0.4.0`; `qs2` is absent
from `scdock-r-dev:v0.5.13` and every chunk reader needs it; `10_gatom_modules.R`'s combined
path still calls `gatom::` directly; `14839-DM-cGAS` is uncommitted with the moved GATOM
numbers in it; the msigdbr re-run's 47,988 → 72,408 MSigDB rows are still unread.

---

## 6. Tradeoffs, stated plainly

**The cost of stabilizing first.** Phase 8 ships no new capability. Its entire output is a
contract, some renames, some tests and four vignettes — and the work it makes possible is work
that could have started immediately instead. If the activity layer were needed next week, this
sequencing would be wrong.

**The cost of not stabilizing first**, which is why it is right anyway. Three consumers are
mid-migration, and every one of them is currently free to depend on a two-day-old
`coresh_match()` signature with no way to know it is provisional. The 20 shims have no removal
date, so they are permanent by default. Adding three subsystems to that would make the
package's own surface the next thing needing a refactor.

**The alternative — stabilize *and* build in parallel.** Rejected. The audits in S2, S3 and S5
are exactly the kind of work that a moving surface invalidates; doing them while adding
`tf_activity()` means doing them twice.

**The alternative — skip Phase 9 and stabilize only.** Tempting, and rejected for one reason:
an audit with no consumer proves nothing. G5 is the step that turns "the names are consistent"
into "someone ran a real analysis through them and the answer matched an independent
implementation".

---

## 7. S1 and G1 as executed (2026-08-13)

Two stages, each written by a delegated agent in an isolated clone, gated here, and read by an
independent reviewer. Package state: **1212 tests passing, golden 20/20, `R CMD check` 0/0/0,
73 exports.**

### S1 — the contract exists

`bulkirna_api()` returns one row per export with `name`, `layer`, `lifecycle`, `frozen`,
`superseded_by`, `removed_in`. 45 stable, 8 experimental, 20 deprecated, 24 frozen, all shims
removed in `1.0.0`.

**The two-axis design was the right call and the numbers prove it:** 20 of the 24
signature-frozen names are also deprecated, so a single tier column would have had to lie about
one or the other. `frozen ∖ deprecated` is exactly `{build_dge, download_gatom_references,
ensure_dir, format_pathway_name}` — the four original toolkit names that remain the real API.

The load-bearing test compares the registry against the `NAMESPACE` in both directions at test
time, so an export cannot be added without classifying it. The reviewer verified this is not
satisfiable vacuously, and flagged that `getNamespaceExports()` would break it under
`load_all(export_all = TRUE)`; it now parses `NAMESPACE`, matching
`test-namespace-hygiene.R`.

**Two defects worth recording.** The hardcoded `superseded_by` strings for the message-only
shims were being written unconditionally, so a shim later gaining a machine-readable target
would have drifted from its own warning in silence — now a build failure. And two shims named
`gs_db` and `gs_result` as replacements, **neither of which is exported**: callers were being
sent to functions they cannot reach. `list_to_term2gene` now names `gsdb_register()`;
`empty_gsea_tibble` admits there is no exported constructor, in the manner `save_gsea_log`
already used.

### G1 — the p-value path works

`coresh_match(pvalues = TRUE)` routes to `fgsea::geseca()`. `pct_var` still divides by the
stored `totalVar` and is unchanged. `log2err` is now a column. `coresh_search()` ranks by
ascending `p_value` when p-values are requested. No `:::` anywhere.

**The RNG finding, which cost the most to diagnose.** `geseca()` takes no seed and draws its
internal C++ seeds from R's RNG, and `BiocParallel::bplapply()` switches `RNGkind()` to
`"L'Ecuyer-CMRG"` inside the task — **even with `SerialParam()`**, so this is not about cores at
all. `set.seed(seed)` alone therefore gives different answers inside and outside a parallel
call. `.coresh_with_seed()` pins the generator as well as the seed and restores the caller's
state. Verified on real data: 1,000 datasets, `n_cores` 1 and 4, identical p-values.

**The defect the test suite could not find.** 0.7% of real datasets carry `NA` Entrez ids — 11
of 1,500, one of them 9,998 rows of 10,000. `geseca()` rejects `NA` rownames and
`make.unique()` handles them erratically, returning `NA` for the first and the string `"NA.1"`
for the second. The full suite passed before this surfaced; only the compendium sweep found it.
**A unit test on a hand-built fixture cannot discover what real data contains**, which is the
argument for G5 running early rather than last.

Five further defects came from the review, none of which fire on the happy path: an unguarded
`[[1L]]` where `geseca()` can return no rows, a duplicated query id silently decoupling
`pct_var` from `p_value`, locale-dependent tie-breaking, `.Random.seed` left behind in a fresh
session, and a per-dataset warning when restoring a caller's legacy `sample.kind`. Three tests
were also weaker than they looked; the ranking one now scales `totalVar`, which moves `pct_var`
without touching the GESECA input and so forces the two orderings to disagree.

### G5, half done — and it is the strongest evidence we have

`HALLMARK_HYPOXIA`, 200 genes, against the **whole human compendium: 44,253 datasets in 87
seconds** on 12 cores, matching upstream's claimed scale. Then the check the vignette itself
prescribes, GEO titles for the top 10:

| rank | accession | title |
|---|---|---|
| 1 | GSE131379 | Hypoxia effects on the transcriptome of HeLa cells |
| 2 | GSE198308 | SRSF6-GFP HeLa under normoxic and hypoxic conditions |
| 3 | GSE59449 | MCF7 breast epithelial cells in normoxic and hypoxic conditions |
| 4 | GSE239892 | Von Hippel Lindau tumor suppressor controls m6A-dependent expression |
| 5 | GSE196274 | Regulation of HIF1α signaling in ER positive breast cancer |
| 6 | GSE179327 | mRNA transcriptome of hypoxia with shSETDB1 in HeLa |
| 10 | GSE71401 | Tumor hypoxia causes DNA hypermethylation by reducing TET activity |

**Seven of the top ten are explicitly hypoxia or HIF experiments**, with VHL — the canonical HIF
degradation pathway — at rank 4. That is the port validated against biology rather than against
itself.

**Outstanding:** the same query through the web UI at <https://alserglab.wustl.edu/coresh>,
compared accession by accession. It needs a browser, so it is the owner's step.

### Left as they were

The mouse sweep, and the `snapshot` field reading `coresh` rather than the tag in one test run —
an artifact of the scratch symlink used for that run, not of `.ref_path()`, which reported
`syn66227307_20260721` correctly whenever pointed at a real refcache layout.

---

## 8. The nick-trimming pass (2026-08-13)

Two fan-outs plus a follow-up round, each reviewed. **1325 tests passing, golden 20/20,
`R CMD check` 0/0/0.** The point of the pass was the deeper causes, not the symptoms.

### Deeper cause 1 — a concern with no owner drifts, even inside this package

Seeding a stochastic library call was implemented **three times, three ways**:

| site | restored the caller's stream | pinned the generator |
|---|---|---|
| `.coresh_with_seed()` | yes | yes |
| `.gs_fgsea()` | yes | **no** |
| `gatom_module()` | **no** | **no** |

This is the founding pain point of the whole refactor — generic machinery retyped and drifted —
reproduced inside the package that exists to end it, in the one area where the symptom is silent:
not a crash, a different draw.

**Both failure modes were measured, not inferred.** `BiocParallel::bplapply()` switches
`RNGkind()` to `"L'Ecuyer-CMRG"` inside the task **even with `SerialParam()`**, so `set.seed()`
alone drives a different generator there. `gs_test()` on one fixed input returned
`0.8963006 0.6011883 0.5638954` in the parent and `0.8981157 0.6039905 0.5644991` inside
`bplapply` — same seed, different numbers, no warning. And `gatom_module()` overwrote the
caller's stream, the exact defect a comment in `gs-test.R` describes being fixed *there*.

`R/rng.R` is now the only place in the package allowed to seed. `gatom_module()` still seeds
twice, in two independent wraps, because collapsing them would change every module ever computed;
the test that guarantees nothing moves is `.with_pinned_seed(42, rnorm(5))` being identical to
`{set.seed(42); rnorm(5)}`. After the change `gs_test()` returns the **parent's** values in both
contexts, so parallel agrees with serial without serial moving.

### Deeper cause 2 — a hand-built fixture contains only what its author imagined

The `NA` Entrez ids in 0.7% of real datasets passed the entire suite and were caught by a
compendium sweep. That is a blind spot in the *method* of testing, not in any one test.

`tests/fixtures/coresh-chunk-micro.rds` is four **real** dataset objects, 300 genes each, carrying
between them every structure that has broken this package: repeated ids, missing ids, and a
matrix whose columns are principal components. Columns are deliberately not subset — the first
attempt cut them to 8 and silently destroyed three of the four properties, which is the same
mistake one level up. `make_coresh_micro.R` reproduces the file byte-for-byte; it was the only
fixture in the tree without a generating script.

### Provenance that can be silently wrong

`.ref_path()` recorded `basename(Sys.readlink(current))`. Two ways that is meaningless while the
call succeeds: `current` pointing outside its source directory records whatever that path is
called, and `current` being a **real directory** records the literal string `"current"`, which
reads as plausible in a `gs_result`'s provenance. Both now warn, once per source, re-armed by
`.clear_ref_resolutions()`. A nested-but-internal snapshot stays silent, and the real refcache
still records `syn66227307_20260721` with no warning.

### On the review

The reviewer's most severe finding was **wrong** — it argued from R's C sources that a bare
`RNGkind()` initialises `.Random.seed`, making a branch dead. Measured in R 4.5.3: `exists()` is
`FALSE` before and after. The lines were left alone. Its second half was right and more useful:
that branch had no test. **A finding worth checking is not the same as a finding worth acting
on**, and the check cost one command.

Three tests were also asserting nothing: a tautology over a whitelist, a fixture branch that
built `c(NA, query)` and stripped the `NA` back out before calling anything, and an
`expected_size` that reduced to `length(query)` so the duplicate-collapsing path was never
reached. All three now assert the property they were named for.

### What the agents got right, and what that cost

Codex stopped and reported rather than guessing **six times across the session**, and was right
every time: an off-by-one in my export count, five shims with no machine-readable successor, three
formals I assumed existed, an undecided RNG policy, an inconsistency between two of my own briefs,
and a clone I had made from a commit that predated the work it was supposed to build on. That
last one was my error and cost a full round. One round died on a `bwrap` sandbox failure and wrote
nothing; `setsid` fixed it.
