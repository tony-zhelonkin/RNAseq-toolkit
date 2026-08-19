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

**The package has 74 exports, 20 of which are deprecated shims, and no stability contract at
all.** Nothing says which names are safe to build on. Meanwhile 24 exports are frozen by an
explicit promise, 5 `coresh_*` functions are two days old and unproven against a consumer, and
13 inherited live exports sit outside every layer prefix. Adding `tf_activity()` on top of that would
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
| Exports | **74** — 54 live, **20 deprecated shims** |
| Live exports under an analysis-layer prefix | 39 (`gs_` 17, `gsdb_` 6, `de_` 6, `gatom_` 5, `coresh_` 5) |
| Other live exports | **15** — the 13 inherited S2 targets plus `bulkirna_api()` and `bulkirna_stochastic()` from S1 |
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
| **S3** | **Argument-name consistency audit.** One spelling per concept across all 54 live exports: `species`, `db`, `contrast`, `seed`, `quiet`, `verbose`, `path`, `min_size`/`max_size`. Today `gs_test()` and `coresh_search()` do not agree on everything, and the species aliases accepted by `gatom_refs()`, `gsdb_msigdb()` and `gene_to_entrez()` were each written separately. | A test enumerates formals across exports and fails on a known-bad spelling; one shared `.species()` resolver |
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

---

## 9. Determinism made enforceable (2026-08-14)

Three fan-outs. **1385 tests passing, golden 20/20, `R CMD check` 0/0/0.** §8 fixed the
mechanism; this pass fixed the policy around it and made the invariant enforceable.

### What was still wrong after the mechanism was fixed

`.with_pinned_seed()` was correct and used at all three sites, and none of that was visible or
defended. Measured on the tree at the time:

- **The stochastic surface was undiscoverable.** Four exports took a seed; `gs_test()` was
  stochastic too but accepted `seed` only through `...`, so it did not appear in its own
  signature. Nothing let a user ask which functions are random.
- **`write_session_provenance()` recorded nothing about RNG** — 76 lines, zero mentions of seed,
  generator or `RNGkind` — in a package whose ADR-001 makes the version the unit of
  reproducibility. `RNGkind` is the single most useful line it could carry, since it is what
  silently changed answers inside `bplapply()`.
- **`gs_test()` discarded `log2err`.** fgsea returns it natively; the adapter dropped it. So the
  same estimator family was treated two ways in one package — CoReSh reported its uncertainty and
  `gs_test()` threw it away. Nothing looks wrong at either call site, which makes it a more
  interesting kind of drift than a retyped function.
- **Nothing prevented a fourth seeding site** from appearing.

### What was deliberately not changed

**The three default seeds stay different.** `coresh_*` uses `1L` because that is the literal value
upstream's reference implementation passes, `gatom_module()` uses `42`, the fgsea adapter uses
`123L`, and `run_gsea()` is signature-frozen. Unifying them would silently move every number
already published. They are **documented deviations, not an accident to tidy** —
`bulkirna_stochastic()` records each one with its reason, which is the fix.

**The legacy shims keep their historical shape.** `log2err` moved two golden baselines through
`run_gsea()` and `normalize_gsea_results()`. Rather than re-capture, `.drop_legacy_extras()` stops
the addition at the shim boundary: a shim exists so unmigrated callers need not change, two of
three consumers are unmigrated, and such a caller may index positionally. Anyone who wants the
uncertainty bound wants `gs_test()`. Golden returned to 20/20 with no baseline touched.

**`gs_to_master()` still returns exactly the 14 ADR-002 columns**, verified rather than assumed —
the whole `neg_log_padj` defect this project began with was a column-union accident.

### The enforcement test, and why it has teeth

It walks the parse tree of every file under `R/` and rejects any RNG-state **mutation** outside
`R/rng.R`: `set.seed()`, `RNGkind()` called with arguments, assignment to or removal of
`.Random.seed`. A bare `RNGkind()` is a **read** and stays legal, which matters because
`write_session_provenance()` queries it on purpose.

**That distinction came from an agent refusing the brief.** I first wrote the rule as "`RNGkind`
appears only in `R/rng.R`", and it stopped to point out that `R/data-io.R` would fail it for a
read-only provenance query. It was right, and the corrected rule — mutation, not access — is the
one that states the invariant that actually matters. **That is the sixth time an agent on this
package has stopped on a false premise in one of my briefs.**

Verified in both directions rather than trusted: a `set.seed()` and an argument-bearing
`RNGkind()` planted in `R/utils.R` fail the test and are named individually; a bare `RNGkind()`
read added to the same file passes. The test also asserts `R/rng.R` *does* contain a `set.seed()`,
so it cannot rot into a pass over an empty set.

### The classification loop, which is the durable part

Every name `bulkirna_stochastic()` declares must be either exercised in the loop — same seed twice
identical, two seeds different, the caller's stream and `RNGkind()` untouched — or listed as
covered elsewhere with a reason. **Declaring a sixth stochastic function turns the suite red until
somebody classifies it.** Verified by declaring `gs_score` and watching six tests go red.

That is the property worth having: coverage cannot silently lag the registry, and the registry
cannot silently lag the code, because a separate test derives the registry from `formals()` in both
directions.

### Two lessons about writing tests for stochastic code

- **A tiny fixture proves nothing.** The first `gs_test` case used a 100-gene ramp and three sets;
  its p-values saturate, so two seeds compared equal and the test passed while asserting nothing.
  It now uses 2,000 genes and 16 sets to keep the estimator in the range where the seed matters.
- **State must be captured after the fixture is built.** Capturing before means the fixture's own
  `set.seed()` is mistaken for the function under test disturbing its caller. The comparison helper
  now brackets only the calls.

Three vacuous assertions were also found in existing tests — two serial-versus-`SerialParam()`
comparisons that would pass on empty p-value vectors, one `all()` that would pass on zero rows —
and one in the new code: `paste0()` recycles a zero-length vector to `""` when another argument is
longer, so the offender message built `": "` on the passing case and failed the test it existed to
explain.

---

## 10. S2, S3 and G2 as executed (2026-08-18)

Two stages in parallel, each written by a delegated agent in an isolated clone, each gated and
reviewed here. Package state: **1,465 tests passing, golden 20/20, `R CMD check` 0/0/0, 77
exports** — 46 stable, 10 experimental, 21 deprecated.

The delegated agents had no Docker socket, so **every gate in this section was run by me**, not by
the agent that wrote the code. Both agents were told to claim no gate as passing, and neither did.

### S2 — the 13 unprefixed exports, and only one moves

`04_NAME_AUDIT.md` records all thirteen decisions with reasons. `download_gatom_references` →
**`gatom_download_refs`**, with the frozen name kept as the 21st deprecated shim: five other
members of its family are prefixed and it was the only one that was not.

`ensure_dir` stays `stable`, and that decision changed on evidence — the consumer inventory still
lists two unmigrated callers, so un-exporting it would break a project mid-migration. The other
eleven read correctly as top-level verbs and are recorded as deliberate keeps rather than
oversights.

### S3 — five species handlers, disagreeing about what a species is

I briefed four sites. The agent found a **fifth** and stopped to say so before writing anything:
`annotate_genes()` in `R/dge.R` had its own `match.arg()` plus two `switch()` calls, accepting only
the scientific spellings.

| Site | Accepted before |
|---|---|
| `.gene_id_species()` | six aliases, case- and separator-insensitive |
| `.coresh_species_code()` | four exact, **case-sensitive**; rejected `"Homo sapiens"` |
| `.gatom_species()` | aliases, but its error named only the scientific spellings |
| `.gsdb_species_label()` | any non-empty string |
| `annotate_genes()` | scientific spellings only |

So `gene_to_entrez(species = "Homo sapiens")` worked while `coresh_search(species = "Homo
sapiens")` errored — same package, same session, same species.

`R/species.R` is now the single owner, in the way `R/rng.R` owns seeding. Every alias any of the
five accepted is still accepted; `annotate_genes()` keeps the partial scientific-name matching
`match.arg()` gave it for free, and `.gsdb_species_label()` keeps a path for a user-supplied
custom species. **One visible consequence, recorded in `NEWS.md`:** a recognised alias is now
normalised into the `gs_db` `species` attribute, so `gsdb_register(species = "mouse")` records
`"Mus musculus"`.

### The landmine the rename exposed, which is the most useful finding of the pass

Making the downloader a shim turned one existing test into a 16 MB download.

`test-api.R` loops over every deprecated export and calls it **bare** to prove it warns. The
downloader's defaults are `dir = "00_data/references/gatom"`, `species = "Mus_musculus"`,
`networks = c("kegg", "combined")` — so the suite fetched real reference files from
`artyomovlab.wustl.edu` into a **repo-relative** directory. `test-api.R` sorts before
`test-gatom.R`, and `gatom_refs()` searches that same relative default third, so an assertion that
an empty directory produces a missing-file error found real files and **returned successfully
instead of erroring**.

Three separate defects in one failure: a test reaching the network unskipped, a test writing into
the working directory, and an unrelated assertion in another file inverting its result because of
it. And the general shape: **a generic "call every shim" loop is only safe while every shim happens
to be harmless when called bare**, which was a property nobody was maintaining.

The loop now evaluates only the `.Deprecated()` call, so it still proves the warning fires while
the delegated work is unreachable. `test-gatom.R` mocks the existing `.gatom_search_dirs()` seam
rather than depending on what is on disk. And `helper-state.R` plus `test-zzzz-workdir.R` fail the
suite if any test leaves a new entry in the working directory, with `testthat::set_state_inspector()`
naming the individual test that leaked. `_snaps/` and `Rplots.pdf` are excluded as testthat's and
the graphics device's own, with the reason written next to the exclusion.

An audit of the other twenty shims found no second one that writes or reaches the network.

### G2 — the set-building layer, verified against the reference rather than itself

`coresh_loadings()` and `coresh_sets()` port `extract_gene_loadings.R`. Four things changed because
a package is not a script: the package-local chunk index instead of a second index memoized into
`options()`; provenance carried as a table instead of recovered by re-parsing `CORESH_<query>_<gse>`
set names; Jaccard deduplication keeping the higher-ranked hit by an explicit rule, so the result
no longer depends on input order; and a run that fails everywhere stopping rather than reporting a
clean empty result.

**The gate: on 63 of 63 comparable real hits, `coresh_loadings()` returns the identical top-50
Entrez ids in the identical order as the reference formula computed inline.** Twelve hits were not
comparable — the accession is absent from the current snapshot, or fewer than three query genes are
present.

`gs_db()` gained `provenance` and `set_provenance`. The agent stopped and reported that my brief's
premise was false — `gsdb_msigdb()` records no provenance, `gs_db()` had no such field, and
`[.gs_db` rebuilds from a fixed attribute list and would silently drop anything else. The decision:
make it generic and make it survive subsetting, with `set_provenance` subset to the sets that
remain. A snapshot tag that vanishes on the first `db[1:10]` is worse than no tag, because the
object still looks provenanced — the same failure `.ref_path()` had when it recorded the literal
string `"current"`.

### Against the consumer's stored GMT: 58 sets, and why 24 differ

| Outcome | Sets |
|---|---|
| Identical | 28 |
| Same membership, different order | 2 |
| Different membership | 24 |
| Missing from ours / extra in ours | 1 / 1 |

The missing/extra pair is the **deliberate dedup rule change**: two sets overlapping at Jaccard
0.923, where the reference dropped the later one and we keep the higher-ranked one.

The 24 are **not** attributable to the port, and this took four measurements to establish rather
than assert:

1. **Not the snapshot.** `syn66227307_20260430_migrated` and `syn66227307_20260721` produce
   byte-identical output. I ran the whole comparison against both.
2. **Not a boundary shift.** Genes the stored GMT ranked in its top 50 rank **59 to 8,933** today.
   No cutoff change does that.
3. **Not the symbol mapping.** On one differing set, the reference formula run on today's chunks
   agrees with `coresh_loadings()` on 50 of 50 and with the stored GMT on 32 of 50. The divergence
   is upstream of both implementations.
4. **Not the compendium score either, and the stored provenance cannot arbitrate.**
   `coresh_provenance.csv` holds **13 distinct `pctVar` values across 58 rows** — one per query,
   repeated across all of that query's hits. It recorded the query's score, not the hit's, so it
   cannot be compared against anything.

The GMT is dated 2026-04-25 and the earliest surviving snapshot is the `_migrated` rewrite of
2026-04-30. **The chunks that produced it no longer exist**, so it is not reproducible from
anything on disk — by us or by the script that wrote it.

### A compendium defect found on the way, not yet fixed

**1,635 of 42,465 human GSE accessions appear more than once**, the same accession on a different
platform (`gplId`), sometimes in a different chunk file. Both the reference script and
`coresh_loadings()` take the first match by `gseId` alone, so *which dataset you get depends on
chunk file ordering*. `coresh_chunks()` already returns a `gpl` column and `coresh_sets()` already
orders by it when present, so the index is not the problem — the lookup is. A GSE is not a unique
key in this compendium and the API currently pretends it is. Scheduled below rather than fixed
here, because it changes a signature.

### The merge, and the conflict git could not see

The two stages merged with three textual conflicts, all counts. The fourth was invisible to git
because the two sides touched **different files**: `R/coresh-sets.R` called `.gene_id_species()`,
which the audit had deleted in favour of `.species()`. `devtools::test()` passed on each branch
*and on the merge*, because `load_all()` had both definitions in scope. Only `R CMD check` on the
built package caught it. **Parallel stages need the built-package gate, not the fast one.**

---

## 11. The platform key and G4 (2026-08-18)

Two more stages, same shape. **1,572 tests passing, golden 20/20, `R CMD check` 0/0/0, 78
exports.** Every gate here was run by me; neither agent had a Docker socket and neither claimed
otherwise.

### A GSE accession was never a unique key, and the numbers are large

§10's finding, fixed. `coresh_loadings()` gains an optional `gpl`, `coresh_sets()` threads it into
loading extraction and into per-set provenance, and `coresh_convergence()` — which had the same
defect and which I had not spotted — now requires and groups by `(gse, gpl)`. `coresh_search()`
already scored every platform as its own row, and `coresh_validate()` performs no lookup, so
neither changed.

**Warning rather than error, on ambiguity.** An error would break `coresh_sets()` on any `top_hits`
table produced before this change, and `coresh_search()` did not emit `gpl` until now. So an
ambiguous accession warns, names every available platform, states which one was used and why, and
records that platform in the result and in the provenance. The choice is radix-first, so it does
not depend on chunk-file order, object order within a chunk, or locale.

**Measured on real data, which is what makes this worth the change.** GSE100112 is one of the 1,635.
Its two platforms:

| platform | chunk | top-5 Entrez by absolute loading |
|---|---|---|
| GPL17556 | `chunk_1` | 5996, 8553, 230, 5210, 54541 |
| GPL11154 | `chunk_36` | 5996, 339122, 8553, 6015, 55814 |

**Their top-50 sets share 19 genes of 50.** First-match was choosing between two materially
different answers by file order, silently. Verified also: a single-platform accession stays quiet,
a requested platform that is absent errors naming what is available, and asking for both platforms
of one accession yields two sets with `GPL`-suffixed names rather than one silently discarded as a
duplicate.

One honest limit, worth writing down. `coresh_loadings()` sees one chunk file, so it can only detect
ambiguity *within* that chunk. Cross-chunk ambiguity — which is the common case, since the duplicate
platforms often live in different files — is detected by `coresh_sets()`, which holds the index. A
direct `coresh_loadings()` call on a cross-chunk duplicate therefore does not warn.

### G4 — `gs_coregulation()`

GESECA as a first-class verb: an expression matrix in, a `gs_result` with
`stat_type = "pct_var"` out, `method = "geseca"`, `log2err` retained, seeded through
`.with_pinned_seed()` and declared in the stochastic registry. Rows are centred by default, because
GESECA's variance-along-a-direction is only the intended quantity on centred data and a general
matrix — unlike a CoReSh chunk — does not arrive that way.

Two decisions were mine, and the agent stopped for both rather than inventing them:

- **`pct_var` joins `gs_stat_types()`.** Adding a value to a controlled vocabulary is a widening.
- **`direction` is `NA`, not `"up"`.** The statistic is unsigned. `gs_direction()` on a strictly
  positive statistic would have labelled every row `"up"` — a claim the method does not make, and
  one a filter written `direction == "up"` would silently accept.

**That `NA` found a live defect in `gs_top()`.** It built its grouping key with `paste()`, which
renders `NA` as the string `"NA"`, so an unsigned result would have been grouped together with any
row genuinely carrying that text. It now uses `interaction(exclude = NULL)`.

`gs_to_master()`'s NES guard is deliberately untouched: a `pct_var` result reaches the master table
only with `stat_as_nes = TRUE`, and a test pins that the guard fires by default. A master table
whose `nes` column silently holds three different statistics is one of the defects this project
began with.

### The golden baseline G4 promised, and why it is not there

`capture_golden.R` says what it is for in its own header: every case calls a **frozen legacy name**,
so the goldens prove the refactor did not change behaviour reached through the shims. A new
experimental function is outside that purpose, and adding it would also pin fgsea's estimator
version into a baseline that exists to detect *our* drift.

So **the roadmap's G4 gate is not met as written, deliberately.** What stands in its place is
stronger for a Monte-Carlo function: a test asserting agreement with the centred, unscaled GESECA
formula computed by hand, invariance to row offsets, a planted coregulated set outranking a
scrambled one, and the usual seed relations including parent-versus-`bplapply` agreement. A
snapshot proves a number has not changed; an independent formula proves the number is right.

### Two process findings from this pass

**A textual merge can be clean and still wrong.** `R/coresh-sets.R` called `.gene_id_species()`,
deleted by the concurrent audit in a different file. `devtools::test()` passed on both branches and
on the merge, because `load_all()` had both definitions in scope; only `R CMD check` on the built
package caught it. Separately, a merge kept pre-merge `man/` files, and `R CMD check` reported the
codoc mismatch and a failing example that `devtools::test()` cannot see. **For parallel stages the
built-package gate is not optional.**

**Ten pushbacks, ten correct.** Across this session the delegated agents stopped rather than guessed
ten times: a fifth species handler I had missed, a `gs_db` provenance premise of mine that was false
in three ways, a `stat_type` vocabulary that would have had to be invented, a golden-case
registration that my own read-only rule forbade, and six earlier. Not one was wrong. The cost was
one round each; the alternative in every case was a plausible answer with no basis.
