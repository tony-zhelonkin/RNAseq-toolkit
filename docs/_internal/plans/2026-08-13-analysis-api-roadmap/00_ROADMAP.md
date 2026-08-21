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

**The package needed a stability contract.** Nothing said which names were safe to build on.
Meanwhile 24 exports were frozen by an explicit promise, the `coresh_*` functions were new and
unproven against a consumer, and inherited live exports sat outside every layer prefix. Adding
`tf_activity()` on top of that would mean a consumer could not tell a settled name from a
provisional one, which was the failure the refactor was meant to end.

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
| Exports | **78** — 57 live, **21 deprecated shims** |
| Live exports under an analysis-layer prefix | 43 (`gs_` 18, `gsdb_` 6, `de_` 6, `gatom_` 6, `coresh_` 7) |
| Other live exports | **14** — 12 retained inherited names plus `bulkirna_api()` and `bulkirna_stochastic()` |
| Tests | Fresh run required after the name-audit review fixes |
| Golden baselines | Fresh run required |
| **`R CMD check`** | Fresh run required |
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
| **S1** | ✅ **done 2026-08-13.** `bulkirna_api()` ships the contract as a table. **Write the stability contract.** Three tiers: `stable`, `experimental`, and `deprecated`, with lifecycle separate from the 24-name signature freeze. Publish it as a vignette and a machine-readable table, and make `bulkirna_api()` return it. | A test asserts every export carries exactly one tier, so a new export cannot be added without classifying it |
| **S2** | **Name audit on the 13 unprefixed exports.** Decide per name: keep as a deliberate top-level verb, or move under a prefix with a shim. `download_gatom_references` → `gatom_download_refs` is the clear one, because it already belongs to a family whose other five members are prefixed. `ensure_dir` is a utility that probably should not be exported at all. | Every decision recorded with a reason in the contract vignette; no silent renames |
| **S3** | **Argument-name consistency audit.** One spelling per concept across every live export: `species`, `db`, `contrast`, `seed`, `quiet`, `verbose`, `path`, `min_size`/`max_size`. Today `gs_test()` and `coresh_search()` do not agree on everything, and the species aliases accepted by `gatom_refs()`, `gsdb_msigdb()` and `gene_to_entrez()` were each written separately. | A test enumerates formals across exports and admits only the curated vocabulary; one shared `.species()` resolver |
| **S4** | ✅ **done with S1** — `removed_in` is `1.0.0` for every shim, stated in `MIGRATION.md`. **Set the deprecation clock.** The shims exist because consumer call sites could not migrate at once. Two of three consumers are still unmigrated, so they stay with a named removal version (`v1.0.0`) and a warning that says it. | `MIGRATION.md` states the removal version; a test asserts every shim warns |
| **S5** | **Return-type and error-message consistency.** Every compute function returns a tibble; every renderer returns a ggplot; every validation error names the fix. Mostly true already — this step is the audit that proves it and the tests that keep it true. | Tests assert the class of every export's return on the shipped fixture |
| **S6** | 🟡 **implemented 2026-08-19; gate awaiting an R-capable runner.** `R CMD check` stays at zero notes, and vignettes build. One vignette per live layer: gene sets, DE, GATOM, CoReSh. **The phrase "each runnable on the shipped fixture" is unmet and superseded; see §16.** | `rcmdcheck` in the image, 0/0/0, vignettes built |

**What Phase 8 deliberately does not do:** rename anything frozen, remove any shim, or add any
function other than `bulkirna_api()`.

---

## 3. Phase 9 — one standard CoReSh and GESECA run

Small, and it is the acceptance test for Phase 8. See
[02_CORESH_METHOD.md](02_CORESH_METHOD.md) for the method and the fidelity gaps.

| Step | Work | Gate |
|---|---|---|
| **G1** | ✅ **done 2026-08-13.** Route `coresh_match(pvalues = TRUE)` to `fgsea::geseca()` — `center = FALSE`, `scale = FALSE`, `make.unique()`d rownames, stored `totalVar` kept for `pct_var`, `log2err` carried into the result. Restore p-value ordering in `coresh_search()`. Delete the guard. | p-values on real chunks; agreement with `gesecaCpp` within summed `log2err` on ≥7 of 8 datasets; both rankings reproduce the vignette's shape |
| **G2** | ✅ **done 2026-08-18.** `coresh_loadings()` and `coresh_sets()` — the rest of C3. | **The gate as written is unmet and superseded.** 28 of 58 set memberships match DC_hum_verse's GMT; the 24 that differ are explained by chunk content that no longer exists, established four ways in §10. The gate it was replaced by: the loading projection is identical to the reference on 63 of 63 comparable real hits |
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
   **This sequencing was not followed.** G2 and G4 landed while G5's web-UI half was still open,
   because that half needs a browser and is the owner's step. The mitigation was to gate G2 on
   agreement with the reference implementation rather than on internal consistency, which is
   weaker than an independent implementation and stronger than nothing. G5 remains the only
   external falsifier and it is still open.
5. Then finish **Phase 4** — STING-JR, then DC-nexus. Two unmigrated consumers are what keep
   the deprecated shims alive, so this is on the critical path to `v1.0.0`.

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
`coresh_match()` signature with no way to know it is provisional. At that point the deprecated
shims had no removal date, so they were permanent by default. Adding three subsystems to that
would make the package's own surface the next thing needing a refactor.

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
independent reviewer. The checkpoint gates were **1212 tests passing, golden 20/20, and
`R CMD check` 0/0/0**. The registry has grown since that checkpoint.

### S1 — the contract exists

`bulkirna_api()` returns one row per export with `name`, `layer`, `lifecycle`, `frozen`,
`superseded_by`, `removed_in`. The current registry has 46 stable, 11 experimental, 21
deprecated, and 24 frozen names; all shims are removed in `1.0.0`.

**The two-axis design was the right call:** 21 of the 24 signature-frozen names are now
deprecated, so a single tier column would have to lie about one axis or the other.
`frozen ∖ deprecated` is exactly `{build_dge, ensure_dir, format_pathway_name}`.

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
reviewed here. The current package surface is **78 exports: 46 stable, 11 experimental, and 21
deprecated**. The checkpoint gates recorded here predate this review-fix pass and require a
fresh run.

The delegated agents had no Docker socket, so **every gate in this section was run by me**, not by
the agent that wrote the code. Both agents were told to claim no gate as passing, and neither did.

### S2 — the 13 unprefixed exports, and only one moves

`04_NAME_AUDIT.md` records all thirteen decisions with reasons. `download_gatom_references` →
**`gatom_download_refs`**, with the frozen name kept as a deprecated shim: five other
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

---

## 12. The independent review of G2, and what it found (2026-08-18)

An Opus reviewer read the merged set-building layer. It confirmed the ported math in detail and
found **two real defects, three assertions that cannot fail, and one uncovered line that matters**.
Recording it here because most of it is about *how* the code was verified, not what it computes.

### What it confirmed, with its own measurements

- **Tie order and the `top_n` boundary are not a defect.** `order(-abs(x), method = "radix")` and
  the reference's `order(abs(x), decreasing = TRUE)` are identical *including tie order* — radix is
  stable in both directions, verified in R 4.5.3 on `c(5, 1, 5, 3, 5)`. `E %*% (profile / norm)` is
  bit-identical to the reference's two-step form.
- **NA and duplicated Entrez handling is better than the reference**, not merely equal: `query` is
  validated non-`NA` and the rownames are integer, so `match()` cannot false-match `NA` to `NA`,
  which the reference's `na.omit(match(...))` could.
- **The compendium numbers reproduce independently.** 42,465 unique accessions, 1,635 on more than
  one platform, and **zero repeated `(gse, gpl)` pairs** — so the composite key is in fact a key.
- **The `gs_db` provenance change is sound**, including the paths the tests miss: a zero-set database
  with a zero-row `set_provenance`, a base-`data.frame` `set_provenance`, and `db[0]`.

### It tried to falsify "the input data no longer exists" and could not

This is the part worth reading. §10's conclusion rested on my 63-of-63 loading agreement, and the
reviewer identified the one deviation that check **structurally cannot detect**: `coresh_loadings()`
de-duplicates the query where the reference double-weighted a repeated id in `colSums()`. Feeding
the same query object to both implementations makes the deviation cancel.

So it tested the deviation directly against the real stored queries: all 16, k = 2–25, **zero
duplicates and zero NAs**, and `sym2ent()` never de-duplicated either. Ruled out. Its second
alternative — the wrong platform for a duplicated accession — explains at most 3 of the 24 differing
sets, since only 3 of the 58 accessions are multi-platform. Query provably identical and loadings
differing means `E1024` differs.

**A verification that cannot see a known deviation is not a verification of that deviation.** I had
the right measurement and the wrong confidence in its scope.

### The defect that matters: the honest-failure hole moved rather than closed

`coresh_sets()` counted **extractions**, not **sets**. A hit that extracted cleanly and then fell
outside `min_size`/`max_size` was dropped by a bare `next` — uncounted, unmessaged, and absent from
`failures`, so the total-failure `stop()` could not fire. Reproduced before acting:

```
coresh_sets(hits, list(q = 1:5L), top_n = 5L, min_size = 15L)
→ RESULT class: gs_db  sets: 0  -- no error, no message
```

**The G2 work closed the reference's `tryCatch`/`"skip:"` hole and reopened it one step downstream**,
which is worse than leaving it, because the code now looks like it handles the case. The ways in are
ordinary: retired Entrez ids, a species mismatch reaching `entrez_to_gene()`, a `top_n` smaller than
`min_size` — which was unsatisfiable and unvalidated — or duplicate-collapse-plus-NA-drop on a real
dataset like the fixture's `na_ids` object.

### And the blind spot repeated itself one layer up

`mapped_ids <- loadings$entrez[!is.na(loadings$entrez)]` had **zero coverage**: both `coresh_sets()`
tests mocked `coresh_loadings()` *and* `entrez_to_gene()`, so no test ever passed an `NA` id to
symbol mapping. Delete the guard and the suite still passes; without it, every hit in an
NA-carrying dataset errors.

§8 recorded "a hand-built fixture contains only what its author imagined" and answered it with a
real-data fixture. **This is the same lesson at the next level: a real fixture proves nothing if the
test mocks past it.** Three assertions were also unfalsifiable — an `expect_identical(first, again)`
on a pure function, an `all(is.finite())` that holds for any implementation, and a
duplicate-and-NA test that asserted neither property in its name, comparing instead against an
inline re-derivation that would have made the same mistake.

### On reviewing

The reviewer named, for each of the tests it did *not* fault, the change that would break it — and
found two code paths working but untested, by hand. That is a more useful shape of review than a
list of suspicions: it says what the suite actually holds.

### What was fixed in response

All of it, in one round, gated here. **1,596 tests passing, golden 20/20, `R CMD check` 0/0/0.**
Size drops are counted and always reported, a zero-set outcome always says so, and
`min_size > top_n` now errors — `min_size` (15) must not exceed `top_n` (5); symbol mapping can only
retain or reduce the extracted genes. The taxonomy separates a lookup miss, an extraction failure and
a name collision. The NA guard is exercised by a test that runs the extraction for real on the
`na_ids` fixture object and stubs only the mapper. Provenance survives the human-to-mouse rebuild,
its `data.frame` branch validates column types, and `gsdb_info()`'s two shapes now share their
common columns.

Verified after the fix, both paths:

```
min_size 15 > top_n 5  → `min_size` (15) must not exceed `top_n` (5); symbol mapping can only ...
a legitimate zero-set  → size filter dropped 1 of 1 hits outside [1, 2] genes.
                       → 1 hits attempted, 0 sets produced.
```

**Four of the reviewer's readings were wrong in detail, and the agent fixing them said so** rather
than working around them: there is no `coresh_score()` in this package (the comparison is
`coresh_match()`), the finite-loading assertion was redundant rather than strictly unfalsifiable,
`coresh_loadings()` retains missing reference ids on purpose so NA removal belongs to `coresh_sets()`,
and there is no GMT round-trip to lose provenance through — the only real reconstruction loss was the
species conversion. The finding underneath each was still valid.

Its process point stands and the tables above now reflect it. §5 said run **G5 before G2–G4**, and
that was not done, because G5's remaining half needs a browser. G2's stated gate reads as passed and
was not met. Both are now recorded as unmet-and-superseded with the evidence attached, rather than
left to look green.

---

## 13. The independent review of S2/S3 (2026-08-18)

The second reviewer's verdict: **the change is sound.** All 14 field reads across the six migrated
call sites are correct, no `NA` can reach a filename, and the alias-by-alias comparison found the
consolidation to be a strict widening with one exception. It also found nine things worth fixing, all
now closed. **1,601 tests passing, golden 20/20, `R CMD check` 0/0/0.**

### The blocking finding was documentation contradicting the code it ships with

`R/api.R:16` said "all **20** deprecated shims are frozen" while the registry in that same file
returned **21**, and `R/api.R`'s roxygen ships as `man/bulkirna_api.Rd`. `MIGRATION.md` contradicted
itself twenty lines apart. Both were touched by the change; these two lines were not.

The lesson is not "update the count". It is that **a hardcoded count in prose is what created the
problem.** Counts in tests now derive from the registry, and the prose says "all deprecated shims"
where the number added nothing.

### Three things the audit should have asserted and did not

| Gap | Measured failure |
|---|---|
| `biomart_dataset` asserted **nowhere** in 1,596 tests | Swap the human and mouse marts in `R/species.R`; every test passes, golden stays 20/20, and `annotate_genes(species = "human")` silently queries the **mouse** mart |
| `test-dge.R:79` compared two spellings of one species to each other | Point human's `orgdb` at `org.Mm.eg.db`; both spellings still agree, so the test passes while every human annotation comes from the mouse database |
| Nothing checked the shim forwards `dest_dir` to `dir` | Drop the forwarding; every test passes, and `download_gatom_references(dest_dir = "/scratch/refs")` writes 16 MB into the repo-relative default — **the same failure this audit round already had to chase** |

All three are the wrong-result-not-error class. A resolver whose purpose is to map one concept to
seven representations now has all seven pinned for both species.

### Two behaviours in the resolver were wrong, both verified before fixing

**Its error message named the wrong fix.** `.gsdb_species_label(NA)` reported that `species` must be
one of six human and mouse aliases, on the one path that genuinely accepts any non-empty string. The
message it replaced was correct. `CONVENTIONS.md` §5 asks an error to name the fix; this one named a
different function's.

**And it silently promoted a genus to a species.** Partial matching ran before the custom branch, so
`.species("Mus", allow_custom = TRUE)` returned `"Mus musculus"` where the old code returned `"Mus"`.
`gsdb_register(species = "Mus")` therefore stored a claim the user did not make, and it reached
`gs_result$species`. Partial matching exists only to preserve `match.arg()`'s contract in
`annotate_genes()`, which never asks for custom species, so it now applies only when
`allow_custom = FALSE`. Verified after: `"Mus"` stays `"Mus"`, `"Homo"` still resolves to human.

`04_NAME_AUDIT.md`'s claim that an unknown label "only has underscores changed to spaces" was false
for that case, and is fixed.

### The narrowing, and it is the one we want

`annotate_genes(species = NULL)` used to return `"Mus musculus"`, because `match.arg()` returns
`choices[1L]` for `NULL`. It now errors. So `annotate_genes(ids, species = cfg$species)` with an
unset config key silently meant mouse yesterday and aborts today, which catches the mistake before it
annotates human genes through the mouse package. `NEWS.md` records it, per the freeze rule that a
corrected result is not a breaking change but is written down.

### The vocabulary test had shallow teeth, and the fix is a different shape

S3's stated goal is "one spelling per concept". The test was a **25-name blacklist**, so
`target_dir`, `sp` and `rand_seed` all passed, and it never asserted that a canonical name was
*present* anywhere. The first replacement was a flat 178-name allowlist. That closed the set of
spellings but did not model concepts: adding a spelling to the list was enough to make it policy.

The vocabulary is now grouped into named argument concepts, each with a reason. Every live formal
must belong to exactly one concept, every recorded spelling must remain in use, and the 11 canonical
cross-layer formals must remain present. A concept may have multiple spellings only when the exact
set is recorded as a reasoned exception. The four non-snake-case formals have a parallel exact
exception record. Adding a fourth spelling to the result-limit concept therefore fails even after
the author classifies it there; accepting it requires deliberately widening the exception.

This remains a curated audit, not a type system for argument names. Domain-specific inputs and the
renderer aesthetic surface are single-use concepts under shared reasons. An author can misclassify
a synonym as a new singleton concept, and no structural test can infer that semantic mistake. That
singleton bucket is where the residual enumeration lives: appending there is still the cheapest
mechanical evasion. The test makes classification the cheapest honest fix; review still has to judge
whether the classification is true.

Three candidate unifications came out of building it and were initially left for an owner decision:
`top`/`top_n`/`n_top`, `color_palette`/`colours`/`palette`, and four non-snake-case formals
(`B_cutoff`, `baseMean`, `log2FC`, `max.overlaps`). On 2026-08-19 the owner chose `top_n` and
`palette`; the four upstream non-snake-case formals remain reasoned exceptions.

The layer-coverage test was checked and left alone: it does have teeth, and `test-api.R`'s
registry-versus-`NAMESPACE` comparison closes the drift hole that would otherwise let an unregistered
export slip past it.

### A gap in what "R CMD check 0/0/0" covers

**Both source-tree enforcement tests skip under `R CMD check`**, because an installed package has no
`.R` files to parse — and that includes `test-stochastic.R`'s parse-tree test, the invariant
`AGENTS.md` advertises most confidently. The skip is correct; the silence about it was not. Both now
say why they skipped, so nobody reads a green check as covering them. **They run under
`devtools::test()` only, and that is now stated where the gates are listed.**

### Smaller things closed

The deprecation warning pointed at `help("bulkiRNA-deprecated")`, which does not exist — it was the
only one of 21 shims passing `package =`. `R/gsdb-gatom.R`, which no longer held any `gsdb_*`
function, is now `R/gatom-download.R`. `04_NAME_AUDIT.md`'s reason for keeping `bulkirna_check_deps`
was an argument for the opposite conclusion and is rewritten. And `sQuote()` in the resolver's error
became escaped quotes, matching `CONVENTIONS.md` §5.

### It also strengthened the audit's own evidence

Checked against `sciagent-rna/docs/data/used-functions.tsv`: `ensure_dir` appears in exactly the two
unmigrated consumers, so that "keep" row's premise is precisely true. **`download_gatom_references` is
absent from the inventory entirely — zero consumer call sites**, which is a better argument for the
rename than the family-consistency one the document gave. `annotate_genes`, `gene_to_entrez`,
`entrez_to_gene` and `filter_confounder_genes` are also absent, so their "keep" rows rest purely on
the conceptual argument, which stands: there is no annotation or gene-id layer in the four-layer
architecture, and a rename would mean inventing one.

### On this review, and on the last one

Both reviewers were wrong in detail and right in substance, and both times the agent fixing the
findings said which: here, that not every `annotate_genes()` test passes `use_biomart = FALSE`, that
the species record has seven fields rather than six, and that the existing `skip_if()` calls did have
reasons — merely ones too vague to expose the coverage gap. **The finding underneath each still held.**
The useful review is the one that names, for each test it does *not* fault, the change that would
break it. Both did that, and between them they found four assertions that could not fail and two code
paths with no coverage at all.

---

## 14. Measuring the concept grouping against the criteria set for it (2026-08-19)

The reviewer set three falsifiable tests for §13's grouped vocabulary before it landed, rather than
judging the structure by eye. All three were run. **One passed, one failed by a wide margin, and
the third found a defect in the exceptions themselves.**

### Test 1 — did the bidirectional check survive? Yes

The concern was that `setdiff(vocabulary, observed)` is one line over a flat vector and needs an
`unlist()` over a grouped structure, so the refactor could trade a real property for a nominal one
and rot silently. It did not: `test-name-audit.R:190-194` flattens the grouping into `audited` and
asserts **both** `setdiff(observed, audited)` and `setdiff(audited, observed)` are empty. A formal
that leaves the API must leave its concept by name.

### Test 2 — what fraction landed in the not-an-audited-concept group? 88%, and the threshold was 67%

The reviewer's own numbers: 30 of 178 would mean the grouping bought something real, 120 would mean
the enumeration had relocated. Measured:

| Formals | Reason group |
|---|---|
| **102** | "API-specific input, one spelling, does not compete with another audited concept" |
| **55** | "renderer or presentation input, one spelling, controls one distinct aesthetic" |
| 11 | the canonical cross-layer concepts |
| 6 | the two multi-spelling concepts |
| 4 | upstream non-snake-case spellings |

**157 of 178 sit in two generic buckets** — 88%, well past the 120 the reviewer named as failure. 174
concepts for 178 formals means the structure is very nearly one concept per formal.

So the honest verdict is the one the reviewer predicted: **the enumeration mostly relocated.** Two
boilerplate reason strings covering 157 names are not reasons, they are a category label. What the
grouping genuinely adds over the flat list is narrow — the multi-spelling exception check, which
guards exactly two concepts — and it costs 234 lines where the flat list cost 30.

It is kept rather than reverted because it is a strict superset of the flat list's properties and
because the exception check is what found the defect below. But **§13's claim that this is "the
difference between a blacklist and a contract" was too strong, and so was my first correction of it.**

I wrote that the cheap escape "moved from appending a name to declaring a false singleton concept".
That is still wrong, and the reviewer measured why. The multi-spelling assertion only bites if an
author volunteers the collision by adding the name *inside* an existing concept. The path an author
actually takes is the `unexpected` failure, whose own message offers the way out — *"Classify each
under an existing concept or add a new reasoned concept."* Appending `"n_features"` to the 102-name
alphabetical vector satisfies it in one line, under a sentence somebody else wrote, and the
multi-spelling assertion never sees it because the name is now its own singleton. **No false claim is
required. The cheapest move is still appending a name**, and the trace it leaves — one line added to an
alphabetical list of 102 — is the same trace the flat vector left.

The mechanism is the `single_use()` helper. `top_level` has teeth because the name is the key and the
reason is the value, so appending means writing a sentence. `single_use(formals, reason)` maps one
shared reason across ~100 names, which removes exactly that property. And the reason-quality assertion
— that every reason is a one-line non-empty string — iterates **174 concepts holding 9 distinct
strings**, so it is satisfied by construction for anything appended. That check is close to
assert-nothing in its present form.

Measured at HEAD, after the unification: **174 concepts for 174 formals, all singletons**, 9 distinct
reasons, **168 of 174 formals under a reason shared with other formals**, and 157 in the two
boilerplate buckets. So the multi-spelling assertion now has **no live subject at all** — it guards
against reintroducing a second spelling, which is worth having, but it cannot fire on the append path.

### What the grouping does buy, measured rather than asserted

Three properties are real and none existed in the flat list:

- **The snake_case check derives its offenders from the data by regex** and cross-checks a four-name
  reasoned allowlist. It cannot be appended to without writing a reason.
- **"Each audited formal belongs to exactly one concept"** is a genuine assertion.
- **The bidirectional anti-rot check survives per-concept.**

Everything else the grouping claims is aspiration. If the 234 lines ever need to shrink, the honest
reduction — the reviewer's, and I agree — is to collapse the two boilerplate buckets back to a flat
vector and keep those three: most of the value at a fraction of the cost.

### Test 3 — the fourth-spelling check passed on the path it covered

Planted `n_features` in `result_limit`:

```
`result_limit` has [n_features, n_top, top, top_n]; its exception records [n_top, top, top_n].
A concept may have multiple spellings only when that exact set is a reasoned exception.
Reject the new spelling or deliberately revise the exception.
```

That was a classification. But it only fired for a name classified into an existing multi-spelling
concept, which was the one path an author trying to get a name through would not choose.

### The defect: both multi-spelling exceptions were argued from a false premise

The reviewer asked for the three unrenamed concepts to be named, on the grounds that **a recorded
exception in an enforcement test is load-bearing in a way a table row is not** — it is the thing that
stops the test failing. That was the right instinct. Measured across every live export:

| Concept | Spellings and their public paths before the decision | Frozen? | Consumer call sites |
|---|---|---|---|
| `result_limit` | `top` (`gs_leading_edge`, `gs_plot_bar`, `gs_plot_dot`, `gs_plot_running`, plus `gs_plot_heatmap` methods), `top_n` (`coresh_convergence`, `de_bfc_plot`, `de_md_plot`, `de_volcano`), `n_top` (`coresh_loadings`, `coresh_sets`) | **none of the 11** | **0** |
| `colour_palette` | `color_palette` (4 `de_*`), `colours` (`gs_plot_bar`, `gs_plot_dot`, plus `gs_plot_heatmap` methods), `palette` (`gs_plot_running`) | **none of the 8** | **0** |
| non-snake-case | `B_cutoff`, `baseMean`, `log2FC`, `max.overlaps` | none | 0 |

Both reasons said "frozen signatures". None of the affected signatures is frozen, and none of the
exports appears in `used-functions.tsv`.

**The defect is isolated to these two exceptions, and saying so matters.** The reviewer checked the
*other* keep decisions against the same standard — the 13 rows in `04_NAME_AUDIT.md` — and they hold.
`ensure_dir`'s "two unmigrated consumers call it directly" is exactly true: 2 files, in
`DC-nexus/DC_Dictionary` and `STING-JR/mouse_anchor`, which are precisely the two unmigrated
consumers. So "the audit cited a constraint that does not exist" would be the wrong summary; the
audit's thirteen rows did not, and the enforcement test's two exceptions did.

**And that asymmetry is the general lesson, independent of these two.** A table row in a plan document
is inert: if its premise is false, the row is merely wrong. **An exception inside an enforcement test
is load-bearing** — a false premise there does active work, holding open the one gap the test exists to
close, under a sentence that reads as though the question were settled. Audit the reasons attached to
exceptions before the reasons attached to decisions. The original measurement found 17 explicit formal/export
occurrences. Source inspection before the rename found two more public argument paths:
`gs_plot_heatmap()` declares `top` and `colours` on its `gs_result` and `gs_matrix` methods, reached
through the exported generic's `...`. The same inventory check returned zero call sites for that
export too.

That finding exposed a limit in the enforcement test: walking `formals()` on exported generics does
not see method formals hidden behind `...`. The test now checks the two heatmap method signatures for
the retired spellings explicitly.

**Zero measured cost is an argument that the unification *can* be done, not that it *should*.**
`top`, `top_n` and `n_top` may each read correctly in their own context, and that is a judgement about
the owner's own API rather than a defect in it. What the measurement settled was only that nothing
prevented it; the preference was still theirs. The owner decided on 2026-08-19 to unify the result
limit on `top_n` and the colour palette on `palette`. All affected formals were renamed in place, including the heatmap methods, so positional
callers remain unaffected. No argument-level deprecation aliases were added. The two multi-spelling
exceptions were removed; only the four upstream non-snake-case exceptions remain.

One more defect surfaced at the moment the exceptions went away, and it is worth recording because of
where it lived. With no multi-spelling concepts left, the test compared an empty **named** list against
an empty **unnamed** one and went red — on exactly the state it exists to reward. Same family as the
`paste0(character(0))` defect §8 records: **the passing path of an enforcement test is itself untested
until the day it passes.** An audit that has always had exceptions has never run its own success case.

**1,606 tests passing, golden 20/20, `R CMD check` 0/0/0.** "One spelling per concept" is now a fact
about the API rather than a goal in a document.

### The counting method, not the count

The `gs_plot_heatmap()` blind spot was in the *method* both the audit and I used, and the reviewer
volunteered that it was theirs too: verifying "exactly 12 live exports expose `species`" and "exactly 4
expose `dir`" by walking `formals()` on exported names uses the same technique with the same hole.

So I checked whether it cost anything, since asserting it did not would repeat the mistake. Across all
**8 registered S3 methods**, only `gs_plot_heatmap`'s two hide formals; every other method takes `...`
plus display arguments, and `gs_test` declares `db` on the generic. **No method conceals a `species` or
`dir` formal.** The 12 and the 4 stand — by the distribution of methods in this package, not by the
soundness of the method used. Those are different claims and only the first was ever true.

Then the same check found a live gap. Of the twelve formals declared only on a method,
**`group` and `samples` were not in the audited vocabulary at all**, because `observed` walked exports
only. The audit claimed to cover the live public surface and missed two names. `observed` now includes
the formals of every S3 method registered for an exported generic, resolved through the namespace
method table rather than by stripping a dot-suffix — `dplyr_reconstruct.gs_result` and
`tbl_sum.gs_result` have dots in the generic. `group` and `samples` are classified, and the
heatmap-specific spelling check the earlier round added is **removed as redundant**: an assertion that
looks like a guard and guards nothing is worse than no assertion.

### Every allowlist had only ever run non-empty

The generalisation of §14's empty-list defect, and the reason it is worth acting on before `v1.0.0`
rather than during it. Every test comparing an observed set against a curated allowlist — `top_level`,
the non-snake-case exceptions, `frozen ∧ stable`, the message-only shim list, and **the deprecated-shim
lists, which empty at `1.0.0`** — had run only with that list non-empty. None had demonstrated it
passes when the thing it guards against is finally absent. At `1.0.0` the 21 shims go and
`bulkirna_api(lifecycle = "deprecated")` returns zero rows, so several comparisons meet zero for the
first time during a release.

Those comparisons now run through **one shape-stable helper, called with both the real allowlist and an
empty one**, so the empty path is exercised by the same code rather than by a parallel
reimplementation — which would have been a second thing to keep in sync.

**Measured rather than predicted**, both directions:

| Check | Result |
|---|---|
| Empty the deprecated tier, simulating `v1.0.0` | **57 assertions pass, 0 fail** |
| Drop one live export from its allowlist | **Fails, naming `theme_bulki`** |

The second matters as much as the first: a helper that tolerates zero could have been written to
tolerate everything.

**1,621 tests passing, golden 20/20, `R CMD check` 0/0/0.**

### The two findings interlock, and the fix is now load-bearing

The reviewer's closing observation, which is a better statement of the lesson than either of us had:

Since the unification there are 174 concepts for 174 formals, **all singletons**, so
`multiple_spelling_exceptions` is literally `list()` and the multi-spelling assertion compares two
empty structures. That is exactly the state that turned the suite red — and it is no longer a transient
condition someone hit while emptying a list. **It is the permanent steady state of the test.**

So the shape-normalisation is **load-bearing rather than incidental**: it is the only reason that
assertion passes today, and it stays that way unless someone reintroduces a second spelling, which is
the thing the assertion exists to prevent. An assertion whose success path is the empty case, guarding
a condition that is currently absent, reads to the next person as a comparison of two empty structures
— which is to say, as the most deletable line in the file.

**The general form:** it is not only that a test with a permanent exception list has never run its
success case. It is that **when the exception list finally empties, the empty case stops being the
untested path and becomes the only path** — so the fix that made it work *is* the whole assertion, and
it looks like dead code.

That is a comment-in-the-code problem, not a documentation problem, so both sites now say so at the
point of temptation: `named_list()` in `test-name-audit.R` and the two allowlist helpers in
`test-api.R` carry "do not delete as dead code", with the reason and the measurement. The empty
`multiple_spelling_exceptions` says it is the goal state rather than an unfinished list.

And it sharpens where the pre-`1.0.0` check belongs. The deprecated-shim lists empty at `1.0.0` **and
stay empty**, so the same inversion happens to them: `test-api.R`'s message-only list, the
`frozen ∧ stable` set and the deprecated-metadata assertions all become permanent empty-versus-empty
comparisons. Doing that check now, while a non-empty case still exists to compare against, is the only
time it can be done properly.

### Two notes on how the empty-case work was checked

**An invalid probe that nearly became a finding.** My first attempt at measuring the empty case emptied
the non-snake-case exception list *while the four offenders still existed*, and the run halted. Reading
that halt as a finding would have been a false positive — that state **should** fail, since the
offenders are still there and no longer excused — and the edit had also broken the file's syntax, so
the halt was not even the failure it appeared to be. It was caught by asking what the valid empty state
actually is: offenders and exceptions both absent.

That is `AGENTS.md`'s "an empty result and a broken run must be distinguishable" from the other
direction. Here a broken *probe* and a real defect were indistinguishable from the exit status alone,
and the same discipline applies to the measurement as to the code being measured.

**And an empty-case check needs a positive control.** "57 assertions pass, 0 fail" after emptying the
deprecated tier cannot on its own distinguish a helper that correctly tolerates zero from one that
tolerates everything. It is the paired negative — dropping `theme_bulki` from its allowlist and
watching the same helper fail, naming it — that makes the first number mean anything. Neither run is
worth much alone.

### What the two reviews cost, and what only measurement found

Four rounds of fixes across two reviews. **Four findings surfaced only because a claim was executed
rather than read**, and the split is worth keeping:

| Finding | Surfaced by |
|---|---|
| `biomart_dataset` asserted nowhere; swapping the marts passes 1,596 tests | reading the change adversarially |
| Two exceptions citing frozen signatures across 17 non-frozen exports | reading the change adversarially |
| 157 of 178 formals in two boilerplate buckets, against a 120 threshold | counting instead of reading the structure |
| `group` and `samples` unaudited, behind a generic's `...` | executing a claim that had only been asserted |

The first two came from the reviewer, the third from taking its threshold literally, the fourth from
checking the reviewer. **The ones that took longest to surface were the ones where a sentence sounded
true** — which is the argument for the whole exercise, and for measuring one's own prose as well as
one's code.

---

## 15. G3 — `gsdb_coresh()` (2026-08-19)

**Phase 9 is now complete except G5's browser half.** `gsdb_coresh()` composes `coresh_search()` and
`coresh_sets()` into one provider returning a `gs_db` with `database = "coresh"`, so a CoReSh-derived
database is obtained and consumed exactly as an MSigDB one is. **1,680 tests passing, golden 20/20,
`R CMD check` 0/0/0, 79 exports.**

### The gate, run end to end on the live compendium

C4's original gate — byte agreement with DC-nexus's stored GMT — is void for the reasons in §10, so the
replacement was mostly compositional, since G2 already pinned the numerics against the reference
implementation on 63 of 63 hits.

Three real queries against `syn66227307_20260721`, and the check that matters:

```
Gate 1: provider and hand composition identical, including provenance: TRUE
<gs_db> CoReSh-derived gene sets [coresh]  (Homo sapiens, symbol)
13 sets, 561 unique genes, size 50-50 (median 50)
Provenance: source=coresh; snapshot=syn66227307_20260721;
            chunk_dir=/refcache/coresh/current/preprocessed_chunks/hsa;
            species=Homo sapiens; n_chunks=89; queries=iron_uptake(5),
            ferroptosis(4), heme(4); top_hits=5; top_n=50; min_size=15;
            max_size=500; jaccard_threshold=0.8
```

`gs_test()` then returned a valid `gs_result` on it with no renderer change. **The provider adds no
computation of its own**, which is what "it enters at the provider layer and bends nothing" has to mean
in practice rather than in prose.

**ADR-004's second real exercise, and it lands where the ADR asked**: the snapshot tag is in the
object's provenance, not the log.

### The defect, and it was my instruction that caused it

The first version recorded the queries **and the entire selected-hits table as `dput()` text** —
several kilobytes of escaped R code in every `print()`, carried through every subset.

The agent's reasoning was sound and is in the code: database provenance is constrained to atomic
scalars, so serializing was the only way to fit a list and a data frame into it. **The error was in my
brief**, which said a `gs_db` whose provenance names a snapshot but not the query set cannot be
reproduced from its own record. Taken together with the scalar constraint, that instruction forces
exactly the blob that was written. **It satisfied the validator by defeating the purpose of the field it
was validating.**

So the fix was to correct the instruction. Provenance now records query **names with their sizes**
(`iron_uptake(5), ferroptosis(4), heme(4)`), the hits table is gone because `set_provenance` already
ties each set to its `query_name`, `gse`, `gpl`, `loading_cutoff` and `rank_in_coresh`, and both
functions document plainly that **the exact Entrez ids stay with the caller**. I declined a digest of
the query ids: base R has no object hash and `digest` is not a dependency, so it would mean a new
`Suggests` for a provenance nicety.

**An honest gap a reader can see beats a blob that looks like completeness.** That is the same principle
as `.ref_path()` warning when `current` is not a symlink rather than recording the literal string.

### What Phase 9 looks like finished

| Step | State |
|---|---|
| G1 p-value path | ✅ |
| G2 set-building layer | ✅ gate superseded, replacement met |
| G3 provider | ✅ |
| G4 `gs_coregulation()` | ✅ golden gate deliberately not met, stronger check substituted |
| G5 sweep against GEO titles | ✅ 7 of top 10 hypoxia or HIF |
| **G5 web UI** | **🚫 needs a browser — the owner's step, and the only remaining external falsifier** |

---

## 16. S6 — four layer vignettes, with the runnable gate corrected (2026-08-19)

Four vignettes and their build plumbing are implemented. `DESCRIPTION` now declares
`VignetteBuilder: knitr`, with `knitr` and `rmarkdown` in `Suggests`. The current development image
contains knitr 1.51 and rmarkdown 2.31, so a future build failure is not expected to be a missing
builder. **No R or Docker was available to the implementing agent, so no vignette has been shown to
build and S6's `R CMD check` gate remains open.**

### The phrase "each runnable on the shipped fixture" is unmet and superseded

Both `tests/fixtures/` and `data/` are excluded by `.Rbuildignore`. A built vignette can reach neither
the count fixture nor the source reference tree. Moving either into the build would violate the
fixture ownership rule, and shipping a new counts payload solely for documentation would permanently
grow the package without an analysis need. The gate as written is therefore impossible.

The replacement is stronger where the package has data and honest where it does not:

| Layer | Build-time state | Data and reason |
|---|---|---|
| Gene sets | **evaluated** | `gsdb_load()` reads real MitoPathways sets from shipped `inst/extdata`; a synthetic rank vector carries a deliberately planted signal in one of those sets. `gsdb_from_file()` and `gsdb_register()` demonstrate the other provider boundaries without network access. |
| DE | **evaluated** | The vignette constructs a labelled synthetic count matrix with 20 four-fold-up and 20 four-fold-down genes, then runs `build_dge()`, limma-trend, `gs_ranks()`, `de_pca()` and `de_volcano()`. No fixture or external file is needed. |
| GATOM | `eval = FALSE` for analysis | The optional packages plus about 24 MB of network reference files are prerequisites. The vignette gives the real download, resolution, DE validation, seeded solve, edge-gene extraction and HTML-render calls without inventing output. |
| CoReSh | `eval = FALSE` for analysis | The search needs the roughly 20 GB chunk tree, `qs2`, organism annotation and, for parallel work, BiocParallel. The current image lacks `qs2`. The vignette gives preflight, search, `(gse, gpl)` convergence, provider construction, enrichment and render calls without pretending the compendium ships. |

`inst/extdata/` is the one genuinely runnable shipped data layer. It is not build-ignored and contains
the processed MitoPathways, mitoXplorer, unified mitochondrial and TransportDB databases plus the
master-table schema and metadata registry.

### Why the two unevaluated guides are still vignettes

GATOM and CoReSh belong in `browseVignettes("bulkiRNA")` even though their analysis chunks cannot run
at package build time. Their documentation value is the order of operations and the boundary between
reference preparation, compute and render; moving them to internal `docs/` would hide that workflow
from installed-package users. They contain copy-pasteable calls, state every prerequisite before the
first dependent call, and show no fabricated output. That makes their limitation visible rather than
turning an unavailable external input into a false package-build claim.

### Expected check risks

The vignette boilerplate includes index entries, engine declarations and UTF-8 declarations. The
evaluated work is intentionally small and offline. The remaining unknowns are execution unknowns:
whether every evaluated call behaves in the built namespace, whether optional-package availability
is identical under the check library, and whether rendered figures or generated vignette artifacts
trigger a check note. Only `devtools::document()`, `devtools::test()`, the golden verifier and
`rcmdcheck::rcmdcheck()` in the image can close those questions.

---

## 16. S5 and S6 — Phase 8 complete (2026-08-19)

**1,697 tests passing, golden 20/20, `R CMD check` 0/0/0 including the re-building of vignette
outputs.** Both steps ran in parallel and each found something the step was not looking for.

### S5 — the audit says the roadmap's wording is wrong, not the code

Measured across all 58 live exports rather than asserted: 17 tabular results, 9 renderers returning
ggplots, 6 writers returning their paths invisibly, and 26 in reasoned exception categories. The full
tables are in `05_RETURN_AUDIT.md`.

**The roadmap's "every compute function returns a tibble" does not hold, and the honest output is the
exception list rather than a diff.** `gs_db` is a named list, `gs_score()` returns a numeric matrix,
`build_dge()` returns a `DGEList` because edgeR consumers require one, `de_pca_3d()` returns a `plotly`
object on purpose, and `theme_bulki()` returns a theme rather than a plot. Those are architectural
contracts; changing them to satisfy a sentence in a plan would be the tail wagging the dog.

Error messages now name the offending argument and the fix, with `call. = FALSE` throughout, and two
enforcement tests keep both properties: every live export's return class asserted against its category,
and every `stop()` in `R/` parsed for `call. = FALSE` and a backticked argument.

**One self-inflicted defect, and where it surfaced is the point.** Rewriting `de_pca()`'s zero-variance
message broke a test matching on the old wording. `devtools::test()` did not show it — only `R CMD check`
on the built package did, where the failure also left a testthat `_problems/` directory that the
working-directory guard then caught, so one message rewrite produced two red assertions in different
files. **Error-text churn breaks matchers**, which my brief warned about and the round did anyway, and
the built-package gate is again the one that sees it.

### S6 — the stated gate was impossible, and the agent stopped rather than reinterpreting it

S6 asked for vignettes "runnable on the shipped fixture". **`tests/fixtures/` and `data/` are both
`.Rbuildignore`d**, so a built vignette can read neither. The agent stopped and reported that before
writing four of them, which is the twelfth correct pushback of this work.

Recorded as **unmet-and-superseded**, not reinterpreted. What replaces it, per layer, decided on what can
actually run with no network and no refcache:

| Vignette | Evaluated? | Why |
|---|---|---|
| Gene sets | **Yes** | `inst/extdata` ships real mitocarta/mitoxplorer/transportdb sets and is not build-ignored |
| Differential expression | **Yes** | self-contained synthetic counts, said to be synthetic in prose, with real signal planted |
| GATOM | No | needs ~24 MB of reference files |
| CoReSh | No | needs the 20 GB chunk tree, and `qs2` is not in the image |

The unevaluated two show real copy-pasteable calls and say what the reader must have first. **A vignette
that prints invented output would be worse than one that admits its prerequisites.**

`knitr` and `rmarkdown` join `Suggests` as **dev-only, beside `testthat`**, rather than the
optional-dependency registry: `bulkirna_check_deps()` answers "what must I install to use this feature",
and that answer never includes knitr. The registry-versus-`Suggests` drift test caught the omission
immediately — ADR-003's seam working exactly as it was built to.

### Phase 8 as a whole

| Step | Outcome |
|---|---|
| S1 stability contract | ✅ `bulkirna_api()`, two axes because 20 of 24 frozen names are also deprecated |
| S2 name audit | ✅ 13 decisions recorded, one rename |
| S3 argument audit | ✅ five species handlers unified; two concepts later unified on `top_n` and `palette` |
| S4 deprecation clock | ✅ all 21 shims removed in `1.0.0` |
| S5 return and error audit | ✅ exceptions reasoned rather than forced |
| S6 vignettes | ✅ gate superseded, two layers genuinely runnable |

**A name now tells you what you are allowed to rely on, and a test fails if that stops being true.**
Which was the whole point of the phase, and is worth stating plainly because none of it shipped a
feature.

---

## §17. Phase 4 resumed — Meta-Aging migrated, DC-nexus scoped

The recipe, the call mapping and the measured pitfalls live in
[06_CONSUMER_MIGRATION.md](06_CONSUMER_MIGRATION.md). This section records what the step decided and
what it cost, not how to repeat it.

**The owner narrowed the scope.** Phase 4's three named consumers were `14839-DM-cGAS`, STING-JR and
DC-nexus. STING-JR was surveyed and then **excluded by decision**; Meta-Aging, which was never on the
list, replaced it. The survey is kept because STING-JR is the largest consumer by an order of
magnitude — 19 files and 94 call sites across four sub-repos — and because it is the TF reference
project, so the argument for migrating it first was that its conventions would land before the
activity layer is designed. That argument is now moot and Phases 10–11 stay parked.

**The consequence for `v1.0.0` is real and needs a decision.** The 21 shims exist because consumers
call the old names. With STING-JR excluded, finishing DC-nexus does not empty the set. So `v1.0.0` is
now either "migrate STING-JR after all" or "ship with the shims present" — and the second reopens S4,
which recorded all 21 as removed in `1.0.0`.

**Scope measured, not estimated.** A bare `grep` for the 21 names across STING-JR returned 1,732
matches; the true figure is 94. The difference is the vendored `RNAseq-toolkit` submodule, which holds
the **definitions** and is checked out six times in that project, plus `.ref/` copies of other
projects and `.slice/` frozen pipelines. **Overcounting by 18× would have made this phase look
intractable.** The exclusion list is in §0 of the migration document.

**A fourth version authority, which ADR-001 did not cover.** The step resumed against `bulkiRNA 0.3.1`
in the shared library — before `gs_to_master()` existed. The development tree was at 0.5.0.9000 with
every gate green, so nothing in the repository pointed at it. ADR-001 says the package version is the
unit of reproducibility; it did not say **which installed copy** a consumer loads. Checking that is
now the first precondition of a migration.

**Meta-Aging: what the verification bought.** `14616-DM/02_analysis/scripts/06_gsea_master.R`, commit
`d6633dd` on `migrate/bulkirna`. Running it found four things that reading it had not:

1. `gs_to_master(db = )` takes the gene-set object, not a label. Passing the label inside the script's
   own `tryCatch(..., NULL)` produced **zero rows and no error**. The dominant pain-point shape, for
   the third time, now inside a consumer's error handling rather than the library's.
2. The `database` column would have silently become `msigdb_C3_TFT_GTRD` instead of the project's
   `CollecTRI_TF`, breaking every downstream join while leaving the column present and populated.
3. `gs_validate_master()` **refused** the project's `neg_log_padj` cap at 16 — *"the retired
   cap-at-16 convention is not allowed."* The attempt to preserve the old numbers failed the gate.
   ADR-002 holding against a live consumer, which is a stronger test than any in the suite.
4. The config's own `msigdbr` TODO was wrong in both directions: under 26.1.0, `CP:KEGG` must become
   `CP:KEGG_LEGACY`, while `CP:WIKIPATHWAYS` did not move, contrary to the note predicting both would.
   `gsdb_msigdb()` raised `Unknown subcollection` rather than returning an empty set, so this surfaced
   loudly.

**The fidelity measurement, and why it needed splitting.** Compared against the pre-migration
`master_unified.csv`, raw NES differences reached 4.7 — alarming until the rows were split by whether
`set_size` had changed. **On MitoPathways and MitoXplorer, version-pinned inside the package and so
carrying no reference drift, NES and padj reproduce exactly: max absolute difference 0 across 272
rows, `genes_full_set` identical 272 of 272.** Median absolute NES difference across all
unchanged-`set_size` rows is 9e-4. The rest tracks MSigDB content drift, where only 11% of GO_BP sets
kept their size. A single global tolerance would have hidden the defects and the drift in one number.

**What Phase 4 exposed about the surface**, both now `v1.0.0` blockers in §4 of the handoff:

- `convert_human_to_mouse`'s registered successor covers MSigDB only, and the real call site is a
  custom GMX. **The deprecation registry is tested for presence and shape, never for whether the
  successor can do the job.**
- No exported accessor returns the master-table columns, so a consumer reads
  `inst/extdata/master-schema-v1.csv` directly — and that record's row order is not the column order
  `gs_validate_master()` requires.

**Two limits on what was verified, stated rather than glossed.** The full job is 1,312 GSEA runs, so
the run was two contrasts against all eight collections. And `OmnipathR` is absent from both images,
so the CollecTRI and PROGENy blocks could not execute at all; they were removed from a scratch copy and
left untouched in the committed script. Those blocks stay hand-written because the activity layer is
parked, which is the recurring cost of that decision: every consumer pays it separately.

---

## §18. The version-identity failure, and what an external review made of it

Found while preparing the image pin for Phase 4. Recorded here because the conclusion was *not* to
refactor, and a decision not to act needs its reasoning kept as much as a decision to act.

**The measurement.** `DESCRIPTION` read `Version: 0.5.0`. So did the `v0.5.0` tag. `HEAD` stood 50
commits and ~3,032 changed lines of `R/` beyond it, adding nine files — `api.R`, `rng.R`, `species.R`,
`coresh.R`, `coresh-sets.R`, `gs-coregulation.R`, `gsdb-coresh.R`, `gene-ids.R`, `gatom-download.R` —
and growing the surface from 64 to 79 exports. Tags `v0.4.0` and `v0.5.0` both carry 64 exports and
neither has `bulkirna_api`.

**How it surfaced, which matters.** I proposed bumping the image pin `v0.4.0 → v0.5.0` as an upgrade
and only caught it by counting exports at each ref. Nothing in the repository objected: all gates were
green, because every gate tests content and none tests identity.

**The consultation.** `gpt-5.6-sol`, given the repository read-only, was asked whether this indicates
an architectural defect and whether the design should be reorganised. Its verdict, adopted:

> a release-identity gate failure, not a defect in the package's computational architecture … That is
> release architecture, not package restructuring.

It advised **against** reorganising `R/`, against moving reproducibility ownership back to the image,
against re-cutting `v0.5.0`, and against making an export count the identity check — *"materially
different implementations can have identical exports."* It also corrected the framing this section
started from: `v0.4.0 → v0.5.0` is a no-op **in export surface only**, since that release does carry
behaviour changes including a formatter and GATOM.

**The pattern it named**, which is the durable part: the project *"is excellent at testing object
content, but weaker at testing identity and deployment transitions"*, because those live between
repositories, tags, builds and installed libraries, where package tests do not reach. Every content
invariant here fails a test when violated; every identity invariant was prose. ADR-001 is amended with
the missing inverse obligation.

**Four more instances of the same asymmetry**, three of which were new:

| Instance | Status |
|---|---|
| `write_session_provenance()` records the version but not the source commit, so it faithfully writes an ambiguous `0.5.0` | ADR-001 amendment 4 |
| `install_core.R` records install failures and continues, so the report is observational rather than a gate | open |
| ADR-004 calls the GATOM ad-hoc download tier retired while `gatom_download_refs()` still implements it | already recorded as debt |
| "Compute never plots" is a convention, not a dependency-direction test | open, not implicated here |

**The one that explains the vanished packages.** The image documentation claims a committed
`renv.lock` makes subsequent builds deterministic. Verified, and it is false twice over:
`install_renv_project.R` restores from `/opt/settings/renv.lock` when present, the Dockerfile copies
three R scripts into that directory and **no lockfile**, so every build takes the unlocked branch;
and the tracked `renv.lock` contains **exactly one package, `renv` itself**. CRAN is pinned to the
RSPM snapshot `2026-04-15`, but Bioconductor, both r-universe remotes and every GitHub install float.
**That is the mechanism by which seven transitive packages disappeared between v0.5.10 and v0.5.13
with no commit recording the loss.** ADR-001 had noted the missing `COPY` in passing, as support for
rejecting a different option.

**The agreed sequence**, replacing "bump the pin":

1. Add the version-identity gate; it fails on contact.
2. Set the development version to `0.6.0.9000` — stricter than the usual `x.y.z.9000`-after-`x.y.z`
   convention, and deliberately so: it names the next release and makes a stale version detectable.
3. Run every gate, then cut `v0.6.0`. The 64 → 79 expansion warrants a minor release, not `0.5.1`.
4. Pin the image to the **full commit SHA** that `v0.6.0` resolves to, presented as `0.6.0`.
5. Include the explicit declarations for `OmnipathR` and the seven recovered packages.
6. Rebuild **once**, with a final verification that fails on a wrong version, a wrong `RemoteSha`, or
   a missing required package.
7. Move consumers to that image deliberately.

Do not spend the ~150 minutes to change `v0.4.0` to `v0.5.0`: it installs the genuine old 64-export
release and preserves the misunderstanding. `v0.5.13` remains usable meanwhile, provided nobody
represents it as carrying the current tree.
