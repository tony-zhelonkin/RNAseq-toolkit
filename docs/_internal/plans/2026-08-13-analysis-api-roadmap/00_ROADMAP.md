# bulkiRNA: the analysis-API roadmap, Phases 8–11

**Date:** 2026-08-13 · **Status:** plan of record for the work after Phase 7
**Parent:** [../2026-08-10-bulkirna-package/00_INDEX.md](../2026-08-10-bulkirna-package/00_INDEX.md)
**Companion:** [01_REFERENCE_PROJECTS.md](01_REFERENCE_PROJECTS.md) — the reference
implementations, their conventions and their defects. Read it before designing anything here.

---

## 0. What this roadmap is for

Phases 0–7 answered "stop the scripts from drifting". This roadmap answers a different
question the owner has now stated plainly:

> *a thin clean API that just does the job the way I like it, without having to re-state to
> agents implementing analysis for me to go search for my previous experience reference
> repos — a zero-token architecture for my personal preferences on how bulk and pseudo-bulk
> should be done.*

That reframes the package's purpose. It is not only a library. It is **the executable record
of a set of methodological preferences**, so that an agent asked for a TF-activity analysis
calls `tf_activity()` and inherits every choice already made, instead of reading four
projects and guessing which one was right.

Two consequences follow, and they set the acceptance bar for everything below.

1. **A function that needs a paragraph of instructions to use correctly has failed.** The
   defaults are the preferences. If the right `minsize` is 5, then 5 is the default and the
   documentation says why, once.
2. **The scope test tightens.** Phase 0's rule was "a table in and a table out is a package
   function; a result in and a decision out stays a skill". That still holds. The addition
   is: **a choice that has been made the same way three times is a default, not a decision.**
   Three projects agreeing on `top = 500` for PROGENy is not a judgement call any more.

---

## 1. The pain points this roadmap addresses

The first four are the original ones, restated because they have not gone away; the last
three are new and come from the survey in the companion document.

1. **Silent wrongness that a green run cannot distinguish from a null result.** Eight
   confirmed instances so far. The recurring shape: an empty result and a broken run take the
   same code path. Every new stage must be able to tell them apart, and must fail loudly.
2. **Machinery retyped per project, then drifted.** ~81% of the surveyed CoReSh lines were
   generic. The TF layer is now measurably worse: see §3.
3. **Reference data fetched live, so a result depends on the day.** msigdbr ortholog cache,
   the GATOM downloads, and now a documented OmniPath breakage.
4. **No version to cite.** Answered by ADR-001; every new stage inherits it.
5. **The same question answered by two estimators in two projects, with no record of which
   was preferred or why.** `run_ulm` in one project, `run_mlm` in another, both for
   CollecTRI TF activity.
6. **Methodological rigour trapped in one project.** The curveball-rewiring and exact
   within-donor permutation nulls exist once, in `19_regulon_nulls.R`, and nothing else can
   reach them.
7. **Parameters justified in a comment.** WGCNA's relaxed R² of 0.70 and `minModuleSize` 25
   are correct for n = 12 and are recorded in a comment block that no other project reads.

---

## 2. The ADRs, and what each one now has to carry

The four existing ADRs were written for the GSEA layer. Each extends to the new stages, and
three of them get harder.

**ADR-001 — the package version is the unit of reproducibility.** *Premise:* four version
authorities floated free and no result could name its code. *Decisive argument:* under
image-as-unit the msigdbr bug would have been invisible — tag unchanged, numbers changed.
*What changes:* nothing structurally, but the release cadence matters more once WGCNA and TF
activity land, because both are stochastic or parameter-sensitive and a moved default moves
published numbers. *Rejected alternatives:* image-as-unit, which hides dependency drift; a
lock file inside the package, which a library cannot impose on its consumers.

**ADR-002 — the package owns the master-table schema and validator, not the file.**
*Premise:* whoever computes a column owns its definition, or the definition lives somewhere
untestable. *What changes:* TF activities, PROGENy activities and WGCNA modules are all new
entity classes. **They use the existing 14-column schema with a new `entity_type`, not new
schemas.** `master_tf_activities.csv` already exists in STING-JR with 13 columns; converging
it on the package schema is part of Phase 9. *Rejected alternative:* a schema per entity
class, which is how three master-append implementations happened the first time.

**ADR-003 — `Suggests` stays.** *Premise:* it is the machine-readable list that
`R CMD check`, `install_deps` and the image build read, and 14 light `Imports` are why
`install_github` works on bare R in seconds. *What changes, and this is the real cost of this
roadmap:* `Suggests` grows by `decoupleR`, `OmnipathR`, `WGCNA` and its dependency tree
(`impute`, `preprocessCore`, `dynamicTreeCut`, `fastcluster`). WGCNA is heavy. The preflight
(`bulkirna_check_deps()`) already has a `"coresh"` feature and gains `"activity"` and
`"network"` ones. *Rejected alternative revisited:* splitting into core plus feature packages
is cleanest and is still rejected — one maintainer, one release train — but WGCNA is the first
dependency that makes the case arguable, and if a second one like it appears the decision
should be reopened.

**ADR-004 — two reference-data tiers, bundled or refcache.** *Premise:* three mechanisms had
grown and a fourth environment variable was about to appear. *What changes:* this ADR is
about to earn its keep for the second time. CollecTRI, PROGENy and DoRothEA are exactly the
case it describes — fetchable in principle, unfetchable in practice at current package
versions — so they become a **refcache source with a recorded snapshot**, resolved by
`.ref_path("omnipath", ...)`. *The insight that made the ADR work still applies:* nobody
needs to own pinning, someone needs to own **recording**. Freshness stays the default; the
drift stops being silent.

**A fifth ADR is needed.** *Proposed ADR-005 — the estimator is a default, not a parameter of
last resort.* When two published estimators answer one question, the package picks one as the
default, documents the reason, exposes the other by name, and provides the comparison as a
first-class function rather than as a script somebody writes again. Write it in Phase 9, when
the ULM-versus-MLM decision is actually made.

---

## 3. The GESECA and CoReSh problems, and how to solve them

This section supersedes §5 and §10 of the CoReSh plan and §11's open decision.

### 3.1 What went wrong in my reasoning

Three claims were made in sequence, and the middle one was wrong.

1. *"Reproduce the p-value through public `fgsea::geseca()` and prove equivalence"* — the
   plan's first-choice option. Correct instinct.
2. *"`geseca()` does not honour `sampleSize`, so the public route fails"* — **wrong, and it
   was my error.** I tested it on a null matrix only, where the p-value comes from the
   pre-permutation screen and the multilevel estimator never escalates, so `sampleSize`
   cannot show an effect. On data with real signal it plainly does: `log2err` falls
   **3.04 → 1.41 → 0.625** for `sampleSize` **21 → 101 → 501**. A one-line measurement would
   have caught this and I did not make it before writing the conclusion down.
3. *"So the only route is `fgsea:::gesecaCpp`, and adopting it is an owner decision"* — the
   conclusion that followed from the wrong premise. **It is withdrawn. No internal call is
   needed and no such decision is owed.**

The owner's pointer is what forced the recheck: GESECA is not a CoReSh add-on, it is a
first-class fgsea method with `R/geseca-multilevel.R`, `R/geseca-simple.R`,
`R/geseca-utils.R`, `R/geseca-plot.R`, `src/geseca.cpp` and a tutorial vignette.
`fgsea:::gesecaCpp` is one internal helper inside a fully public method, and the vendored
CoReSh kernel reached past the front door for no reason beyond convenience.

### 3.2 What was separately true, and stays true

- **`alserglab/coresh` v0.1.0 exports nothing.** No `R/` directory: DESCRIPTION, NAMESPACE
  (`exportPattern`), README, two licences, one vignette, three commits from April 2025.
  `coreshMatch()` and `queryGSE()` cannot be called from it. This is unaffected by the
  geseca correction, and it means the skill's documented R path and
  `validate_coresh_install.R` are both broken — a C5 defect, not a package one.
- **`pctVar` should use the stored `obj$totalVar`.** It was computed before the 1024-step
  quantization and is the more faithful value; `geseca()` recomputes `sum(E^2)` from the
  quantized matrix, which differs by ≤1.03e-5 relative. Verified on real chunks. So the
  screen stays package arithmetic and needs no fgsea call at all.
- **`gesecaSimple()` cannot serve the p-value.** It is a plain `nperm` permutation test whose
  p floors near `1/nperm`, and CoReSh reports values below 1e-28. It is still useful as an
  independent cross-check at moderate p.
- **Real chunks have duplicate Entrez rownames**, which `checkGesecaArgs()` rejects. Index by
  `match()`; call `make.unique()` before handing a matrix to `geseca()`. It changes no value.

### 3.3 The fix

**`coresh_match(pvalues = TRUE)` routes to `fgsea::geseca()`**, with
`center = FALSE, scale = FALSE`, `make.unique()`d rownames, `sampleSize` exposed and
defaulting to 21 to match the published CoReSh behaviour, and `log2err` carried into the
result as a column. The `stop()` currently guarding that path is removed, `p_value` stops
being always-`NA`, and the message that named three dead ends is deleted along with the
premise that produced it.

**Agreement is checked within the error bound, not by equality.** Measured on real chunks at
`sampleSize = 21`, `fgsea:::gesecaCpp` and `fgsea::geseca()` agreed within their summed
`log2err` in **7 of 8 datasets**, with `|log2(ratio)|` up to 2.55 against a median summed
bound of 1.53. That is what two Monte-Carlo estimators of one quantity look like. The
regression test asserts the bound, and `sampleSize` is raised where precision matters.

**No `:::` enters the package.** ADR-003's `R CMD check` cleanliness is preserved, an fgsea
release cannot silently break us, and §5 option 2 is finally dead rather than deferred.

### 3.4 The opportunity this opens

GESECA needs an expression matrix and no contrast, which makes it complementary to
DE-driven GSEA rather than an alternative. It belongs in the compute layer beside `gs_test()`
and `gs_score()`:

```r
gs_coregulation(E, db, sample_size = 21L, center = TRUE, min_size = 10L,
                max_size = 500L)                                  -> gs_result
```

It returns a `gs_result` with `stat_type = "pct_var"`, so every existing `gs_plot_*`
renderer works on it unchanged — which is the four-layer architecture paying for itself a
third time. `plotGesecaTable()` is fgsea's own renderer and is **not** what we ship; the
owner's visual style goes through `gs_plot_*` and `gs_theme()` like everything else.

The remaining question, deliberately left open: **does `gs_coregulation()` become a `method`
of `gs_test()` or its own verb?** It takes a matrix where `gs_test()` takes ranks, so a
shared signature would lie about the input. Recommendation: its own verb. Decide before
writing it.

---

## 4. The phases

Phases 0–7 are in the parent plan: 0 inventory and freeze, 1 skeleton, 2 dev loop, 3 image,
4 consumers, 5 skills, 6 retire the submodule, 7 distribution. Current position: **Phase 3
done, Phase 4 in progress, 5–6 blocked on the SciAgent-toolkit refactor.**

### Phase 8 — GESECA, and finishing CoReSh

Smallest phase, biggest unblocking effect, and it closes an open decision rather than opening
one.

| Step | Work | Gate |
|---|---|---|
| **G1** | Wire `coresh_match(pvalues = TRUE)` to `fgsea::geseca()` per §3.3; delete the guard and its message | p-values on real chunks; agreement within summed `log2err` on ≥7 of 8 datasets |
| **G2** | `coresh_loadings()` and `coresh_sets()` — the rest of C3 | Set names and memberships match DC_hum_verse's existing GMT |
| **G3** | `gsdb_coresh()` (C4) on `.ref_path("coresh")`, recording the snapshot tag | Byte-level agreement on set contents |
| **G4** | `gs_coregulation()` as a first-class verb returning `gs_result` | `gs_plot_dot()` and `gs_plot_bar()` render it with no renderer change; golden baseline added |
| **G5** | Cross-check `gesecaSimple()` against `geseca()` at moderate p as an independent guard | Guard asserts on the result, per the house rule |

### Phase 9 — TF and pathway activity

The phase with the most drift to absorb and the most methodological value to rescue.

| Step | Work | Gate |
|---|---|---|
| **T1** | `.ref_path("omnipath", ...)` tier for CollecTRI, PROGENy and DoRothEA, with the snapshot and the human→mouse mapping counts recorded | Networks resolve with OmniPath unreachable; STING-JR's local rebuild reproduced |
| **T2** | Write **ADR-005** and settle ULM versus MLM as the default, with the other exposed by name | The reason is written down before the code is |
| **T3** | `tf_activity()` and `pathway_activity()` returning the ADR-002 schema with a new `entity_type` | STING-JR's `master_tf_activities.csv` reproduced from the package |
| **T4** | `activity_consensus()` — the `decouple(consensus_score = TRUE)` ladder as one call | Rank cascade from `03_decoupler_tf.R:243-364` reproduced |
| **T5** | `regulon_null()` — curveball rewiring preserving size and in-degree, and exact within-block label permutation | `19_regulon_nulls.R`'s published numbers reproduced, including the `1/64` resolution floor |
| **T6** | `gs_plot_*` renderers for activity results | Reuse, not new plotting code |

**T5 is the reason this phase is worth its cost.** Everything above it is convenience;
the nulls are the part that stops a TF result from being a confident artefact, and today they
exist in exactly one project.

### Phase 10 — WGCNA

| Step | Work | Gate |
|---|---|---|
| **W1** | `wgcna_network()` — soft threshold, `blockwiseModules`, eigengenes, hub genes, with the small-n defaults from `2.8` and the reasoning in the documentation | AdaW-eWAT-WL run1 modules reproduced |
| **W2** | `wgcna_sweep()` over `mergeCutHeight`, `deepSplit`, `minModuleSize` | The project's parameter sweep reproduced from one call |
| **W3** | `wgcna_stability()` — eigengene co-correlation, Jaccard overlap, hub stability, preservation statistics | `2.11`'s 893 lines reproduced; the report is the deliverable |
| **W4** | `wgcna_modules_to_gsdb()` — modules become a `gs_db` | `gs_test()` and `gsdb_*` work on modules with no new machinery |
| **W5** | Renderers for eigengene heatmaps and module-trait relationships | `gs_plot_*` and `gs_theme()` |

**W4 is the architectural payoff.** A WGCNA module is a gene set. Once it is a `gs_db`, the
per-module GO analysis of `2.14` is `gs_test()` and needs no new code at all.

### Phase 11 — the uniform surface

Only worth doing once 8–10 exist, and it is what makes the package a zero-token
architecture rather than a collection of wrappers.

| Step | Work | Gate |
|---|---|---|
| **U1** | One vignette per stage, each runnable on the shipped fixture | `R CMD check` builds them |
| **U2** | A `contrast` object carrying design and caveat, from the DC_hum_verse four-family model | F3 cannot be read without its `CROSS_STUDY` flag |
| **U3** | `bulkirna_manifest()` — resolved snapshots, package version, seeds, in one provenance record per run | ADR-001 and ADR-004 become one printable object |
| **U4** | Rewrite the skills onto the package; retire the `source()` instructions | Link-integrity test; blocked with C5/C6 |

---

## 5. Sequencing, and what to do first

**Phase 8 first, and G1 before anything else.** It is a deletion plus a redirect, it closes a
decision I wrongly opened, and it unblocks both the CoReSh p-value path and `gs_coregulation()`.

Then **finish Phase 4** — STING-JR, then DC-nexus — before Phase 9. Reason: STING-JR *is* the
TF reference project. Migrating its GSEA stages first puts its conventions in front of us
while the TF API is being designed, and its `00c_prepare_networks.R` is the T1 specification.

Then Phase 9, then 10. Phase 11 last, and partly blocked.

**Not blocking anything, still open, all recorded elsewhere:**

- The image pin still says `bulkiRNA@v0.4.0` (`scbio-docker/docker/base/R/install_core.R:174`).
- `qs2` is absent from `scdock-r-dev:v0.5.13`; any chunk reader needs it.
- `10_gatom_modules.R`'s combined path still calls `gatom::` directly.
- `14839-DM-cGAS` is uncommitted, with the moved GATOM numbers in it.
- The msigdbr re-run's 47,988 → 72,408 MSigDB rows are still unread.

---

## 6. The tradeoffs of not doing this

Stated plainly, because a roadmap that only lists benefits is not a plan.

**Cost.** Four heavy optional dependencies, WGCNA the heaviest. Roughly three times the
current package surface. Every new default is a new thing that can be wrong in one place
instead of wrong in four — which is better, but it is also *load-bearing*, so a mistake
propagates further and faster than a mistake in a script does.

**Alternative — leave TF, PROGENy and WGCNA as scripts and keep the package to gene sets.**
Honest case for it: those stages are thin wrappers over decoupleR and WGCNA, and the
boundary rule arguably sends thin wrappers to a skill. Honest case against, and the reason
this is rejected: they are demonstrably **not** being done identically, the estimator choice
is already inconsistent across three projects, the networks already failed to fetch, and the
nulls that make the results defensible exist in one place and cannot be reused. That is the
same evidence that justified extracting CoReSh, one layer up.

**Alternative — a sibling package per stage.** Cleanest dependency isolation, and rejected
for the same reason as ADR-003's split: one maintainer, one release train, and a result that
cites three versions is harder to reproduce than one that cites one.

**Alternative — a config file of preferences instead of defaults in code.** Rejected: a
default in code is tested, versioned and visible in `help()`. A YAML file of preferences is
prose that drifts, which is the problem this whole refactor exists to fix.
