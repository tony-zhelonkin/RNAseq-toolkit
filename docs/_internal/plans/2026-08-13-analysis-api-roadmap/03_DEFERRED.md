# Deferred: the activity and network layers

**Date:** 2026-08-13 · **Status:** parked by decision, not abandoned. Do not start these.
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

TF activity, PROGENy and WGCNA were drafted as Phases 9–11 and are **parked** so that the
current API surface can stop moving first. The reasoning is in the roadmap §0 and §6.

This document keeps the survey findings, because they cost real reading to obtain and they
will be needed unchanged when the work resumes. The reference implementations and their
parameter reasoning are in [01_REFERENCE_PROJECTS.md](01_REFERENCE_PROJECTS.md) §3 and §4.

**Nothing here is scheduled. No ADR is owed for any of it until it is.**

---

## Why it will eventually be worth doing

The drift these layers already carry, measured across three projects:

| | DC_hum_verse `1.7`/`1.8` | AdaW-eWAT-WL `2.5`/`2.6` | STING-JR mouse_anchor `03` |
|---|---|---|---|
| TF network | `get_collectri("human")` live | `get_dorothea("mouse", ABC)` | local rebuild |
| TF estimator | `run_ulm` | `run_mlm` | `run_ulm(.mor = "mor", minsize = 5)` + consensus |
| PROGENy | `top = 500`, MLM | `top = 500`, MLM | local rebuild, MLM |
| `minsize` | not set | 5 | 5 |

Two estimators for one question, three network provenances, one inconsistent `minsize`.

**And the upstream breakage is documented, not hypothetical.** STING-JR
`02_analysis/scripts/00c_prepare_networks.R:12-15`: decoupleR 2.16.0 with OmnipathR 3.18.4 —
both current Bioconductor release — **fail to fetch CollecTRI or PROGENy**. OmnipathR's web
wrapper falls back to a static table whose loader errors, and `OmnipathR::collectri()` dies on
an internal `ncbi_tax_id` join. That project reconstructs the networks locally, maps human to
mouse through `babelgene` orthology with a title-case fallback, records the mapping counts, and
caches the result. That is ADR-004's refcache tier, hand-built once.

**The part worth the most is the nulls, not the estimator.**
`STING-JR/human_treg_arthritis/02_analysis/scripts/19_regulon_nulls.R`:

- **Curveball rewiring** of the TF→target bipartite graph, preserving every regulon's size
  *and every target's in-degree*, so a drawn regulon oversamples promiscuous high-|t| genes
  exactly as much as the observed one does. The null's centre therefore sits away from zero and
  beating it tests something a size-matched null cannot.
- **Exact within-donor label permutation** — six donors with both arms, all `2^6 = 64`
  configurations enumerated, the contrast refitted end to end, no seed entering the p-value.
  The finest attainable one-sided p is `1/64 = 0.0156`, and the script publishes that limit
  rather than hiding it.

Those exist in one file and nothing else can reach them. That is the strongest single argument
for eventually doing this work, and it will still be true later.

---

## The shape it would take, when it resumes

Recorded at the level of detail that survives without being a commitment.

**Networks first.** A `.ref_path("omnipath", ...)` tier for CollecTRI, PROGENy and DoRothEA,
recording the snapshot and the human→mouse mapping counts. Reproducing STING-JR's local rebuild
is the acceptance test. Nothing else can be trusted until the networks stop depending on
whether OmniPath was up that day.

**Then the two verbs**, `tf_activity()` and `pathway_activity()`, returning the ADR-002 schema
with a new `entity_type` rather than a new schema. STING-JR's `master_tf_activities.csv` is 13
columns and would need converging.

**Then `activity_consensus()`** — the `decouple(consensus_score = TRUE)` ladder, and the rank
cascade forensics from `03_decoupler_tf.R:243-364`, as one call.

**Then `regulon_null()`** — the two nulls above.

**WGCNA separately**, and the ratio in the reference project is the specification: 893 lines of
robustness against 563 of network construction. A single run is a hypothesis. `wgcna_sweep()`
and `wgcna_stability()` are the deliverable; `wgcna_network()` is the easy part. The
architectural payoff is `wgcna_modules_to_gsdb()` — a module is a gene set, so once it is a
`gs_db` the per-module GO analysis of `2.14_wgcna_go_comparison.R` is `gs_test()` and needs no
new code.

---

## Open questions to answer then, not now

1. **ULM or MLM as the default for TF activity?** Two projects disagree. The answer wants a
   comparison on real contrasts, not a preference. This was going to be ADR-005; the ADR is
   withdrawn until the work is scheduled, because an ADR for unscheduled work rots.
2. **Does the WGCNA dependency tree justify a sibling package?** `WGCNA` plus `impute`,
   `preprocessCore`, `dynamicTreeCut` and `fastcluster` is the heaviest thing either layer
   would add. ADR-003 rejected splitting on the grounds of one maintainer and one release
   train, and that still holds — but WGCNA is the first dependency that makes the case
   arguable, and a second one like it should reopen the decision.
3. **Does a `contrast` object belong in the package?** The DC_hum_verse four-family model
   (`coresh-slice/README.md`) carries a design and a caveat with each contrast, so that F3
   cannot be read without its `CROSS_STUDY` flag. That is a good idea whose natural home is
   unclear: it touches DE, GSEA and activity equally, and it may be a skill rather than a type.
