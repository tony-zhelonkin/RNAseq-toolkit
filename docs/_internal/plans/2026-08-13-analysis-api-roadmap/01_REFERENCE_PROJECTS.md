# Reference implementations: what the API has to reproduce

**Date:** 2026-08-13 · **Status:** durable reference; update when a better implementation appears
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

This document exists so that nobody — human or agent — has to go looking for how these
analyses were done before. Each section names **one canonical reference implementation**, the
conventions it established, the defects it carries, and the package surface that should
replace it.

Read this before designing an API. Do not re-derive a convention that is already recorded
here; if you disagree with one, say so and change it here first.

**How to use a reference project.** Read it. Do not run it, do not edit it, and do not commit
in it. These are live research repositories holding real data, several of them with
uncommitted work. Mount them read-only.

---

## 0. Conventions that hold across every stage

These came out of the projects below and are already load-bearing in `bulkiRNA`.

| Convention | Value | Where it came from |
|---|---|---|
| Compute and render are separate scripts | `NN_thing.R` + `NN_thing_viz.R` | Every STING-JR stage; the reason `gs_plot_*` never computes |
| One long master table per entity class | `master_gsea_table.csv`, `master_tf_activities.csv`, `master_de_table.csv` | DC_hum_verse `1.5`, STING-JR `03` |
| `contrast` is the join key across every stage | a string, stacked vertically | `coresh-slice/README.md` |
| Rank metric is the moderated t-statistic | `gsea.rank_metric` in project config | STING-JR `03:28`, WL `config.R` |
| DE input is limma-trend or edgeR QL | never raw fold change | STING-JR `02_de_limma_trend.R`, DC_hum_verse `1.1` |
| Everything expensive is cached to `checkpoints/` | keyed by stage, read back on re-run | all three projects |
| A stage prints what it did, with counts | `message()` per step, not silence | all three projects |
| Seeds are set from config, never hardcoded per call | `GSEA_SEED %||% 123` | STING-JR `03:34` |

The master-table schema itself is owned by
[ADR-002](../../adr/ADR-002-master-table-schema.md) and implemented by `gs_to_master()` /
`gs_validate_master()`. **A new entity class gets a new `entity_type` value, not a new
schema.**

---

## 1. CoReSh coregulation search

**Canonical reference:**
`/scratch/current/antonz/projects/DC-nexus/DC_hum_verse/02_analysis/scripts/coresh-slice`

The most complete pipeline in the estate, and the one whose `README.md` is worth reading in
full before touching anything CoReSh. 15 scripts, 5,125 lines counting the frozen `.slice/`
copy.

| Script | Lines | What it owns |
|---|---|---|
| `0_feasibility_matrix.py` | — | Scans `obs`, emits `contrasts_manifest.csv`, drops underpowered contrasts before any compute |
| `1.1_pseudobulk_de.R` | 435 | edgeR QL per contrast, four contrast families |
| `1.2_msigdb_gsea.R` | 260 | Hallmark, KEGG, Reactome, GO_BP |
| `1.7_decoupler_tf_analysis.R` | 238 | CollecTRI ULM |
| `1.8_progeny_analysis.R` | 236 | PROGENy MLM |
| `1.10_custom_iron_dc_gsea.R` | 269 | Custom programs + MitoPathways |
| `1.11_coresh_queries.R` | 475 | CoReSh Q1–Q16 |
| `1.12_coresh_derived_gsea.R` | 382 | CoReSh-derived sets → GSEA |
| `1.5_finalize_master_tables.R` | 256 | Power tiers, `caveat_flags`, GO_BP filter |
| `lib/{coresh_batch,extract_gene_loadings,symbols_to_entrez}.R` | 269 | The engine, md5-identical in four checkouts |

**The four contrast families, and why F4 exists.** This is the single best piece of design
thinking in the estate and the package should not flatten it. F1 cluster-vs-rest, F2 paired
tumour-vs-adjacent, F3 tumour-vs-healthy with batch unidentifiable, F4
adjacent-vs-healthy as an explicit estimate of the study-effect direction. Per gene,
`logFC_F3 ≈ logFC_F2 + logFC_F4`, so plotting F2 against F4 separates defensible biology
from study artifact. **A contrast carries a design and a caveat, not just a name.**

**Chunk tree.** `/data2/users/shared/refcache/coresh/current/preprocessed_chunks/{hsa,mmu}`,
snapshot `syn66227307_20260721`, 174 chunk files, ~21 GB. The project resolves it through
`CORESH_CHUNKS` with a legacy in-tree fallback. **The package resolves it through
`.ref_path("coresh", "preprocessed_chunks", <hsa|mmu>)` instead**
([ADR-004](../../adr/ADR-004-reference-data-tiers.md)), which needs no new environment
variable and records the snapshot tag.

**Official web UI:** <https://alserglab.wustl.edu/coresh>. This is the smoke test — run a
query there before believing a local sweep, which is what the skill already tells you to do.

**Already ported** (`R/coresh.R`, 2026-08-13): `coresh_chunks()`, `coresh_match()`,
`coresh_search()`, `coresh_convergence()`, `coresh_validate()`. `pct_var` verified
bit-identical to hand computation on 12 of 12 real datasets.

**Still to port:** `coresh_loadings()` and `coresh_sets()` from `extract_gene_loadings.R`,
then `gsdb_coresh()`. See the roadmap.

---

## 2. GESECA

**Upstream:** `alserglab/fgsea` — `R/geseca-multilevel.R`, `R/geseca-simple.R`,
`R/geseca-utils.R`, `R/geseca-plot.R`, `src/geseca.cpp`. Tutorial:
`vignettes/geseca-tutorial.Rmd` at commit `d466383`.

GESECA asks a different question from GSEA and it is the question CoReSh is built on: not
"is this set enriched at one end of one ranking" but **"is this set coregulated across the
samples of this experiment"**. It needs an expression matrix, not a ranked vector, and it
therefore needs no contrast at all — which makes it the natural companion to a DE-driven
GSEA rather than a competitor.

`fgsea` exports `geseca()` (multilevel Monte Carlo) and `gesecaSimple()` (plain `nperm`
permutation), plus `plotGesecaTable()`. Both are public. There is **no reason for
`bulkiRNA` to touch `fgsea:::gesecaCpp`** — see the roadmap for the measurement that
established this and corrected an earlier wrong conclusion of mine.

**Mechanics that bite, all verified on real chunks:**

- `center = TRUE` is the default. CoReSh chunks are already centred, so pass
  `center = FALSE` or the variance is computed twice.
- `checkGesecaArgs()` rejects duplicate `rownames(E)`. Every real CoReSh chunk has them,
  because Entrez ids repeat. `make.unique()` first; it changes no value.
- `sampleSize` controls the multilevel estimator's precision and **is honoured**:
  `log2err` went 3.04 → 1.41 → 0.625 for `sampleSize` 21 → 101 → 501. It looks ignored if
  you only test a null matrix, because there the answer comes from the pre-permutation
  screen and multilevel never escalates.
- Two runs never give the same p-value. Compare within `log2err`, never for equality.
- `log2err` is `Inf` for the most extreme sets, meaning the estimate is past reliable
  resolution. Report it rather than hiding it.

---

## 3. TF activity and pathway activity — the drift is already here

Three projects, three different answers to the same question. This is the §2-of-the-CoReSh-plan
pattern repeating in a layer the package has not reached yet, and it is the strongest argument
for doing this stage next.

| | DC_hum_verse `1.7`/`1.8` | AdaW-eWAT-WL `2.5`/`2.6` | STING-JR mouse_anchor `03` |
|---|---|---|---|
| TF network | `get_collectri("human")` live | **`get_dorothea("mouse", ABC)`** | **local rebuild**, `00c_prepare_networks.R` |
| TF estimator | **`run_ulm`** | **`run_mlm`** | `run_ulm(.mor = "mor", minsize = 5)` + `decouple(consensus_score = TRUE)` |
| PROGENy | `get_progeny("human", top = 500)`, `run_mlm` | `get_progeny("mouse", top = 500)`, `run_mlm` | local rebuild, MLM |
| `minsize` | not set | `5` | `5` |
| Multiple testing | BH within contrast | per contrast | BH within contrast |

**Two estimators for one question, three network provenances, and one upstream breakage.**

**The upstream breakage is real and documented.** STING-JR
`02_analysis/scripts/00c_prepare_networks.R:12-15`: decoupleR 2.16.0 with OmnipathR 3.18.4
— both current Bioconductor release — **fail to fetch CollecTRI or PROGENy**. OmnipathR's
web wrapper falls back to a static table whose loader errors, and
`OmnipathR::collectri()` dies on an internal `ncbi_tax_id` join. The project reconstructs
the networks locally, maps human to mouse through `babelgene` orthology with a title-case
fallback, records the mapping counts, and caches to `net_{collectri,progeny,dorothea}_mouse.rds`.

That is exactly the refcache tier from ADR-004, hand-built in one project. **The networks
must become a recorded reference-data source, not a live fetch**, or every TF result depends
on whether OmniPath was up that day.

**Nulls worth more than the estimator.**
`/scratch/current/antonz/projects/STING-JR/human_treg_arthritis/02_analysis/scripts/19_regulon_nulls.R`
is the most rigorous file in the estate. Two nulls the size-matched ones cannot express:

- **Curveball rewiring** of the TF→target bipartite graph, preserving every regulon's size
  *and every target's in-degree*, so a drawn regulon oversamples promiscuous high-|t| genes
  exactly as much as the observed one does. The null's centre therefore sits away from zero,
  and beating it means something a size-matched null cannot test.
- **Exact within-donor label permutation**: six donors carrying both arms, all `2^6 = 64`
  configurations enumerated and the contrast refitted end to end. No seed enters the
  p-value. The price is resolution — the finest one-sided p is `1/64 = 0.0156` — and the
  script publishes that rather than hiding it.

Companion stage 18 (`18_tf_activity.R`, 41 KB) establishes the ladder those nulls complete:
does the activity survive a network swap, an estimator swap, and does it scale with regulon
size. **That ladder is the deliverable, not one ULM call.**

Also worth keeping: `03_decoupler_tf.R`'s rank-cascade forensics (`03:243-364`) tracing one
TF's rank across regulon swap → MLM → consensus, and `03d_ulm_mechanic.R`, which is
lookup-only by design — it explains a score using the existing t-vector and refits nothing.

---

## 4. WGCNA

**Canonical reference:**
`/data1/users/antonz/projects/AdaW-eWAT-WL-bulkRNAseq/02_analysis/wl_subset`

| Script | Lines | What it owns |
|---|---|---|
| `2.8_wgcna_network.R` | 563 | `pickSoftThreshold` → `blockwiseModules` → eigengenes → hub genes |
| `2.9_wgcna_tables.R` | 282 | Module tables |
| `2.10_wgcna_comparison.R` | — | Cross-run comparison |
| `2.11_wgcna_robustness.R` | **893** | Eigengene co-correlation, Jaccard overlap, hub stability, WGCNA preservation statistics |
| `2.14_wgcna_go_comparison.R` | 517 | Per-module GO |
| `3.11_wgcna_viz.R`, `3.12_heatmap_wgcna_eigengenes.R` | 398 + 264 | Render only |

**Parameters, and the reasoning behind them** (`2.8:8-15`): n = 12 samples, below the
recommended n ≥ 15, so — signed network to preserve directionality, R² threshold relaxed to
**0.70** from the standard 0.80, `minModuleSize` **25** from the standard 30, `deepSplit` 2,
`mergeCutHeight` 0.20, `maxBlockSize` 15000, VST or logCPM input, and **effect sizes
(|r| > 0.6) over p-values**. The script takes `RUN_NAME MERGE_CUT DEEP_SPLIT MIN_SIZE` on
the command line precisely so the sweep is a first-class operation.

**The important part is 2.11, not 2.8.** A single WGCNA run is a hypothesis; the module
structure is only real if it survives the parameter sweep. 893 lines of robustness against
563 of network construction is the correct ratio and the API must preserve it: **a sweep and
a stability report are the deliverable, one `blockwiseModules()` call is not.**

---

## 5. Where each pain point in this document already has a package answer

| Reference-project pain | Package answer |
|---|---|
| Engine vendored four times by `cp` | `R/coresh.R`, one tested copy |
| Two `parse_gmt`, two ranked-vector builders | `gsdb_from_file()`, `gs_ranks()` |
| Three master-append implementations | `gs_to_master()`, `gs_validate_master()` |
| `neg_log_padj` capped at 16 in one place, absent in another | ADR-002, one expression |
| Chunk path via `CORESH_CHUNKS` plus in-tree fallback | `.ref_path()`, snapshot recorded |
| Confounder regexes retyped in prose and in a script | `filter_confounder_genes()` |
| `sym2ent`/`ent2sym` in the global environment | `gene_to_entrez()`, `entrez_to_gene()` |
| Unseeded `scoreGraph()` | `gatom_module(seed = )` |
| **Two TF estimators, three network sources** | **not yet — see roadmap Phase 9** |
| **WGCNA parameters justified in a comment** | **not yet — see roadmap Phase 10** |
