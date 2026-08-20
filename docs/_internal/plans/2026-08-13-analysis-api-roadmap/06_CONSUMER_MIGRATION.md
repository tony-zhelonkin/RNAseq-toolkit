# Migrating a consumer to the package API

**Date:** 2026-08-20 · **Status:** live; update as each consumer lands
**Parent:** [00_ROADMAP.md](00_ROADMAP.md) · **Phase:** 4

This is the recipe for Phase 4, written so the next consumer does not rediscover what the first two
cost. Read §1 before touching a script, §2 for the call mapping, and §3 for the four ways a migrated
script fails while looking fine.

**Why Phase 4 is the critical path.** The 21 deprecated shims exist only because consumers still call
the old names. `v1.0.0` is defined as the release where they go. So every consumer migrated is a
direct step toward that release, and the last one unblocks both Phase 6 and the shim removal.

---

## 0. Scope, measured rather than estimated

`grep` for the 21 deprecated names across a project overcounts by roughly two orders of magnitude,
because the vendored `RNAseq-toolkit` submodule contains the **definitions** and is checked out once
per sub-project. Exclude these before counting:

| Exclude | Why |
|---|---|
| any path containing `01_modules/` or `01_scripts/RNAseq-toolkit` | the vendored submodule; these are definitions, retired by Phase 6, not call sites |
| `.ref/` | read-only copies of *other* projects; migrating them edits a copy of someone else's work |
| `_archive/` | frozen |
| `.slice/` | frozen prior-session pipelines, so marked in their own README |
| `SciAgent-toolkit` | out of bounds by standing instruction |

Measured with those exclusions:

| Project | Sub-repo | Files | Call sites | State |
|---|---|---|---|---|
| Meta-Aging | `14616-DM` | 1 | 10 | ✅ migrated, `d6633dd` |
| Meta-Aging | `14761-DM`, `14782-DM`, `neuroimmune-receptor-atlas` | 0 | 0 | nothing to do |
| DC-nexus | `DC_Dictionary` | 1 | 19 | ⬜ branch created |
| DC-nexus | `DC_hum_verse` | 3 | 18 | ⬜ |
| DC-nexus | `DC_mouse_cancer` | 0 | 0 | nothing to do |
| STING-JR | 4 sub-repos | 19 | 94 | 🚫 excluded by the owner |

STING-JR was surveyed and then excluded by decision. The survey is kept because it is the largest
consumer and the numbers will be wanted if that decision changes: `mouse_anchor` 10 files / 56 sites,
`human_treg_arthritis` 6 / 30, `human_pbmc_febrile` 2 / 4, `human_ra_synovium` 1 / 1.

**A project is not a repository.** All three consumers are superprojects of git submodules, each its
own repo, several on `dev` rather than `main`, and all with uncommitted work. "Migrate DC-nexus" means
committing in two live research repositories. Work on a `migrate/bulkirna` branch per sub-repo, stage
only the files you changed, and never `git add -A`.

---

## 1. Before editing anything: the two preconditions

**The installed library must be current.** This is the precondition that silently invalidates
everything else. The shared library at `/data1/users/antonz/pipeline/.rlib-bulkirna` held **0.3.1**
when Phase 4 resumed, which predates `gs_to_master()` and `gs_validate_master()` — the exact verbs
`normalize_gsea_results()` migrates to. A migration written against a stale library either fails
confusingly or, worse, silently uses an older definition. Check first:

```r
packageVersion("bulkiRNA"); setdiff(c("gs_ranks","gs_test","gs_to_master",
  "gsdb_load","gs_validate_master"), getNamespaceExports("bulkiRNA"))
```

**The consumer's own image must ship it.** This is where the first attempt went wrong, so the reasoning
is worth keeping.

`scdock-r-dev:v0.5.10` — which **8 of 9 projects** used — ships no `bulkiRNA` at all, so the migrated
script failed at `library(bulkiRNA)`. The first fix mounted the installed library read-only and
prepended it with `R_LIBS`. It worked and was still wrong: it depended on a host path outside the
project, and because the library is built for one R minor version it could not be applied to the
`dev-archr` service (R 4.4.1 against 4.5.3), so the two services in one compose file diverged. **A fix
that cannot be applied uniformly is a warning that it is the wrong fix.**

**Bump the image instead.** `scdock-r-dev:v0.5.13` ships `bulkiRNA`, and for these consumers the bump
is safe and checked rather than assumed:

- Seven packages are dropped between v0.5.10 and v0.5.13 — `EnhancedVolcano`, `GPArotation`, `mnormt`,
  `psych`, `reactome.db`, `sankey`, `simplegraph`. **Check these against the consumer**; neither
  Meta-Aging nor DC-nexus references any.
- `ggpubr` 0.6.3 → 1.0.0 and `rstatix` 0.7.3 → 1.1.0 are major bumps. Neither project uses either.
- 55 packages are added, including `bulkiRNA`, `gatom`, `WGCNA`, `homologene`, `org.Mm.eg.db` and
  `orthogene` — the last two of which DC-nexus needs.

**The residual version skew is the real open problem.** v0.5.13 ships **0.4.0** while the package tree
is at 0.5.0, so a consumer runs a different version from the one its migration was verified against —
which is the drift ADR-001 exists to prevent, reappearing one level down. Measured rather than
assumed: running the migrated script under both versions gives identical row counts, identical numeric
columns, identical `core_enrichment` and `genes_full_set`, and a byte-identical
`gsea_summary_stats.csv`. The only difference is **12 `pathway_name` values**, where 0.4.0 lowercases a
leading small word (`"the Nlrp3 Inflammasome"`) and 0.5.0 capitalises it. No published number moves.

Closing it properly is one line — `scbio-docker/docker/base/R/install_core.R:174`, currently
`"tony-zhelonkin/bulkiRNA@v0.4.0"` — plus a ~150-minute rebuild affecting every project on the shared
image. That is a standing owner decision, not a migration step.

**Known environment gap the bump does not fix:** `OmnipathR` is absent from both v0.5.10 and `v0.5.13`, and
the consumers' devcontainers add it in `postcreate`. Any script with a decoupleR
`get_collectri()`/`get_progeny()` block therefore cannot be run end to end in a bare container. This
is pre-existing and unrelated to the migration, but it caps what can be verified — say so rather than
reporting a partial run as a pass.

---

## 2. The call mapping

| Legacy | Package | Notes |
|---|---|---|
| `run_gsea(DE_results=, rank_metric=)` | `gs_ranks(x, metric=)` then `gs_test(ranks, db, contrast=)` | ranking and testing are separate verbs now; `gs_ranks()` defaults `genes` to `rownames(x)`, which matches the legacy preamble |
| `normalize_gsea_results()` | `gs_to_master(res, db=, universe=, entity_type=)` | see §3.1 — `db` is the **object** |
| `load_reference_db(database=, toolkit_dir=)$T2G` | `gsdb_load(database, species=)` | `mitopathways`, `mitoxplorer`, `mito_unified`, `transportdb` all ship in `inst/extdata`; no `toolkit_dir`, no network |
| `parse_gmx(file=, prefix=)` | `gsdb_from_file(path, prefix=, database_label=)` | |
| `parse_mitoxplorer()` | `gsdb_load("mitoxplorer")` | |
| `filter_by_size()` | `min_size`/`max_size` on the provider or on `gs_test()` | not interchangeable; see §3.4 |
| `gsea_dotplot()` | `gs_plot_dot()` | |
| `gsea_barplot()` | `gs_plot_bar()` | |
| `gsea_running_sum_plot()` | `gs_plot_running()` | |
| `plot_all_gsea_results()` | `gs_plot_*()` plus `gs_save()` | no single replacement, by design |
| `convert_human_to_mouse(T2G=)` | **no successor for arbitrary sets** — see §4 | |

Legacy fgsea knobs map onto `gs_test()`'s dots: `pvalue_cutoff = 1.0` becomes nothing at all, because
`gs_test()` returns every tested set; `nperm` becomes `n_perm_simple`; `eps` and `score_type` keep
their names; `seed` is passed through and pinned.

---

## 3. The four ways a migrated script fails while looking fine

Each of these was hit for real on the first consumer. Three were caught only by running the script.

### 3.1 `gs_to_master(db = )` takes the gene-set object, not a label

`normalize_gsea_results()` took `database = "Hallmark"`, a string. `gs_to_master()`'s `db` is the
`gs_db` itself, used with `universe` to derive `genes_full_set`. Passing the name raises
``` `db` must be a named list of gene vectors or `NULL` ```, which in a `tryCatch(..., NULL)` wrapper
becomes **zero rows and no error** — the project shape this whole effort exists to remove. Pass the
object.

### 3.2 The `database` column changes value, and downstream joins are keyed on it

`gsdb_msigdb()` labels itself from the collection: `msigdb_C3_TFT_GTRD`, not the project's
`CollecTRI_TF`. Every figure and table that joins on `database` breaks silently — the column is
present and populated, just with different strings. Relabel the **provider**, not the master rows, so
the per-cell checkpoints carry the project's key too:

```r
db <- gsdb_msigdb(species = "Mus musculus", collection = spec$collection,
                  subcollection = spec$subcollection)
attr(db, "database") <- db_name
```

### 3.3 Checkpoint files change format under unchanged names

`gs_test()` returns a `gs_result` tibble where `run_gsea()` returned a `clusterProfiler` `gseaResult`.
Every consumer caches these. A cached rerun then feeds the old format into `gs_to_master()`. Rename
the checkpoints — a `_bulkirna` suffix — which makes the change explicit and leaves the old files for
comparison. Do not delete them; they are the only baseline available for §5.

### 3.4 Size filters are applied at different points

`minGSSize`/`maxGSSize` on `clusterProfiler::GSEA()` filtered **after** intersecting sets with the
ranked genes. `gs_test(min_size=, max_size=)` filters at the same point; `gsdb_load(min_size=)`
filters the **raw** sets. They are not interchangeable, and moving the argument to the provider
quietly changes which pathways are tested. Keep the bounds on `gs_test()`.

---

## 4. `convert_human_to_mouse()` has no successor for this use

Recorded because it changes the `v1.0.0` removal plan.

`bulkirna_api()` gives its successor as `gsdb_msigdb(species = "Mus musculus", db_species = "HS")`,
which is correct for MSigDB and only for MSigDB. `DC_Dictionary` calls it on a T2G parsed from a
custom GMX, and **no export ortholog-maps an arbitrary `gs_db`.** The successor named in the registry
does not cover the call site.

For `DC_Dictionary` specifically this does not bite: the set being converted is MitoPathways, and
`gsdb_load("mitopathways")` ships 119 sets for both `Homo sapiens` and `Mus musculus`, so the
parse-then-convert step disappears rather than being translated. But the general gap is real, and one
of two things must be true before `v1.0.0`: either an ortholog verb exists for arbitrary sets, or
the registry's `superseded_by` for `convert_human_to_mouse` is narrowed to say it covers MSigDB only.

---

## 5. How to verify, and what the verification is worth

Running the migrated script is not optional. Three of the four defects in §3 pass a reading and pass
`parse()`.

The full job is usually too large: `14616-DM` is 164 contrasts × 8 collections = 1,312 GSEA runs. Cut
the **input**, not the script. Build a scratch tree that mirrors the project layout, mount it at the
path the script hardcodes, and reduce the DE checkpoint to two contrasts. The script stays
byte-identical to what is committed. Never point a verification run at the live `03_results/`.

**Then compare against the pre-migration master table, and separate the two causes of disagreement.**
Raw comparison on `14616-DM` looked alarming — NES differing by up to 4.7. It resolved cleanly once
the rows were split by whether `set_size` changed:

- **On version-pinned sets — MitoPathways and MitoXplorer, bundled in the package, so no reference
  drift is possible — NES and padj reproduce exactly: max absolute difference 0 across 272 rows, and
  `genes_full_set` identical 272 of 272.** This is the measurement that establishes the translation is
  faithful.
- Across all rows whose `set_size` was unchanged, median absolute NES difference is 9e-4.
- The remaining differences track MSigDB content drift under `msigdbr 26.1.0`, where only **11%** of
  GO_BP sets kept their previous size.

Pick a database with no reference drift and require exact agreement there. A global tolerance over
everything would have hidden both the real defects and the real drift in the same number.

---

## 6. Behaviour changes to expect, and to write down

These are not defects. They are the package's definitions replacing the project's, which is the point
of ADR-002 — but they move published numbers, so they belong in the commit body.

1. **`neg_log_padj` loses any cap.** `gs_validate_master()` rejects the capped form outright:
   *"the retired cap-at-16 convention is not allowed."* Rows with `padj < 1e-16` stop reading exactly
   `16` and carry their true value, up to ~308. `padj` itself does not move. Attempting to preserve
   the project's cap fails the gate, which is the ADR holding.
2. **`core_enrichment` becomes fgsea's leading edge** rather than `clusterProfiler`'s. On the 272
   exactly-reproducing rows it differs on 52 — smaller on 46, larger on 6, a subset of the old value
   on 266 — while NES and padj are identical on all 272. So the statistic is unchanged and the
   leading-edge *definition* is not.
3. **Column order is not the schema record's row order.** `inst/extdata/master-schema-v1.csv` lists
   `entity_type` last because it is the optional column; `gs_to_master()` emits it **first** and
   `gs_validate_master()` requires that. Read the columns from the CSV rather than a literal, then
   move `entity_type` to the front.
4. **Stale MSigDB collection strings surface here, loudly.** `gsdb_msigdb()` raises
   `Unknown subcollection` rather than returning an empty set. On `14616-DM` this caught a config
   comment that was wrong in both directions: under `msigdbr 26.1.0`, `CP:KEGG` must become
   `CP:KEGG_LEGACY` (186 sets), while `CP:WIKIPATHWAYS` did **not** move (925 sets), contrary to the
   note in the file predicting both would.

**No exported accessor returns the schema columns.** Consumers must read the shipped CSV or hardcode
the names. Worth closing.

---

## 7. What stays hand-written after migration

`bulkiRNA` has no TF-activity or pathway-activity verb, so the CollecTRI/ULM and PROGENy/MLM blocks
stay as they are and are still aligned to the master schema by hand. Leave them alone and validate the
assembled table with `gs_validate_master()`, which is the guard that makes hand-built rows safe. This
is the concrete cost of parking Phases 10–11, and the clearest argument for eventually unparking them:
every consumer pays it separately.
