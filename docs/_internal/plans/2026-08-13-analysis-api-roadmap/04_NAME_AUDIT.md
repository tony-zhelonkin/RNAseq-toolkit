# S2/S3 name and argument audit

**Date:** 2026-08-19 · **Status:** implemented; R gates pending
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

## Measured surface

The registry currently contains 78 exports: 46 stable, 11 experimental, and 21
deprecated. The signature freeze is a separate axis with 24 names. Three frozen
names remain live: `build_dge()`, `ensure_dir()`, and
`format_pathway_name()`. The other 21 are deprecated shims.

## S2: the 13 inherited unprefixed exports

| Export | Decision | Reason |
|---|---|---|
| `read_counts_matrix()` | Keep stable | The specific I/O verb already says what it reads; an `io_` prefix adds no useful distinction. |
| `read_metadata()` | Keep stable | It is the companion specific I/O verb to `read_counts_matrix()`. |
| `ensure_dir()` | Keep stable | Its signature is frozen and the inventory names exactly two unmigrated consumers. Deprecating it would create migration work with no replacement benefit. |
| `write_session_provenance()` | Keep stable | The full verb describes a cross-layer reproducibility artifact rather than one analysis layer. |
| `build_dge()` | Keep stable | This frozen constructor creates the object consumed by the `de_*` renderer layer; it is not itself a renderer. |
| `annotate_genes()` | Keep stable | Annotation precedes and serves more than differential expression, so `de_annotate_genes()` would claim ownership the DE layer does not have. It has no inventory call sites, so this decision is conceptual. |
| `gene_to_entrez()` | Keep experimental | The direction is explicit and the conversion is cross-layer. It has no inventory call sites. |
| `entrez_to_gene()` | Keep experimental | It is the symmetric inverse of `gene_to_entrez()` and is already explicit. It has no inventory call sites. |
| `filter_confounder_genes()` | Keep experimental | Confounder filtering is a gene-level preparation step used across analysis layers. It has no inventory call sites. |
| `format_pathway_name()` | Keep stable | This frozen, general presentation verb is not computation and is used by multiple renderers. |
| `theme_bulki()` | Keep stable | A short house-theme name follows the established `ggplot2::theme_*()` vocabulary. |
| `bulkirna_check_deps()` | Keep stable | The `bulkirna_` prefix identifies a package-wide dependency check rather than an analysis-layer operation. |
| `download_gatom_references()` | Replace with `gatom_download_refs()` | The inventory contains no consumer call sites, so the rename creates no migration work. It also makes all six GATOM names use the `gatom_` prefix. The frozen old signature remains as a deprecated shim through 1.0.0. |

These call-site statements were checked directly against
`sciagent-rna/docs/data/used-functions.tsv`: `ensure_dir` is the only one of
the six named functions present, with exactly the two consumers recorded
above.

`bulkirna_api()` and `bulkirna_stochastic()` landed in S1 after the inherited
13-name inventory. They are deliberate package-metadata names and are covered
by the same enforcement allowlist.

## S3: measured argument vocabulary

This table was measured from the formals of all 57 live exports **plus the
formals of every S3 method registered for an exported generic**, resolved
through the namespace method table. Counts are functions exposing the formal
after S2; deprecated shims are excluded.

The method matters, and the first version of this sentence recorded the one that
had a gap: walking `formals()` on exported names alone cannot see an argument
declared on an S3 method behind a generic's `...`, which is how `group` and
`samples` on `gs_plot_heatmap.gs_matrix` went unaudited. **Every count below was
rechecked under the corrected method and none moved** — no method conceals a
`species`, `dir` or other cross-layer formal; only `top_n` and `palette` gain
occurrences, by two each, and neither appears in this table. So the numbers were
right while the technique was not, and those are separate claims.

| Concept | Canonical formal | Live exports | Finding |
|---|---|---:|---|
| Target species | `species` | 13 | One spelling, but six normalization implementations disagreed. They now use `.species()`. |
| Gene-set object | `db` | 5 | Consistent. `database` remains separate: it is a registry key or result filter, not a `gs_db` object. |
| Contrast | `contrast` | 2 | Consistent. |
| Random seed | `seed` | 4 | Consistent on explicit formals. `gs_test()` retains its established seed-through-`...` adapter interface. |
| Quiet reporting | `quiet` | 3 | Consistent. |
| Progress reporting | `verbose` | 8 | Consistent. |
| File path | `path` | 7 | Consistent. Frozen `ensure_dir(path)` remains the one directory operation using `path`. |
| Minimum set size | `min_size` | 7 | Consistent. `gs_test()` accepts method-specific size bounds through `...`. |
| Maximum set size | `max_size` | 7 | Consistent. `gs_test()` accepts method-specific size bounds through `...`. |
| Worker count | `n_cores` | 1 | Consistent. |
| Directory | `dir` | 4 | `download_gatom_references(dest_dir)` was the only live variant. New `gatom_download_refs(dir)` is canonical; the frozen shim retains `dest_dir`. |

`db_species` is deliberately not renamed: it is the source MSigDB organism
code, distinct from the target `species`. Likewise `chunk_dir` is a specific
CoReSh directory, not a generic output directory.

The enforcement test now collects every explicit formal from every live export and
requires it to belong to the curated vocabulary. It also requires the canonical
cross-layer names in the table above to remain present. A new spelling such as
`target_dir`, `sp`, or `rand_seed` therefore fails until the audit explicitly
adopts or rejects it. The public formals `B_cutoff`, `baseMean`, `log2FC`, and
`max.overlaps` remain reasoned exceptions to the package's snake-case convention
because they follow upstream vocabulary.

## S3 owner decision: one spelling for limits and palettes

On 2026-08-19, the owner chose `top_n` for the result-limit concept and
`palette` for the colour-palette concept. Named calls using the previous
spellings are breaking changes; no argument-level compatibility aliases were
added.

The decision rested on two measurements. None of the affected explicit
signatures is frozen, and `used-functions.tsv` contains zero call sites for
every affected export. The result-limit rename covers `gs_leading_edge()`,
`gs_plot_bar()`, `gs_plot_dot()`, `gs_plot_running()`, `coresh_loadings()`, and
`coresh_sets()`. The `coresh_sets()` database-provenance field changed from
`n_top` to `top_n` with its formal. The palette rename covers `de_bfc_plot()`,
`de_md_plot()`, `de_pca()`, `de_volcano()`, `gs_plot_bar()`, and
`gs_plot_dot()`.

Source inspection found one export the original `formals()` measurement missed:
`gs_plot_heatmap()` accepts both concepts through `...`, with `top` and
`colours` declared on its `gs_result` and `gs_matrix` methods. It is stable,
not frozen, and also has zero inventory call sites. Both method formals were
renamed in place to `top_n` and `palette`. The enforcement test now checks these
method formals explicitly so the generic's `...` cannot hide a retired spelling.

## Shared species contract

`.species()` accepts `human`, `Homo sapiens`, `hsa`, `mouse`, `Mus musculus`,
and `mmu`, case-insensitively and with spaces or underscores normalized. It
returns seven pinned fields: common and scientific names, three-letter code,
organism annotation package, biomaRt dataset, and two GATOM representations.

`annotate_genes()` previously inherited partial matching from `match.arg()`.
That behavior is retained deliberately for case-sensitive scientific-name
prefixes, so `species = "Homo"` continues to resolve. The package does not
extend partial matching to lower-case abbreviations that never worked there.

`.gsdb_species_label()` calls `.species(allow_custom = TRUE)`. Known human and
mouse aliases normalize to scientific names; an unknown non-empty label remains
valid and only has underscores changed to spaces. This preserves user-supplied
custom species in `gs_db` providers while keeping validating callers strict.
Custom labels do not use partial matching: `"Mus"` remains `"Mus"`. Partial
matching exists only on the strict path used by `annotate_genes()`, where
`"Homo"` continues to resolve to `"Homo sapiens"`.

The audit brief named five sites. Inspection found a sixth: the downloader had
its own direct `switch()` in addition to `.gatom_species()`. It was migrated too;
leaving it would have made the claimed single owner false.

## Enforcement scope

The layer, argument-vocabulary, and RNG-mutation tests parse the package's `R/`
sources. Those files are unavailable when `R CMD check` runs the tests against
the installed package, so the three checks skip there with an explicit
source-tree-only reason. They execute under source-tree test runs. The API
registry still checks itself against the installed namespace in both
directions under `R CMD check`.

## File organization

The GATOM downloader is not a gene-set provider and no `gsdb_*` function shares
its file. Its implementation now lives in `R/gatom-download.R`, mirrored by
`tests/testthat/test-gatom-download.R`, so the file name follows the exported
`gatom_` family.

## Deliberately unchanged

- No signature-frozen formal changed. In particular,
  `download_gatom_references(dest_dir, species, networks, overwrite)` remains
  callable as written, and `ensure_dir(path)` remains stable.
- The general `database` key/filter formal was not conflated with a `db` object.
- The target `species` formal was not conflated with MSigDB's `db_species`.
- No RNG-state mutation was added; `R/rng.R` remains the sole owner.
- No consumer project, golden baseline, fixture, generated manual page, or
  hand-written `NAMESPACE` entry was changed.
