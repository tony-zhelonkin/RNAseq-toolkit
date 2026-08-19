# S2/S3 name and argument audit

**Date:** 2026-08-18 · **Status:** implemented; R gates pending
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

## Measured surface

The audit started from 74 exports: 54 live and 20 deprecated. This corrects the
52-live snapshot in the roadmap, which predated S1 and G1. S2 adds one canonical
export and moves its frozen predecessor onto the deprecation clock, so the
expected post-documentation surface is 75 exports: 54 live and 21 deprecated.
The signature-frozen live names at the start were `build_dge()`,
`download_gatom_references()`, `ensure_dir()`, and `format_pathway_name()`.

## S2: the 13 inherited unprefixed exports

| Export | Decision | Reason |
|---|---|---|
| `read_counts_matrix()` | Keep stable | The specific I/O verb already says what it reads; an `io_` prefix adds no useful distinction. |
| `read_metadata()` | Keep stable | It is the companion specific I/O verb to `read_counts_matrix()`. |
| `ensure_dir()` | Keep stable | Its signature is frozen and two unmigrated consumers call it directly. Deprecating it would create migration work with no replacement benefit. |
| `write_session_provenance()` | Keep stable | The full verb describes a cross-layer reproducibility artifact rather than one analysis layer. |
| `build_dge()` | Keep stable | This frozen constructor creates the object consumed by the `de_*` renderer layer; it is not itself a renderer. |
| `annotate_genes()` | Keep stable | Annotation precedes and serves more than differential expression, so `de_annotate_genes()` would claim ownership the DE layer does not have. |
| `gene_to_entrez()` | Keep experimental | The direction is explicit and the conversion is cross-layer. |
| `entrez_to_gene()` | Keep experimental | It is the symmetric inverse of `gene_to_entrez()` and is already explicit. |
| `filter_confounder_genes()` | Keep experimental | Confounder filtering is a gene-level preparation step used across analysis layers. |
| `format_pathway_name()` | Keep stable | This frozen, general presentation verb is not computation and is used by multiple renderers. |
| `theme_bulki()` | Keep stable | A short house-theme name follows the established `ggplot2::theme_*()` vocabulary. |
| `bulkirna_check_deps()` | Keep stable | Package metadata already carries the package name as its namespace. |
| `download_gatom_references()` | Replace with `gatom_download_refs()` | It is the sole member of a six-function GATOM family without the `gatom_` prefix. The frozen old signature remains as a deprecated shim through 1.0.0. |

`bulkirna_api()` and `bulkirna_stochastic()` landed in S1 after the inherited
13-name inventory. They are deliberate package-metadata names and are covered
by the same enforcement allowlist.

## S3: measured argument vocabulary

This table was measured from the formals of all 54 live exports. Counts are
functions exposing the formal after S2; deprecated shims are excluded.

| Concept | Canonical formal | Live exports | Finding |
|---|---|---:|---|
| Target species | `species` | 12 | One spelling, but six normalization implementations disagreed. They now use `.species()`. |
| Gene-set object | `db` | 4 | Consistent. `database` remains separate: it is a registry key or result filter, not a `gs_db` object. |
| Contrast | `contrast` | 2 | Consistent. |
| Random seed | `seed` | 3 | Consistent on explicit formals. `gs_test()` retains its established seed-through-`...` adapter interface. |
| Quiet reporting | `quiet` | 3 | Consistent. |
| Progress reporting | `verbose` | 6 | Consistent. |
| File path | `path` | 7 | Consistent. Frozen `ensure_dir(path)` remains the one directory operation using `path`. |
| Minimum set size | `min_size` | 5 | Consistent. `gs_test()` accepts method-specific size bounds through `...`. |
| Maximum set size | `max_size` | 5 | Consistent. `gs_test()` accepts method-specific size bounds through `...`. |
| Worker count | `n_cores` | 1 | Consistent. |
| Directory | `dir` | 4 | `download_gatom_references(dest_dir)` was the only live variant. New `gatom_download_refs(dir)` is canonical; the frozen shim retains `dest_dir`. |

`db_species` is deliberately not renamed: it is the source MSigDB organism
code, distinct from the target `species`. Likewise `chunk_dir` is a specific
CoReSh directory, not a generic output directory.

## Shared species contract

`.species()` accepts `human`, `Homo sapiens`, `hsa`, `mouse`, `Mus musculus`,
and `mmu`, case-insensitively and with spaces or underscores normalized. It
returns the common and scientific names, three-letter code, organism annotation
package, biomaRt dataset, and GATOM representations.

`annotate_genes()` previously inherited partial matching from `match.arg()`.
That behavior is retained deliberately for case-sensitive scientific-name
prefixes, so `species = "Homo"` continues to resolve. The package does not
extend partial matching to lower-case abbreviations that never worked there.

`.gsdb_species_label()` calls `.species(allow_custom = TRUE)`. Known human and
mouse aliases normalize to scientific names; an unknown non-empty label remains
valid and only has underscores changed to spaces. This preserves user-supplied
custom species in `gs_db` providers while keeping validating callers strict.

The audit brief named five sites. Inspection found a sixth: the downloader had
its own direct `switch()` in addition to `.gatom_species()`. It was migrated too;
leaving it would have made the claimed single owner false.

## Deliberately unchanged

- No signature-frozen formal changed. In particular,
  `download_gatom_references(dest_dir, species, networks, overwrite)` remains
  callable as written, and `ensure_dir(path)` remains stable.
- The general `database` key/filter formal was not conflated with a `db` object.
- The target `species` formal was not conflated with MSigDB's `db_species`.
- No RNG-state mutation was added; `R/rng.R` remains the sole owner.
- No consumer project, golden baseline, fixture, generated manual page, or
  hand-written `NAMESPACE` entry was changed.
