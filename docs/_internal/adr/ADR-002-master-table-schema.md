# ADR-002 — `bulkiRNA` owns the master-table schema and validator, not the file

**Status:** Accepted 2026-08-12 · **Decided by:** Anton Zhelonkin

## Context

The master GSEA table is a 14-column contract crossing every analysis project
(`pathway_id`, `pathway_name`, `database`, `contrast`, `nes`, `pvalue`, `padj`, `set_size`,
`leading_edge_size`, `gene_ratio`, `core_enrichment`, `genes_full_set`, `direction`,
`neg_log_padj`) — and it was owned by nobody. The consequences, all confirmed:

- **three** drifted master-append implementations across the consumer projects;
- **one column carrying two incompatible saturation conventions** — `normalize_gsea.R:214`
  caps `neg_log_padj` at 16, while the migrated `05` and the renderer `12_gsea_viz.R` reach
  ~307.65 via `pmax(padj, .Machine$double.xmin)`, and `12_gsea_viz.R:487` *sizes points* by
  that column, so one figure can mix both scales;
- **four columns silently dropped** by a ten-name allowlist at
  `08_coresh_derived_gsea.R:215`, which is why 492 CoReSh rows carry `neg_log_padj = NA`
  next to a `padj` of 2.6e-33.

A data contract cannot be fixed by a code convention. It had drifted precisely because it
lived as prose in three places.

## Decision

**`bulkiRNA` owns the columns and their semantics. It does not own the file.**

| Owner | Owns | Rationale |
|---|---|---|
| **`bulkiRNA`** | the versioned schema (`inst/extdata/`), `gs_to_master()`, `gs_validate_master()` | it computes the columns, and it is the only layer with a test suite |
| **SciAgent-toolkit skill** | which contrasts to assemble, where the file goes, how to read it | judgement, no return value |
| **analysis repo** | the file itself | a deliverable, not a library artifact |

This is the boundary rule from the CoReSh plan §3 applied without exception: *a table in and
a table out is a package function; a result in and a decision out stays a skill.* A schema
validator has a return value, so it cannot live in a skill — and prose is exactly how it
drifted three ways.

The deciding argument is authorship: **whoever computes a column should own its
definition**, or the definition lives somewhere that cannot be tested. `bulkiRNA` has 810
tests. The toolkit has none.

**The `neg_log_padj` convention is settled here:**

```r
neg_log_padj <- -log10(pmax(padj, .Machine$double.xmin))
```

Information-preserving, and already what two of the three call sites do. The cap at 16 is
retired.

## Consequences

- **C0 is revised.** Do **not** lift `gs_to_master()` into
  `02_analysis/helpers/de_gsea_helpers.R`. Lift it into the package. That retires all three
  drifted implementations at once, and the saturation convention becomes a tested line
  rather than a rule three files each remember differently.
- The schema is **data in `inst/extdata/`, not code**, and it is versioned. A consumer with
  a different downstream format registers a different schema rather than forking the
  function.
- `gs_to_master()` **returns a tibble and writes nothing.** No directory layout, no append
  strategy, and no file path enters the package.
- `genes_full_set` still obeys MADR-008 — computed by the caller from the `gs_db`, which is
  why the `db` argument exists.
- `entity_type` is a consumer concept (pathway / metabolite / CoReSh hit), so it is an
  argument, not a hardcoded value.

## The honest counterargument

The master table is an *analysis deliverable*, and putting its schema into a
general-purpose library couples the library to one lab's downstream format. Two things make
that tolerable: the schema is registrable data rather than hardcoded logic, and the
function returns a tibble rather than writing a file. If a second consumer ever needs a
materially different shape, the schema registry — not a fork — is the seam.

## Rejected alternatives

- **SciAgent-toolkit owns it.** Rejected: a skill cannot be `library()`d, so the schema
  would again be prose or a vendored script, which is the failure being fixed.
- **The analysis-repo template owns it** (a schema file per project via `init-project.sh`).
  Rejected: it makes per-project drift the *default* rather than an accident.
