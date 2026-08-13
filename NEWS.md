# bulkiRNA 0.5.0

The package version is the unit of reproducibility (ADR-001), so the changes
below get a version of their own: they moved numbers in a consumer project.

## API

* `bulkirna_api()` adds a machine-readable registry of all public functions,
  keeping lifecycle status and the historical signature freeze as independent
  axes.

## GATOM

* `gatom_module()` gains a final argument, `gene2reaction_extra`. It reaches
  `gatom::makeMetabolicGraph(gene2reaction.extra = )` and is what the combined
  network needs to see all of its edges. On the 14839 combined graph it takes
  the edge count from 7,608 to 8,179 against an unchanged 6,495 vertices. A
  value that is neither `NULL` nor a data frame with `gene` and `reaction`
  columns is rejected with a message naming the fix.

* `gatom_module()` records `attr(m, "solution_weight")`, the weight
  `mwcsr::solve_mwcsp()` reports for the module it returned. The three-call
  sequence discarded it, so consumers reaching for the solver's own objective
  had nothing to read.

* `gatom_refs()` accepts `network = "combined"` alongside `"kegg"`, resolves
  `network.combined.rds` and `met.combined.db.rds`, passes the flavour through
  to `download_gatom_references()`, and names the network in both the argument
  check and the not-found message.

## Formatting

* `format_pathway_name()` capitalizes a function word when it opens a label:
  `VIA_NEGATIVE_REGULATION` reads `Via Negative Regulation` rather than
  `via Negative Regulation`. Mid-label function words stay lowercase, so
  `HALLMARK_TNFA_SIGNALING_VIA_NFKB` still reads
  `TNF-alpha Signaling via NF-kappaB`. Only the nine function words are
  affected; Greek letters, chemical prefixes, roman numerals and lowercase gene
  symbols such as `mTOR` are untouched.

  This is a behaviour change inside a frozen signature. The export freeze
  covers formals, not output.

## Reproducibility note for existing GATOM results

`gatom::scoreGraph()` is stochastic. `gatom_module()` has seeded it since it was
written, but hand-rolled three-call sequences generally do not, and results
produced that way are one unrepeatable draw — the same input scored 165.85 and
then 169.32 in a single session. Modules recomputed through `gatom_module()`
will differ from such caches. The new values are the reproducible ones.

# bulkiRNA 0.4.0 and earlier

These releases predate this file. See the git tags and `MIGRATION.md` for the
legacy-name shims and the 24 frozen exports.
