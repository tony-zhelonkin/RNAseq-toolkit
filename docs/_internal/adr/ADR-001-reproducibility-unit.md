# ADR-001 — The package version is the unit of reproducibility

**Status:** Accepted 2026-08-12 · **Decided by:** Anton Zhelonkin

## Context

Four version authorities were floating free: the git tag, `DESCRIPTION: Version`, the
`scbio-docker` image `VERSION`, and whatever `write_session_provenance()` happened to
record. The toolkit's own ADR-D1 says "the submodule is the scientific lock" — but
`bulkiRNA` is deliberately no longer a submodule, so that clause stopped covering it and
nothing replaced it. A result could not name the code that produced it without ambiguity
about *which* of the four numbers was the claim.

## Decision

**The `bulkiRNA` package version is the unit of reproducibility.**

An `scbio-docker` tag is a promise about which tool versions the image contains, and
nothing more. Provisioning the image is `scbio-docker`'s own responsibility, on its own
clock. Analysis repositories carry their own version lock for anything that may have
drifted when a runtime install happened.

## Consequences

Two mechanical obligations follow, and they are the whole cost of this decision:

1. **`Version` must not move without a tag.** The install pin is
   `remotes::install_github("tony-zhelonkin/bulkiRNA@v0.3.1")`. A commit that bumps
   `DESCRIPTION: Version` without creating the matching tag makes that pin resolve to code
   that does not match its own version string — the one failure mode that would make the
   unit meaningless. A release check should refuse the bump unless the tag exists.

2. **`write_session_provenance()` is what analysis repos lock against**, so it must emit
   the full package-version closure, not only `bulkiRNA`. The drift being guarded against
   is `fgsea` or `msigdbr` moving under a fixed image tag.

## Why this is the right cut

The `msigdbr` 26.1.0 ortholog-cache bug (plan index §12) would have been **invisible**
under an image-as-unit model: the image tag never changed, and the numbers did. Under
package-as-unit with a dependency closure in provenance, a re-run whose numbers move has
something concrete to point at.

## Rejected alternatives

- **The image tag as the unit.** Rejected: it hides dependency drift, as above, and it
  couples the scientific claim to an infrastructure artifact rebuilt for unrelated reasons.
- **A lock file inside `bulkiRNA`.** Rejected: a library cannot lock its own consumers'
  environments, and `renv.lock` is never `COPY`'d into the image anyway.
