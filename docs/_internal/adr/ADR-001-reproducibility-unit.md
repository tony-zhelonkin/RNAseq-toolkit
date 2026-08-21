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

## Amendment 2026-08-20 — the inverse obligation

Obligation 1 is one-directional, and that was the defect. It forbids a version moving without
a tag. It never forbade **code moving while still claiming an already-released version.**

Measured: `DESCRIPTION` read `Version: 0.5.0` while `HEAD` stood 50 commits and ~3,032 changed
lines of `R/` past the `v0.5.0` tag, adding nine files and growing the surface from 64 to 79
exports. Two materially different code states both called themselves `0.5.0`, so the unit
identified nothing and every downstream record inherited the ambiguity. `NEWS.md` already read
`0.5.0.9000`, so the discipline existed in prose rather than as a gate.

3. **A release version means exactly one commit.** If `Version` is `A.B.C` then `HEAD` is the
   commit `vA.B.C` points to and the tree is clean. Any commit after a release carries
   `A.B.C.9000` with `A.B.C` strictly above every released tag, which makes the next release
   number explicit and a stale version detectable. A tag `vA.B.C` points at a commit whose
   `DESCRIPTION` says `A.B.C`. Enforced by a test that fails rather than a release checklist:
   this accumulated over 50 commits, which is past any checklist's reach.

4. **Provenance records the version and the source commit**, not one or the other. For a
   release the commit witnesses the version-to-tag mapping; for a `.9000` build it is the only
   thing identifying the code, since every development commit shares one version string.

**The general lesson, worth more than the fix.** Every *content* invariant here is enforced by a
test that fails: the master-table schema, RNG ownership, name spellings, return types, the
registry against `NAMESPACE` in both directions. Identity and deployment invariants were left to
prose, because they live between repositories, tags, builds and installed libraries, where
package tests do not reach. That asymmetry is where to look for the next instance.

One such instance is visible in this document. The rejected-alternatives note below observes that
`renv.lock` is never `COPY`'d into the image — recorded in passing to support a different
decision, never treated as a defect. It is the mechanism by which seven packages silently
vanished between two image builds.

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
