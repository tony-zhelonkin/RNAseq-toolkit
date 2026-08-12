# ADR-004 — Two reference-data tiers: bundled or refcache, nothing else

**Status:** Accepted 2026-08-12 · **Decided by:** Anton Zhelonkin

## Context

Three different mechanisms had grown for "reference data the package needs":

1. **Bundled** in `inst/extdata/` — MitoCarta 3.0, mitoXplorer 3.0, transportdb,
   mitochondria_unified. **168 KB total**, with `inst/extdata/METADATA.yaml` already
   recording `name`, `description`, `source_url`, `license`, `version`, `bundled`, and
   `citations` per entry.
2. **refcache** — the shared snapshot cache (`tony-zhelonkin/refcache`): a generic
   `refcache.sh` driver (snapshot → verify → `MANIFEST.json` → atomic `rename(2)` flip →
   prune), one file per source, an ephemeral fetcher image holding credentials, and 115 GB
   of bytes at a deployment path exposed to containers as `REFCACHE_ROOT=/refcache` (`:ro`).
   Two live sources: cisTarget 73 G, CoReSh 42 G.
3. **Ad-hoc download** — `download_gatom_references()`, ~24 MB over public HTTPS, with a
   hardcoded `dest_dir = "00_data/references/gatom"` default.

Left alone, `gsdb_coresh()` would have invented a fourth via `Sys.getenv("CORESH_CHUNKS")`.

## Decision

**Two tiers only.**

| Tier | Criterion | Mechanism | Examples |
|---|---|---|---|
| **Bundled** | no stable API or fetchable versioned source; small | `inst/extdata/` + `METADATA.yaml`, versioned by the package version | MitoCarta, mitoXplorer, transportdb — 168 KB |
| **refcache** | *has* a fetchable versioned source; shared across projects | `$REFCACHE_ROOT/<source>/current` | CoReSh 42 G, cisTarget 73 G |

> **The rule: if it can be re-fetched from a versioned source, it is refcache; if it cannot,
> it is bundled. Size is the tiebreaker, not the criterion.**

Three further consequences:

**Bundling MitoCarta was correct and stays.** A source with no API cannot be re-fetched
reproducibly, so there is nothing for refcache to lock — **vendoring is the archive**. At
168 KB the size objection does not exist, and `METADATA.yaml` is already the provenance
registry that makes it defensible.

**Tier 3 is retired.** GATOM's references match the refcache source profile exactly —
public HTTPS, versioned files, sha-verifiable — the same profile as `sources/cistarget.sh`.
The driver's design claim is that adding a source is one file, so: write
`sources/gatom.sh`, and `download_gatom_references()` becomes a deprecated shim. This also
removes the last hardcoded host-path default from the package.

**One internal resolver, not per-provider environment variables.** `.ref_path(source, ...)`
resolves in a fixed order — explicit argument → `$REFCACHE_ROOT/<source>/current` → error
naming the exact `refcache.sh` command — and no provider ever knows a host layout.
`Sys.getenv("CORESH_CHUNKS")` in the CoReSh plan §4 is superseded by it.

## The pinning gap, and how this closes it

The refcache design note states the tension plainly: *"Always resolve `current/` trades
reproducibility for freshness. A refresh silently changes what a re-run of last month's
analysis reads… pinning an analysis to `cistarget_20260811` is a project-level concern that
nothing currently owns."*

With ADR-001 settled, **nobody needs to own pinning; someone needs to own recording.**
`current` is a symlink, so the snapshot tag is knowable at resolve time.

> **`.ref_path()` must resolve `current` and report the tag it landed on, and
> `write_session_provenance()` must record it.**

A re-run that reads different bytes then produces a lock differing from last month's, and
the change becomes *visible* instead of silent. Freshness stays the default — which is the
reason the cache was built this way — and reproducibility becomes detectable after the
fact. Pinning to a specific snapshot remains available through the explicit argument, for
the rare analysis that needs it.

This is the `msigdbr` guard's lesson one layer out: **do not prevent the drift; make it
impossible for the drift to be silent.**

## Residual items

- Neither live `MANIFEST.json` carries `produced_by`, so today's 115 GB cannot name the
  driver revision that produced it. Accepted debt — recording the snapshot *directory name*,
  which does exist, recovers most of the traceability now.
- `license: "Academic (Broad Institute)"` inside an MIT-licensed package needs one `NOTICE`
  paragraph clarifying that the MIT grant covers the code, not the bundled datasets. Not
  urgent; required before anything is published.
- `inst/extdata/METADATA.yaml` carries `version: "1.0.0"` for the registry, which is not yet
  surfaced in result provenance. Fold into the same provenance change.

## Rejected alternatives

- **A refcache source for MitoCarta.** Rejected: no stable URL or API means the fetch is not
  reproducible, so this is strictly worse than vendoring.
- **A separate `bulkiRNA.data` package.** Rejected: the standard Bioconductor pattern, and
  it would decouple data updates from code releases — but at 168 KB it buys a second release
  train for nothing. Revisit only if bundled data passes single-digit MB.
- **Per-provider environment variables** (`CORESH_CHUNKS`, `CISTARGET_DIR`, …). Rejected:
  one code path per source is how three mechanisms became three, and each invents its own
  error message.
