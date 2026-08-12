# Architecture decisions of record

Dated, numbered, and **durable**: unlike `plans/`, which records what we were about to do,
this directory records what we decided and why, so a future reader can tell a deliberate
choice from an accident.

One file per decision. A decision that is later reversed is not deleted — it is marked
`Superseded by ADR-NNN` and left in place, because the reasoning is the valuable part.

| ADR | Decision | Status |
|---|---|---|
| [ADR-001](ADR-001-reproducibility-unit.md) | The package version is the unit of reproducibility | Accepted 2026-08-12 |
| [ADR-002](ADR-002-master-table-schema.md) | `bulkiRNA` owns the master-table schema and validator, not the file | Accepted 2026-08-12 |
| [ADR-003](ADR-003-optional-dependencies.md) | `Suggests` stays; add a bulk preflight and `Remotes` | Accepted 2026-08-12 |
| [ADR-004](ADR-004-reference-data-tiers.md) | Two reference-data tiers: bundled or refcache, nothing else | Accepted 2026-08-12 |

Earlier decisions predate this directory and live in the plan index
(`plans/2026-08-10-bulkirna-package/00_INDEX.md`) — notably MADR-008 (`genes_full_set` is
computed by the caller), the four-layer architecture, and the export freeze. They are not
back-filled here; the plan index remains their home.

## Upstream

[`../upstream/`](../upstream/) holds bug reports filed against, or drafted for, our
dependencies — with the reproducible example alongside the prose, so the claim can be re-checked
rather than taken on trust.

| Report | Status |
|---|---|
| [msigdbr ortholog cache](../upstream/2026-08-12-msigdbr-ortholog-cache.md) | drafted 2026-08-12, not yet filed |
