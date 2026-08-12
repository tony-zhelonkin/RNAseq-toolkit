# Legacy pre-package documentation

These four documents describe the **pre-package** RNAseq-toolkit: `source()`d scripts,
`clusterProfiler::GSEA()`, `gseaResult` S4 objects, and the master-table conventions of the
12868-EH project specifically.

They are kept for provenance — they explain why several design decisions look the way they
do, and they are the written record of the pipeline the golden baselines were captured
against. **None of their code examples run against bulkiRNA.**

For current documentation see [../../../README.md](../../../README.md),
[../../API_REFERENCE.md](../../API_REFERENCE.md), [../../WORKFLOWS.md](../../WORKFLOWS.md),
and [../../../MIGRATION.md](../../../MIGRATION.md) for the old → new mapping.

| File | Describes |
|---|---|
| `01-core-pipeline-and-toolkit.md` | counts → master CSV flow, config system, checkpoint caching |
| `02-msigdb-gsea-pipeline.md` | the clusterProfiler MSigDB path and its normalized schema |
| `03-custom-database-gsea.md` | T2G/T2N model, external database parsing, append pattern |
| `04-output-artifacts-and-visualization.md` | output layout, CSV schemas, pathway-explorer inputs |
