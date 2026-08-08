# Golden-harness fixture

Synthetic bulk RNA-seq data with **planted gene-set signal**, used to capture golden
outputs before the bulkiRNA refactor and to verify them after. Fully synthetic — no real
biology, no identities.

Regenerate deterministically (`set.seed(20260807)`):

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/pkg -w /pkg \
  scdock-r-dev:v0.5.11 Rscript tests/fixtures/make_fixture.R
```

The `--user` flag is required; without it the container writes as root and `saveRDS` fails
on a permission error.

## Files (272 KB total)

| File | Contents |
|---|---|
| `counts.rds` | integer matrix, 2000 genes × 8 samples, Ensembl-style rownames |
| `metadata.rds` | `sample`, `group` (WT/KO ×4), `batch` (A/B) |
| `toptable.rds` | limma-voom `topTable`, `~ 0 + group + batch`, contrast `KO − WT`; carries `ensembl_id` + `gene_name` |
| `gene_annotation.rds` | Ensembl ID → symbol map |
| `genesets.rds` | `list(alpha = <15 sets>, beta = <8 sets>)`, keyed on symbol |
| `de_intermediates.rds` | filtered counts, voom `logcpm`, `design`, `contrast` |

## What is planted, and why

| Set | Design intent |
|---|---|
| `ALPHA_SIGNAL_UP` (50) | Strong positive. Must come out significant and **up** |
| `ALPHA_SIGNAL_DOWN` (50) | Strong negative. Must come out significant and **down** |
| `ALPHA_NULL_SET` (50) | **Negative control** — a real set with no planted effect. Must stay non-significant |
| `ALPHA_BACKGROUND_01..12` | Random sets, so FDR has something to correct against |
| `BETA_PARTIAL_UP` | 25 genes shared with `SIGNAL_UP` + 15 random — partial overlap across databases |
| `BETA_MIXED` | 12 up + 12 down — **should cancel** and stay non-significant |
| `BETA_BACKGROUND_01..06` | Random |

`beta` deliberately shares genes with `alpha`. Independent databases are the unrealistic
case, and the overlap is what exercises cross-database plotting.

Also planted: **20 duplicated gene symbols** (unique Ensembl IDs) to exercise duplicate-ID
collapsing, and 80 scattered non-set DE genes so the volcano is not two clean blobs.

## Verified behaviour at `752481f` / `scdock-r-dev:v0.5.11`

1872 of 2000 genes survive `filterByExpr`; 27 DE at `adj.P.Val < 0.05` (16 up / 11 down).

| Engine | `SIGNAL_UP` | `SIGNAL_DOWN` | `NULL_SET` | `BETA_MIXED` |
|---|---|---|---|---|
| fgsea | NES +3.32, padj 9e-30 | NES −2.97, padj 1e-15 | padj 0.66 | padj 0.21 |
| fora | FE 7.56, padj 7e-31 | FE 5.88, padj 2e-15 | ns | — |
| GSVA → limma | t +10.9, padj 2e-10 | t −10.5, padj 3e-10 | ns | — |

All three engines agree on direction and significance while reporting different
statistics — which is what the `stat` / `stat_type` contract exists to carry.

Per-gene DE is deliberately modest (27 genes). GSEA is more sensitive than per-gene
testing, and a fixture where everything is significant would not exercise the
non-significant rendering paths.
