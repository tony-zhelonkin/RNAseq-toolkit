# What CoReSh's method actually is, and where it lives

**Date:** 2026-08-13 · **Status:** durable reference; every claim below is measured or quoted
**Parent:** [00_ROADMAP.md](00_ROADMAP.md)

The question this answers: CoReSh advertises "a comparison-independent score for ranking
datasets from the compendium by their relevance to the query gene set. Inspired by PCA, this
score captures the dataset variance explained by a certain direction in the sample space."
Does the source implement that method, and can we run the analysis they intended?

**Short answers.** The method is real and fully implemented — **in `fgsea`, not in `coresh`**.
The `coresh` package implements nothing. The score is `fgsea`'s GESECA score, line for line.
We can run exactly the analysis the authors intended, and one of the two rankings they
document already works and is verified bit-identical.

---

## 1. The score, traced to its implementation

**The vignette says so itself.** `vignettes/coresh-local.Rmd`, immediately before defining
`coreshMatch()`:

> Let's define a function to match an object against the query gene set, which will
> rely on `geseca` method from the `fgsea` package

The wrapper computes:

```r
curProfile <- colSums(E[queryIdxs, , drop = FALSE])
queryVar   <- sum(curProfile ** 2)
pctVar     <- queryVar / k / obj$totalVar * 100
```

And `fgsea`'s own scoring function, in full:

```r
fgsea:::calcGesecaScores <- function(indxs, E) {
  return(sum(colSums(E[indxs, , drop = FALSE])^2))
}
```

**These are the same expression.** Verified numerically: wrapper `317.9887245020`, fgsea
`317.9887245020`, ratio exactly 1.

The normalization is also fgsea's. From `fgsea::gesecaSimple`'s body:

```r
totalVar <- sum(rowSums(E^2))
pvals[, pctVar := pctVar / size / totalVar * 100]
```

`queryVar / k / totalVar * 100` is `pctVar / size / totalVar * 100`. **`coreshMatch()` is
`fgsea::geseca()` restricted to one gene set, with the denominator read from storage instead
of recomputed.** Nothing else.

The p-value is fgsea's too: `fgsea:::gesecaCpp` in the vignette, which is the internal helper
behind the public `fgsea::geseca()` multilevel estimator. `src/geseca.cpp`,
`R/geseca-multilevel.R`, `R/geseca-simple.R`, `R/geseca-utils.R` and `R/geseca-plot.R` are all
fgsea files.

### So what does the `coresh` package contain?

`DESCRIPTION`, `NAMESPACE` (`exportPattern("^[[:alpha:]]+")` with nothing to match),
`README.md`, two licence files, and one vignette. **No `R/` directory.** Three commits, all
April 2025. It is a dependency declaration plus a documented recipe — installing it gives you
`fgsea`, `qs2`, `data.table`, `fastmatch`, `msigdbr`, `dplyr` and `BiocParallel`, and a
vignette to copy from. That is a legitimate way to publish a method built on an existing
implementation; it is just not a library, and nothing can call into it.

---

## 2. Why "inspired by PCA" is accurate, and what the score really measures

Let `E` be the stored matrix, rows genes, columns the sample-space coordinates, rows centred.
For a query set `S` of size `k` present in the dataset, let

- `p = Σ_{g∈S} e_g` — the summed profile, a vector in sample space
- `score = ‖p‖²`, and `pctVar = score / k / totalVar × 100`

**It is a variance along a direction, up to a bound.** Take the set's own mean direction
`u = p/‖p‖`. By Cauchy–Schwarz,

```
score / k  =  ‖p‖² / k  ≤  Σ_{g∈S} (e_g · u)²
```

with equality exactly when every gene projects equally onto `u`. Verified: `84.406` against
`102.085` on a correlated fixture. So `score/k` is a **closed-form lower bound on the variance
the query genes carry along their own mean direction** — the same quantity PCA would obtain
from a leading eigenvector, without an eigendecomposition. "Inspired by PCA" rather than "is
PCA" is the honest description, and the paper's wording is fair.

**What it detects is coregulation, and the scale is interpretable.** Expanding,

```
score / k  =  mean_g ‖e_g‖²  +  (1/k) Σ_{g≠h} ⟨e_g, e_h⟩
```

so the score is the average per-gene variance *plus* the cross-gene covariance. Divided by the
mean per-gene variance it is a coregulation multiplier bounded by `k`. Measured on a latent
factor model with `k = 20`:

| latent loading ρ | `score/k` ÷ mean gene variance |
|---|---|
| 0.0 | 1.14 |
| 0.3 | 2.02 |
| 0.6 | 8.57 |
| 0.9 | 16.21 |
| 0.999 | 19.93 |

1 for independent genes, `k` for perfectly coregulated ones. That is the whole idea.

**"Comparison-independent" is the real contribution.** The score is a property of one
`(dataset, gene set)` pair. There is no contrast, no case and control, no ranking of genes,
and nothing to choose — which is what lets one query sweep 40,000 datasets whose designs are
unknown and mutually incomparable. Dividing by `k` removes the query-size dependence and
dividing by `totalVar` removes the dataset's depth and scale, leaving a number comparable
across the compendium. **The innovation is the normalization plus the compendium, not a new
statistic.**

---

## 3. The compendium, which is the part that is genuinely theirs

Per the vignette: `E1024` "is centered, potentially reduced with a Principal Component
Analysis, multiplied by 1024 and rounded to the nearest integer."

Measured on `chunk_1` of the mouse tree (500 datasets, snapshot `syn66227307_20260721`):

| Property | Value |
|---|---|
| Datasets per chunk file | ~500 |
| Chunk files, `hsa` + `mmu` | 174 |
| Genes per matrix | up to 10,000 |
| Sample counts | 4 to 437 |
| **Stored matrix columns** | **4 to 20 — capped at 20** |
| **PCA-reduced (`ncol(E) < nsamples`)** | **96 of 500** |
| Fields per object | `gseId`, `gplId`, `E1024`, `rownames`, `samples`, `nsamples`, `totalVar`, `wordMatrix` |

GSE10415 has 437 samples stored as 20 columns. So for a fifth of the compendium **the columns
are principal components, not samples** — which is why the marketing phrase is "a direction in
the sample space" rather than "a sample".

**Why the reduction is lossless for this statistic, and this is the elegant part.** Both the
numerator and the denominator are invariant under an orthogonal change of basis on the
columns. If `E' = EQ` with `Q` orthogonal, then `colSums(E'[S,]) = colSums(E[S,])·Q`, whose
squared norm is unchanged, and `sum(E'²) = sum(E²)` because the Frobenius norm is invariant.
Verified: score `317.9887245020` before and after a random orthogonal rotation, total variance
`4498.56616747` before and after, both bit-identical.

So a PCA rotation costs nothing and truncation costs only the variance in the dropped
components. Storing 20 columns instead of 437 is a 20× saving on a statistic that cannot tell
the difference. **`totalVar` is stored because it describes the matrix as stored**: measured
against `sum((E1024/1024)²)` the relative gap is ~1e-6 for reduced matrices and ~3e-6 median
for unreduced ones, maximum 6.1e-5 — so the gap is the 1024-step quantization, **not** the
truncation. My earlier note attributing it to quantization alone was right, for a reason I had
not checked at the time.

---

## 4. Can we run the analysis they intended? Yes

The vignette documents **two rankings**, and treats both as first-class:

```r
varRanking  <- varRanking[order(pctVar, decreasing = TRUE)]   # 10-20 seconds, whole compendium
pvalRanking <- pvalRanking[order(pval)]                        # a couple of minutes
```

> Ranking by p-value usually is a bit more specific, but takes about couple of minutes to
> calculate.

| Their step | Our status |
|---|---|
| `pctVar` per dataset | ✅ `coresh_match()` / `coresh_search()`, **bit-identical to hand computation on 12 of 12 real datasets**, difference exactly 0 |
| Chunk iteration with `BiocParallel` over files | ✅ `coresh_search(n_cores = )` |
| Rank by `pctVar` descending | ✅ |
| GESECA p-value per dataset | ⏳ **G1** — route to `fgsea::geseca()`; agrees with `gesecaCpp` within summed `log2err` in 7 of 8 datasets |
| Rank by `pval` ascending | ⏳ **G1** — `coresh_search()` currently ranks by `pct_var` unconditionally and must restore this branch |
| Query from MSigDB `ncbi_gene` as integer | ✅ `gsdb_msigdb()` + `gene_to_entrez()`; the vignette's own example is `HALLMARK_HYPOXIA` |
| Dataset metadata | ⚠️ **not in the snapshot.** The vignette fetches titles from GEO per accession. Ours should too, or say plainly that it does not. |

**Two fidelity items G1 must fix, both introduced by me and both traceable to a premise that
turned out to be false:**

1. `pvalues = TRUE` currently errors. Upstream treats the p-value ranking as the *more
   specific* of the two, so it is not an optional extra.
2. `coresh_search()` dropped the p-value ordering branch. That was correct only while the
   p-value did not exist.

**The one thing we should keep doing differently.** Use the stored `obj$totalVar`, not
`geseca()`'s recomputed `sum(E²)`. It was computed before quantization, it is what the authors
divide by, and it costs nothing to prefer. Where `geseca()` computes its own `pctVar` we take
its p-value and keep ours.

**Reference smoke test, in this order:** the web UI at <https://alserglab.wustl.edu/coresh>
first, then `HALLMARK_HYPOXIA` locally through `coresh_search()`, then compare the top
accessions. Two independent implementations of one method agreeing on the top hits is worth
more than any unit test in the suite.
