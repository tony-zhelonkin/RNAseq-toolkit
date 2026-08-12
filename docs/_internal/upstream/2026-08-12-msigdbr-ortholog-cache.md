# Upstream bug report — draft

**Where:** https://github.com/igordot/msigdbr/issues (the `BugReports` field in
`packageDescription("msigdbr")`). Maintainer: Igor Dolgalev.

**How:** a single GitHub issue. No email, no CRAN contact — the package is actively developed
on GitHub and this is a plain reproducible bug. Include the reprex; it is short and it does the
arguing for us.

---

## Title

Ortholog cache is keyed on species only, so every collection after the first in a session is
silently reduced to its intersection with the first

---

## Body

`msigdbr()` caches the ortholog mapping table under a key derived from the target species
alone, but builds that table from **only the genes of the collection being queried**. Every
later call in the same session reuses it, so the result is silently reduced to the intersection
of the two collections' gene spaces. No warning is emitted, and the returned object looks
entirely normal.

Because the effect depends on call *order* rather than call *arguments*, a script gives
different answers depending on which collection it happens to query first.

### Reproducible example

Two fresh R sessions, identical final call:

```r
# Session 1 --- query Reactome alone
library(msigdbr)
r <- msigdbr(db_species = "HS", species = "mouse",
             collection = "C2", subcollection = "CP:REACTOME")
length(unique(r$gene_symbol))
#> 10762
nrow(r)
#> 105806

# Session 2 --- query Hallmark first, then the exact same Reactome call
library(msigdbr)
invisible(msigdbr(db_species = "HS", species = "mouse", collection = "H"))
r <- msigdbr(db_species = "HS", species = "mouse",
             collection = "C2", subcollection = "CP:REACTOME")
length(unique(r$gene_symbol))
#> 3688          # <-- 66% of the mapped genes are gone
nrow(r)
#> 44635
```

### The effect is symmetric — whichever collection runs second is truncated

```r
# Hallmark first
H = 4393   then Reactome = 3688     # Reactome loses 7074 genes
# Reactome first
Reactome = 10762   then H = 3688    # Hallmark loses 705 genes
```

`3688` in both directions is the intersection of the two collections' mapped gene spaces. The
first call is always correct; the second is always reduced to that intersection.

### Cause

In `msigdbr()`:

```r
species_id <- babelgene::species(species)$taxon_id
orthologs_key <- paste0("orthologs", species_id)
if (exists(orthologs_key, envir = pkg_env, inherits = FALSE)) {
  species_genes <- pkg_env[[orthologs_key]]
} else {
  species_genes <- babelgene::orthologs(genes = unique(mdb$db_ensembl_gene), species = species)
  ...
  pkg_env[[orthologs_key]] <- species_genes
}
```

`mdb` has already been filtered by `collection` — `load_gene_sets(collection =)` reads only
that collection's RDS files — and by `subcollection` immediately above. So
`unique(mdb$db_ensembl_gene)` is one collection's gene space, while `orthologs_key` records only
the taxon. The subsequent
`inner_join(mdb, species_genes, by = "db_ensembl_gene")` then drops every gene absent from the
cached table.

### Why this is easy to miss

Nothing fails. The call returns a well-formed tibble with plausible set sizes, so the only
visible symptom is downstream: gene sets test as smaller than they are, and in GSEA the
resulting p-values and FDR are inflated. We found it only by noticing that one script reported
677 Reactome pathways where an equivalent one reported 1,050, and the discrepancy tracked the
order of collections in a loop.

### Suggested fix

The cache key is not really the problem — the *scope of the cached table* is. Building it from
the full human gene space of the database, rather than from the queried subset, would make a
species-only key correct:

```r
species_genes <- babelgene::orthologs(genes = <all db_ensembl_gene for the species>, ...)
```

Then the per-call `inner_join` subsets it naturally, and the cache is reusable by construction.
That costs one larger `babelgene::orthologs()` call on the first query of a session and removes
the order dependence entirely.

Alternatives, if that is unattractive: include the queried gene set in the cache key, or drop
the memoisation.

### Workaround, for anyone hitting this before a fix

```r
key <- paste0("orthologs", babelgene::species("mouse")$taxon_id)
if (exists(key, envir = asNamespace("msigdbr")$pkg_env, inherits = FALSE)) {
  rm(list = key, envir = asNamespace("msigdbr")$pkg_env)
}
```

Dropping the cache before each query restores correct results (`H = 4393`, then
`Reactome = 10762`). It reaches into the package's internals, so it is a stopgap rather than a
fix — and it is worth pairing with an independent check on the returned gene coverage, since a
workaround that silently stops working restores exactly the silent failure it was meant to
prevent.

### Session info

```
msigdbr   26.1.0
babelgene 22.9
dplyr     1.2.1
R version 4.5.3 (2026-03-11)
```

Happy to open a PR for the fix if the approach above looks right to you.
