# Upstream bug report

**Where:** https://github.com/igordot/msigdbr/issues (the `BugReports` field in
`packageDescription("msigdbr")`). Maintainer: Igor Dolgalev.

**How:** a single GitHub issue. No email, no CRAN contact. The package is actively developed on
GitHub and this is a plain reproducible bug. Include the reprex; it does the arguing.

**Filed:** see the `## Filed` section at the end.

---

## Title

Ortholog cache is keyed on species only, so every collection after the first in a session is
silently reduced to its intersection with the first

---

## Body

`msigdbr()` caches the ortholog mapping table under a key derived from the target species alone,
but builds that table from only the genes of the collection being queried.

Every later call in the same session reuses it. The result is silently reduced to the
intersection of the two collections' gene spaces. No warning is emitted, and the returned object
looks normal in principle, except that the number of sets tested is smaller than you'd expect.

Because the effect depends on call order rather than call arguments, the same script gives
different answers depending on which collection it queries first.

### Reproducible example

Two fresh R sessions, identical final call:

```r
# Session 1 - query Reactome alone
library(msigdbr)
r <- msigdbr(db_species = "HS", species = "mouse",
             collection = "C2", subcollection = "CP:REACTOME")
length(unique(r$gene_symbol))
#> 10762
nrow(r)
#> 105806

# Session 2 - query Hallmark first, then the exact same Reactome call
library(msigdbr)
invisible(msigdbr(db_species = "HS", species = "mouse", collection = "H"))
r <- msigdbr(db_species = "HS", species = "mouse",
             collection = "C2", subcollection = "CP:REACTOME")
length(unique(r$gene_symbol))
#> 3688          # 66% of the mapped genes are gone
nrow(r)
#> 44635
```

### Whichever collection is queried first, the second gets truncated

```r
# Hallmark first
H = 4393   then Reactome = 3688     # Reactome loses 7074 genes
# Reactome first
Reactome = 10762   then H = 3688    # Hallmark loses 705 genes
```

The second call is reduced to the intersection of the two gene spaces.
That is exact at the join key. Measuring `db_ensembl_gene` from clean sessions:

```r
H = 4362   Reactome = 10847   intersection = 3655
Reactome queried after Hallmark: 3655
```

The mouse-symbol count above (3688) differs from a naive symbol-level intersection (3689) by one,
because ortholog mapping is many-to-many. The truncation itself happens on the join key.

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

`mdb` has already been filtered by `collection`, since `load_gene_sets(collection =)` reads only
that collection's RDS files, and by `subcollection` immediately above.

So `unique(mdb$db_ensembl_gene)` is one collection's gene space, while `orthologs_key` records
only the taxon.

The subsequent `inner_join(mdb, species_genes, by = "db_ensembl_gene")` then drops every gene
absent from the cached table.

### Why this is easy to miss

Nothing fails. The call returns a well-formed tibble with plausible set sizes, so the only
symptom is downstream, when gene sets test as smaller than you'd expect. The resulting p-values
and FDRs are inflated too.

I found it by noticing that one script reported 677 Reactome pathways where an equivalent one
reported 1,050.

### Suggested fix

Building the table from the full human gene space of the database would make a species-only key
correct:

```r
species_genes <- babelgene::orthologs(genes = <all db_ensembl_gene for the species>, ...)
```

Then the per-call `inner_join` subsets it correctly. On the first query of a session that costs
one larger `babelgene::orthologs()` call, and the order dependence goes away.

If that is unattractive: include the queried genes in the cache key, or drop the memoisation.

### Workaround

```r
key <- paste0("orthologs", babelgene::species("mouse")$taxon_id)
if (exists(key, envir = asNamespace("msigdbr")$pkg_env, inherits = FALSE)) {
  rm(list = key, envir = asNamespace("msigdbr")$pkg_env)
}
```

Dropping the cache before each query restores correct results.

It reaches into package internals, so it needs an independent check on the returned gene
coverage. A workaround that silently stops working restores exactly the failure it was meant to
prevent.

### Session info

```
msigdbr   26.1.0
babelgene 22.9
dplyr     1.2.1
R version 4.5.3 (2026-03-11)
```

Happy to open a PR if the approach above looks right to you.

---

## Filed

Opened with `gh issue create --repo igordot/msigdbr`. URL recorded below once filed.
