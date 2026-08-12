# ADR-003 — `Suggests` stays; add a bulk preflight and `Remotes`

**Status:** Accepted 2026-08-12 · **Decided by:** Anton Zhelonkin
**Consulted:** codex `gpt-5.5`, high reasoning, 2026-08-12 — concurred independently

## Context

`DESCRIPTION` carries 14 light `Imports` (`dplyr`, `fgsea`, `ggplot2`, `ggrepel`,
`msigdbr`, `rlang`, `scales`, `stringr`, `tibble`, plus base-adjacent) and 20 `Suggests`.
The concern raised was that `Suggests` looks like "redundancy for no reason" given that
every optional dependency is already guarded at its use site, and that the primary
distribution is a container where the heavy dependencies are present anyway.

The guards are real and already correct — `.require_pkg(pkg, what, install =)` in
`R/utils.R:66`, used at **14 sites**, with source-accurate install hints:

```r
.require_pkg("edgeR", "PCA of a `DGEList`", 'BiocManager::install("edgeR")')
.require_pkg("gatom", "gatom_module()",     'BiocManager::install("gatom")')
.require_pkg("mwcsr", "gatom_module()")     # CRAN — the default hint is right
```

## Decision

**Keep `Suggests`. No packages move.** Add two things:

1. **An exported preflight**, `bulkirna_check_deps(features = ...)`, reporting every
   optional package's presence, version, and correct install command, grouped by feature
   area.
2. **`Remotes:`** for dependencies that a plain CRAN mirror cannot resolve.

And run `bulkirna_check_deps("all")` as a **post-build smoke test in `scbio-docker`**.

## Rationale

`Suggests` is not a second copy of the guards — it is *metadata about* the guards. The
guards make an absent package fail legibly at call time; `Suggests` is the machine-readable
list that `R CMD check`, `devtools::install_deps(dependencies = TRUE)`, `pak`, and the
image build read. Deleting it keeps the late failure and removes the ability to enumerate
what is missing: the worst of both.

The decisive constraint is the install surface. `gatom`, `mwcsr`, `GSVA` and `org.*.eg.db`
do not resolve from a plain CRAN mirror. Today
`remotes::install_github("tony-zhelonkin/bulkiRNA@v0.3.1")` succeeds on bare R in seconds
**because** the hard floor is 14 light packages. Promoting the optional 20 to `Imports`
turns every optional feature into a hard install blocker — including for users who never
touch those features. That lean floor is the asset; `Suggests` is what protects it.

R's own guidance matches: packages used only in examples, tests, vignettes, or optional
function bodies belong in `Suggests`, consumed conditionally via `requireNamespace()` plus
`pkg::fun()`. That is exactly the implemented pattern.

## What was actually missing

Only one thing: **nothing reports the state of the whole optional closure at once.** The
five packages absent from the image (`org.Mm.eg.db`, `org.Hs.eg.db`, `homologene`, plus
`gatom`/`mwcsr`) were found by grep, not by a tool. The preflight closes that; the guards
remain the correctness mechanism, and the preflight is convenience.

## Consequences

- Examples, tests and vignettes touching `Suggests` must use
  `testthat::skip_if_not_installed()` / conditional evaluation.
- `_R_CHECK_DEPENDS_ONLY_=true` is the way to *prove* the core installs without the
  optional closure. Worth adding to the check routine.
- Suggested-package **version** requirements are not enforced at `library()` time; where a
  version matters, `.require_pkg()` must check it.
- Anything optional used unconditionally at load time, in an S4 class definition, or in
  namespace registration belongs in `Imports` — none currently are.

## Rejected alternatives

- **Collapse everything into `Imports`.** Rejected: hard install blocker, per above.
- **Split into core + feature packages** (`bulkiRNA.gsea`, `bulkiRNA.network`, …).
  Rejected: cleanest isolation, but it buys a second release train and API surface for a
  ~60-export internal package with one maintainer and one audience.
- **`Additional_repositories` instead of `Remotes`.** Not either/or — `Remotes` is honoured
  by `remotes`/`devtools` and not by base `install.packages()`, so repository metadata is
  complementary, not a substitute.
