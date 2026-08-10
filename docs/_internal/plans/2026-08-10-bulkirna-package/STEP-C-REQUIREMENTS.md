# Step C — requirements collected from the B3–B6 batch

Written during the batch, deliberately as a separate file: `00_INDEX.md` §5a is B6's live
spec and §3 is the shared brief, so editing them while agents run repeats finding 12. Fold
this into `00_INDEX.md` once the batch is closed.

Status of this file: **B4 complete and accepted. B5 code verified, report pending. B6 code
verified with one revision requested. B3 still running.**

---

## 1. Shim requirements that are easy to lose

### `run_gsea()` shim must stash two attributes — from B4

A `gs_result` carries neither the ranked vector nor set membership, but the old S4
`gseaResult` carried both. Rather than change a Step A contract file, `gs_plot_running()`
falls back to `attr(x, "ranks")` and `attr(x, "gene_sets")` when `ranks`/`db` are absent.

**So Step C's `run_gsea()` shim must set both attributes on what it returns.** One line, but
if it is missed the deprecated running-sum path breaks with a confusing "no ranks" error
rather than an obvious failure. Covered by B4's test *"ranks and gene sets can arrive as
attributes of x"*.

### `gsea_running_sum_plot()` — B4's verified recipe

All 14 old formals stay reproducible in order with old defaults. Mapping:

| Old | New | Note |
|---|---|---|
| `gsea_obj` | `x` | positional |
| `gene_set_ids` | `pathways` | both old keyings work: integer row indices and character ids |
| `palette`, `labels` | same | both keyed by pathway id |
| `legend_pos`, `base_size`, `max_name_length`, `title`, `panel_heights`, `base_theme` | identical | same semantics |
| `legend_position` | same | new `c("right","bottom","inside","none")` is a superset of old |
| `es_ylim`, `xticks`, `rug_ylabels` | **absorbed** | accepted and ignored; see below |

`es_ylim` needs per-facet zoom (`ggh4x`, not worth a dependency). `xticks`/`rug_ylabels` are
now structurally true — x ticks only on the bottom panel, rug y index blanked — so both old
values collapse onto the good one. `attr(p, "grs_restyle")` is dead: it existed only because
you cannot `+` onto a patchwork, and the return value is now a plain ggplot.

---

## 2. Golden re-captures — re-capture, never loosen

| Case | Change | Why |
|---|---|---|
| `gsea_running_sum_plot` | class loses its `patchwork` prefix (`ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg`); data goes from 3 sub-plots to 5 layers in one plot | B4's deliberate redesign; enrichplot dropped |

Every other case must stay byte-identical. **A golden change on any B branch before the
harness migration is a bug, not a feature.**

---

## 3. Design-doc corrections owed (all in `sciagent-rna/docs/`)

1. **`08` §4 (B4 brief) and `07` §2** describe the running ES as "a cumulative sum over the
   ranked list". Literally true but it invites the hand-rolled implementation the dispatch
   note forbade. Name `fgsea::plotEnrichmentData()` directly.
2. **`07` §2** undersells the rewrite as "drawing three stacked panels in ggplot". Everything
   except **unequal panel heights** is that easy; panel heights are not expressible in stock
   ggplot2 facets at all, and that is where the module earns its complexity. The doc gives no
   warning.
3. **`07` §7** lists 4 names in the Render group but describes 5 — `gs_save` is missing.
4. **`07` §7's "~30 exports"** predates B6. The real figure is ~35 with the metabolic-network
   layer. Without this, the next reader reads B6 as scope creep.

---

## 4. Open verification items

- **Visual check of one running-sum figure once b3+b4 are both merged.** B4's
  `.grs_default_theme()` resolves `theme_bulki` via `get0(..., asNamespace("bulkiRNA"))` and
  currently falls back to `theme_minimal()` because B3 has not merged. The merge is the first
  time the real foundation theme is exercised. Chrome is applied after the foundation and is
  designed to survive it; confirm rather than assume.
- **`grDevices`** — B4 calls `grDevices::colorRampPalette()` for >9 pathways. Base-bundled and
  unlisted; some `R CMD check` configurations NOTE it. Add to Imports only if
  `check --as-cran` actually complains.
- **`plots[[i]]$coordinates <- shared_coord`** in B5's `de_volcano_grid()` assigns an
  internal field of an S7 ggplot. Works under 4.0.3 and is guarded by an `expect_silent` test
  plus a limits assertion. Confirm the assertion reads the **built** plot's panel ranges, not
  `$coordinates` read back — only the former catches a silent accept-and-ignore under a future
  ggplot2. Awaiting B5's judgement on whether a public-API route was rejected for cause.
- **B6's BUM degeneracy.** The original pipeline test fed 20 all-significant genes, so gatom
  reported `Gene p-value threshold: 1.000000` with the fit pinned to the parameter-space
  boundary — there is no uniform/null component to fit. The plumbing and seed-stability are
  proven; the method's behaviour is not. Realistic-distribution test requested.

---

## 5. Dependency ledger decisions pending agent input

| Package | Status |
|---|---|
| `ggrepel`, `scales`, `stringr` | assigned to B3 — expect the unused-Imports NOTE to shrink |
| `tidyr` | assigned to B5 — confirm actually used |
| `FactoMineR`, `factoextra` | **not in DESCRIPTION**; B5 to request if FAMD stays |
| `readxl`, `org.Hs.eg.db`, `org.Mm.eg.db` | unused Suggests; remove unless B5 claims them |
| `grDevices` | see §4 |
