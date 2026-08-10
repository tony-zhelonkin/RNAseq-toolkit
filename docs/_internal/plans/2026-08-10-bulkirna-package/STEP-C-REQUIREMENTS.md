# Step C — requirements collected from the B3–B6 batch

Written during the batch, deliberately as a separate file: `00_INDEX.md` §5a is B6's live
spec and §3 is the shared brief, so editing them while agents run repeats finding 12. Fold
this into `00_INDEX.md` once the batch is closed.

Status: **B3, B4, B5, B6 all complete and merged.** 41 exports · `check --as-cran` 0/0/0 ·
609 tests (dev) / 635 (gatom image), 0 failures · golden 20/20. Reports received from B4 and
B5 in full; **B3 and B6 never reported** — their reasoning was reconstructed by the integrator
from the source, so treat §6 below as the record.

---

## 0. Frozen-name shim recipes — CONFIRMED

### From B5 (verbatim, verified against the source)

`create_standard_volcano` → `de_volcano` and `create_MD_plot` → `de_md_plot` are **1:1 with no
argument renaming**, so both shims are literal pass-throughs. Every formal keeps its name,
order and default. `build_dge` and `ensure_dir` **keep their own names** — no shim needed, and
B5 confirmed it never redefined `ensure_dir`.

The design decision that protects this: `de_volcano`'s one new argument, `orientation`, sits
**after `...`**, so R can only match it by name and no positional call against the frozen
signature can shift onto it. (`orientation = "vertical"` is where `create_vertical_volcano`
went; `combine_volcano_row` became `de_volcano_grid`.)

### From B3's source (reconstructed — B3 did not report)

- `theme_bulki(base_size = 14, base_family = "", grid = FALSE)`; `custom_minimal_theme_with_grid`
  → `theme_bulki(grid = TRUE)`.
- `plot_all_gsea_results` body → `R/gs-plot-all.R`; `save_gsea_log` body → `R/gs-save.R`.
  Both files carry a roxygen line naming the deprecated function they back.
- `format_pathway_name` keeps its own name and stays exported — **the golden harness self-test
  wraps precisely this function in `toupper()`**, so its signature must not move.
- `gsea_dotplot`/`gsea_dotplot_facet`/`gsea_barplot` → `gs_plot_dot`/`gs_plot_bar`. Argument
  mapping NOT confirmed by an agent report; **Step C must diff the formals itself.**

---

## 0a. The `--as-cran` clean-up already applied at integration

| Fix | Commit | Why it only appeared on merge |
|---|---|---|
| `patchwork` → Suggests | `4ac6539` | `de_volcano_grid()` used it via `::` behind a guard; declared nowhere |
| `tidyr` out of Imports | `4ac6539` | assigned to B5, honestly unused — B5 declined to fake a use |
| `grDevices`, `tools` → Imports | `34000b7` | used via `::` in three modules, never declared |
| ASCII-only R sources | after `0fe0cd4` | B5 cleared its 7 files; 14 lines remained in Step A/B1/B2 files |

**Non-ASCII rule for future work:** in **code strings** use `\uXXXX` escapes so the literal
stays byte-identical (`≥`/`≤` appear in `sprintf()` labels on golden-tested plots — an ASCII
substitute would move rendered text and the golden data). In **roxygen prose** use ASCII,
since an escape renders literally in a comment.

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
| `readxl`, `org.Hs.eg.db`, `org.Mm.eg.db` | **CLAIMED by B5 — KEEP.** `read_metadata()` dispatches `.xlsx`/`.xls` to `readxl`; `annotate_genes()` maps ENSEMBL→SYMBOL/ENTREZID, and `species` is now a real argument so the human db is genuinely reachable. Finding 17 is **resolved: do not remove.** |
| `FactoMineR`, `factoextra` | **DO NOT ADD.** B5 deleted the FAMD module per `07` §8 rather than porting it, which makes `08` §4's `save_plot` → `famd_save_plot` rename moot. |
| `grDevices`, `tools` | added to Imports, `34000b7` |
| `patchwork` | added to Suggests, `4ac6539` |
| `tidyr` | removed from Imports, `4ac6539` |

---

## 6. Findings from the B3–B6 batch

19. **`data.table` was an undeclared hard dependency** of `io_helpers.R` and `provenance.R`
    (`fread`/`fwrite`) — absent from DESCRIPTION *and* from `07` §8's ledger. B5 rewrote both
    readers on `utils::read.delim` with header-based separator sniffing rather than add an
    unbudgeted dependency. Strictly better than the old extension-based dispatch, which
    mishandled semicolon CSVs and tab-delimited `.txt`.
20. **`08` §4 contradicts `07` §8 on FAMD** — the third doc-level contradiction this pass, same
    class as the one that made B1 delete a frozen export. `07` §8 says FactoMineR/factoextra are
    gone; `08` §4 orders a rename inside the module `07` deletes. B5 followed `07` per the
    "where they disagree" rule and said so.
21. **`set.seed()` before `solve_mwcsp()` is NOT sufficient** — §5a specified it and it was
    wrong. BioNet's BUM fit *inside* `gatom::scoreGraph()` draws random starts, so seeding only
    the solver left the **scores** irreproducible: one seed produced 40- and 53-node modules on
    consecutive runs. `gatom_module()` now seeds before `scoreGraph()` **and** before the
    solver. Only a realistic p-value distribution exposed this; the original all-significant
    toy input made scoring effectively deterministic and the seed-stability test passed
    vacuously. **A test that passes on degenerate input is not evidence.**
22. **BUM degeneracy is the GATOM smell test.** All-significant input → `Gene p-value threshold:
    1.000000` plus "parameters on the limit of the defined parameter space", because there is no
    uniform/null component to fit. Realistic input (18,766-symbol enzyme universe, mostly
    `rbeta(0.7, 1)`, planted kynurenine signal) → threshold `0.029485`, no warning, recovers
    IDO1/KYNU/KYAT3. Both regimes are kept as tests, the degenerate one explicitly labelled.
23. **B4 and B5 reached opposite conclusions on `patchwork` in the same batch** — B4 avoided it
    (undeclared) and returned a plain ggplot; B5 used it for `de_volcano_grid()`. Both are right
    locally: B4 composes one figure with a shared x axis, which a facet expresses natively; B5
    composes N independent caller-supplied ggplots with a collected legend, which has no
    ggplot-native equivalent. Parallel agents cannot see this class of divergence — **the
    integrator must look for it.**
24. **Golden proves the OLD code still works, not that the NEW code matches it.**
    `capture_golden.R` sources `scripts/`, not `R/`. B5 made this explicit and tested the new
    FDR-boundary path directly against a hand-computable fixture, including a **negative**
    assertion that the line is *not* at `-log10(p_cutoff)`. Step C's harness migration is what
    finally ties the two together — until then, a green golden run says less than it appears to.
25. **Three of four agents went idle without reporting** (B3 and B6 never did; B4 and B5 did
    after being asked twice). Same failure mode as B2. The gates are verifiable independently;
    the *reasoning* is not. B3's and B6's shim mappings in §0 above are integrator
    reconstructions, not agent statements — treat them accordingly.
26. **Two integrator misreads worth recording.** (a) A filtered `grep` dropped a
    `══ Warnings ══` header and made three gatom warnings look like newly-skipped tests; an
    unfiltered re-run showed one skip, correctly skipped. (b) B3's branch appeared to delete
    `STEP-C-REQUIREMENTS.md`, which was just divergence from a base commit made after it
    branched. **Check the file's history before concluding an agent went out of scope.**
