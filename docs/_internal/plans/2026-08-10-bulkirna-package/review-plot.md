# Renderer-layer review (R/gs-plot-*.R, R/gs-theme.R) — HEAD fc016ff

All findings were produced by building plots in the container (`pkgload::load_all("/pkg")`,
ggplot2 4.0.3) and inspecting `p$labels`, `levels()`, `attr(p, "gs_source")`,
`ggplot_build()` layer data / panel params, plus one rendered PNG. Old-vs-new quotes are
from `git show ff80de2^:scripts/...`.

### 1. `gs_plot_dot()` silently drops the *most significant* pathway when `padj == 0`  [severity: high]
**R/gs-plot-utils.R:173** (`df$neg_log_padj <- -log10(df$padj)`) feeding
**R/gs-plot-dot.R:114** (`aes(size = .data$neg_log_padj)`) — `padj == 0` gives `Inf`, the
size scale cannot map it, and the point is dropped with **no warning and no error**. The
significance outline layer (gs-plot-dot.R:120) drops it too.

**Evidence** — ran in container:
```
> pp <- gs_plot_dot(r0)   # padj = c(0, 1e-3, 0.2)
> ggplot_build(pp)$data[[1]][, c("x","y","size","fill")]
     x y size    fill
1  2.0 3  Inf #B35806      <- P1, padj = 0
2 -1.0 1   10 #98ABD2
3  0.5 2    2 #EECEBA
> ggsave("/out/dot0.png", pp)   # no warning, no error
```
The rendered PNG (`/tmp/dot0.png`) shows the "Set 1" row completely empty while its axis
tick is still drawn — a blank row that reads as "not tested".

**Failure scenario** — a `limma`/`gs_score` path or any p-value that underflows to 0 yields
`padj == 0`; that pathway is the strongest hit in the analysis and it is the one missing
from the published dotplot.

**Suggested fix** — clamp in `.gs_plot_frame()`:
`df$neg_log_padj <- -log10(pmax(df$padj, .Machine$double.xmin))` (or the smallest non-zero
`padj` in the frame), and document the clamp.

### 2. Two distinct pathways with the same formatted name collapse onto one axis row  [severity: high]
**R/gs-plot-utils.R:233-241** (`.gs_order_labels()`): the categorical axis is
`factor(df$label)`, and `label` is the *formatted, prefix-stripped, wrapped* pathway name.
`format_pathway_name(strip_prefix = TRUE)` maps `HALLMARK_APOPTOSIS` and `KEGG_APOPTOSIS`
to the identical string, so the two rows overplot and one pathway becomes invisible. Every
renderer that uses `.gs_order_labels` is affected (dot, bar, heatmap).

**Evidence** — ran in container:
```
> s3 <- attr(gs_plot_dot(r3, top_n = 4), "gs_source") # names: HALLMARK_APOPTOSIS, KEGG_APOPTOSIS, A_X, B_Y
  pathway_id     label
1         P1 Apoptosis
2         P2 Apoptosis
> levels(s3$label);  ggplot_build(p3)$layout$panel_params[[1]]$y$get_labels()
[1] "Apoptosis" "A X" "B Y"      # 4 pathways -> 3 rows
```
No test asserts label uniqueness (`grep label tests/testthat/test-gs-plot-dot.R` only checks
`labels$x` and column presence), and `ggplot_build()$data` carries no axis text, so the
golden gate cannot see it.

**Failure scenario** — a multi-database `gs_result` (`facet = "database"` is an advertised
use case, and `compare = "database"` more so) with the same set name in two collections:
one pathway's dot is drawn on top of the other's; the figure shows 19 of 20 pathways.

**Suggested fix** — key the factor on `pathway_id` and supply the display text via
`scale_y_discrete(labels = )`, or de-duplicate labels (append the database label) when
`anyDuplicated(df$label[!duplicated(df$pathway_id)])`.

### 3. `gs_plot_running()` legend shows raw MSigDB ids; every other renderer formats them  [severity: medium]
**R/gs-plot-running.R:316-340** (`.grs_labels()`) takes `x$pathway_name` and only *wraps*
it — `format_pathway_name()` is never called (`grep format_pathway_name R/gs-plot-running.R`
→ no hits). The dot/bar/heatmap renderers all route labels through
`format_pathway_name()` (gs-plot-utils.R:169).

**Evidence** — same object, two renderers, in container:
```
dot:     "p53 Pathway" || "TNF-alpha Signaling via NF-kappaB"
running: "HALLMARK_P53_PATHWAY" || "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
```
Invisible to tests because `test-gs-plot-running.R` fixtures already carry pretty names
("Alpha response").

**Failure scenario** — a figure panel pairing a dotplot with a running-sum plot of the same
pathways labels them two different ways; the running plot leaks the snake_case id the
package elsewhere promises never to show a reader.

**Suggested fix** — in `.grs_labels()`, pass `base` through
`format_pathway_name(..., strip_prefix = TRUE)` before `.grs_wrap()`, with a `strip_prefix`
argument for opt-out.

### 4. Heatmap x axis prints the literal placeholder `"contrast"` for a single-contrast result  [severity: medium]
**R/gs-plot-heatmap.R:85** — `by = "contrast"` is the default and `df$column <- df[[by]]` is
drawn verbatim. `gs_test()` defaults `contrast = "contrast"` (R/gs-test.R:41), so the common
single-contrast case gets an axis tick reading `contrast`. The function's own roxygen example
(gs-plot-heatmap.R:22, `gs_plot_heatmap(res)`) produces exactly this.

**Evidence** — ran the documented example in the container:
```
> res <- gs_test(ranks, db, min_size = 1, max_size = 10)
> res$contrast                                      # "contrast" "contrast"
> ggplot_build(gs_plot_heatmap(res))$layout$panel_params[[1]]$x$get_labels()
[1] "contrast"
```
(With a hand-built result whose `contrast` is `NA`, the tick reads `NA` — same class of
defect.)

**Failure scenario** — a user runs one contrast, calls `gs_plot_heatmap()`, and publishes a
heatmap whose single column is labelled "contrast".

**Suggested fix** — when the column values are all identical (or the placeholder / `NA`),
blank the x axis text for that dimension, or error asking for `by = "database"` — and change
the roxygen example to a multi-contrast object so `R CMD check` renders something sane.

### 5. Barplot fill limits changed from a fixed ±3.5 to per-figure symmetric limits, undocumented  [severity: medium]
**R/gs-plot-bar.R:75** — `if (is.null(limits)) limits <- .gs_symmetric_limits(df$stat)`.
The old renderer used a *fixed* scale:
```r
# scripts/GSEA/GSEA_plotting/gsea_barplot.R
nes_limits = c(-3.5, 3.5), ... oob = scales::squish
```
so the same NES had the same colour in every barplot of a project. Now each figure rescales
to its own max, and this divergence is not called out in the `@param limits` doc (which only
says "`NULL` uses symmetric limits"). `gs_plot_dot()` and both heatmap methods share the
behaviour.

**Evidence** — old source quoted above vs. new line 75; verified numerically: a result with
`max|stat| = 2` produces fills `#B35806 / #98ABD2` at ±2 (probe output, finding 1), i.e. the
saturated end of the ramp, where the old code would have shown a mid-tone.

**Failure scenario** — two barplots side by side in one figure (e.g. Hallmark vs GO from
`.gs_plot_all()`, which renders one per database): a pathway at NES 1.5 is pale in one panel
and saturated in the other, and a reader compares the colours.

**Suggested fix** — keep the data-driven default but say so in `@param limits`, and have
`.gs_plot_all()` pass one shared `limits` across the databases it renders.

### 6. `gs_source` row order does not match the order the bars/dots are drawn in  [severity: low]
**R/gs-plot-bar.R:92 / R/gs-plot-utils.R:284** — the attached table keeps *selection* order
(`.gs_select_top()` sort), while the drawn order comes from the factor levels set by
`.gs_order_labels()`. `gs_save()` writes this table next to the figure as its companion.

**Evidence** — ran in container (`sort_by = "stat"` default):
```
> sb <- attr(gs_plot_bar(r2, top_n = 4), "gs_source"); sb[, c("pathway_id","stat")]
  P1 -2.00 ; P4 2.00 ; P2 -0.67 ; P3 0.67      # table order
> levels(sb$label)   # figure order, bottom -> top
[1] "Set 1" "Set 2" "Set 3" "Set 4"            # P1, P2, P3, P4
```
**Failure scenario** — none for correctness; a reader diffing the CSV against the figure
row-by-row is misled. Maintainability/UX.

**Suggested fix** — reorder `df` by `levels(df$label)` (reversed for `coord_flip`) inside
`.gs_attach_source()`, or document that the table is in selection order.

### 7. Shared facet plumbing lives in `gs-plot-dot.R`, not `gs-plot-utils.R`  [severity: low]
**R/gs-plot-dot.R:157-205** — `.gs_check_compare()`, `.gs_facet_columns()` and
`.gs_add_facets()` are defined in the dotplot file but `gs_plot_bar()` calls
`.gs_facet_columns()`/`.gs_add_facets()` (gs-plot-bar.R:77, 91). Meanwhile
`gs_plot_heatmap.gs_result()` re-implements the same idea inline
(`df$column <- if (by == "database") df$database_label else df[[by]]`, heatmap:85) and
`gs_plot_heatmap.gs_matrix()` hand-writes its own `facet_grid()` (heatmap:187).

**Evidence** — grep/quote of the three call sites above; the file header of gs-plot-utils.R
states it holds "internal helpers shared by the gs_plot_* renderers", which these are.

**Failure scenario** — none; maintainability only (a facet-behaviour fix now has to be made
in two places, and the heatmap's third copy will drift).

**Suggested fix** — move the three helpers into `gs-plot-utils.R` and route the heatmap's
`column`/facet handling through `.gs_facet_columns()`.

### 8. Two small contract wrinkles: variance selection inside the matrix heatmap, and a duplicated 14 pt floor  [severity: low]
**R/gs-plot-heatmap.R:143** — `spread <- apply(as.matrix(m), 1L, stats::var, ...)`: the
renderer derives a statistic (row variance) to select rows, which is the one place in this
layer where "renderers never compute" bends. Related: the docs say `top_n` selects "by score
variance", but rows are then *ordered* by mean score (`.gs_order_labels(by = "score")`,
heatmap:177), so the visual order is not the selection order.
**R/gs-plot-running.R:507** — `base_size <- max(base_size, 14)` re-applies the floor that
`theme_bulki()` already enforces (gs-theme.R:26), i.e. the clamp the `.theme_bulki(floor =)`
split exists to keep in one place; it also silently overrides a `base_theme` the caller
passed at a smaller size.

**Evidence** — ran in container: `gs_plot_heatmap(gm, top_n = 3)` kept `S1,S3,S4` for row
variances `111.77, 0.05, 0.84, 2.77, 0.35` (correct top-3 set, but drawn in mean-score
order). Floor duplication is a direct quote of the two lines.

**Failure scenario** — none for the variance path (the selection is right); the duplicated
floor means `gs_plot_running(base_size = 11, base_theme = my_11pt_theme)` renders at 14 pt
chrome over an 11 pt base, i.e. mismatched text sizes.

**Suggested fix** — note in `@param top_n` that selection is by variance and display order by
mean score; drop the `max(base_size, 14)` in `.grs_theme()` and let `theme_bulki()` own the
floor (skip it entirely when `base_theme` is supplied).

## Categories with nothing to report
- **Colour/fill values**: the diverging ramp `#2166AC / #F7F7F7 / #B35806` and
  `oob = scales::squish` match the old `gsea_barplot()`/`gsea_dotplot()` defaults exactly;
  verified by quoting both. No finding beyond the limits issue (#5).
- **Factor level direction**: bar and dot both put the largest value at the top after
  `coord_flip` / discrete-y, matching the old `reorder(Description, NES)` behaviour —
  checked by printing `levels()` (finding 6's output). Nothing wrong.
- **`gs_plot_running()` colour keying**: the id-keyed `scale_colour_manual(values, breaks,
  labels)` really is permutation-proof; the reverse-palette test at
  test-gs-plot-running.R:74 exercises it and passes. Nothing found.
- **Theme element sizes**: `theme_bulki()` differs from the old
  `custom_minimal_theme_with_grid()` only by the intentional 14 pt floor, the optional
  `grid` argument, the `strip.*` settings, and `colour = "transparent"` in place of
  `color = NA` (correct for ggplot2 4.0). Nothing found.

## Not covered
I did not run the full `devtools::test(".")` suite or the golden gate — this review was
scoped to exercising the renderers directly. All eight findings above were reproduced in the
container or are direct old-vs-new source quotes; none are inferred from reading alone.
