# DE-layer review — bulkiRNA @ fc016ff (read-only)

Everything below was either **executed** in `scdock-r-dev:v0.5.11` (ggplot2 4.0.3,
`pkgload::load_all("/pkg")`) or is a **side-by-side quote** of old vs new source. Each
finding says which.

---

### 1. `de_bfc_plot()` prints invented gene names ("Gene11", "Gene51") whenever `highlight_gene` overlaps the top-N set  [severity: high]

**R/de-bfc.R:86-89, 134** — highlighted genes are appended with `rbind()`, which
*uniquifies* duplicate row names (`Gene1` -> `Gene11`, `Gene5` -> `Gene51`). The very next
line, `lab_df[!duplicated(rownames(lab_df)), ]`, therefore never removes anything — the
names are no longer duplicates. `lab_df$.gene_label <- rownames(lab_df)` then puts the
mangled string on the figure, and `rownames(lab_df) %in% highlight_gene` is `FALSE` for it,
so the duplicate is also *not* bold/black.

**Evidence** — run in container:
```r
de <- data.frame(logFC=c(-3,-1,0.2,1.5,4,2.5,-2.2), B=c(4,-1,-3,0.5,6,3,2),
                 row.names=paste0("Gene",1:7))
p <- de_bfc_plot(de, fc_cutoff=1, top_n=3, highlight_gene=c("Gene5","Gene1"))
# repel layer data$.gene_label:
[1] "Gene5"  "Gene6"  "Gene4"  "Gene1"  "Gene7"  "Gene2"  "Gene11" "Gene51"
```
`de_volcano()` on the same input is correct (`"Gene6" "Gene4" "Gene7" "Gene2"` +
`"Gene1" "Gene5"` in the bold layer) because it *subtracts* highlights from `lab_df`
instead of appending; `de_md_plot()` is safe because it dedups on the `gene_id`
**column**, not on rownames. `de_bfc_plot` is the only one that dedups on rownames after
an `rbind`.

**Failure scenario** — `de_bfc_plot(topTable, highlight_gene = c("Ryr2","Atp2a2"))` where
those genes are also top-N: the published panel carries labels `Ryr21`, `Atp2a21` next to
the real ones. A reader looking up "Ryr21" finds nothing.

**Suggested fix** — adopt the volcano's split (`extra <- df[rownames(df) %in%
highlight_gene, ]`, remove those rows from `lab_df`, draw two layers), or dedup on an
explicit `gene_id` column added *before* the `rbind`. Add a test with
`highlight_gene` overlapping the top-N and assert `all(labels %in% rownames(de))`.

---

### 2. `de_volcano_grid()` silently drops the threshold caption from every figure  [severity: high]

**R/de-volcano.R:376** — `first_cap <- ggplot_build(plots[[1]])$layout$plot$labels$caption`.
Under ggplot2 4.0.3 the built layout no longer carries `$plot`, so `first_cap` is always
`NULL`. Line 388 then sets `caption = NULL` on **every** panel (the `keep_first_caption &&
i == 1` branch assigns `first_cap`, i.e. `NULL`, too), and the promotion at line 395 is
dead code. Both documented modes lose the caption.

**Evidence** — run in container:
```r
pa <- de_volcano(de, fc_cutoff=1, orientation="vertical")
pa$labels$caption
#> "Dashed lines: vert. – FDR ≤ 0.05 (p ≤ 0.01); horiz. – |log2FC| ≥ 1.0"
b <- ggplot_build(pa); !is.null(b$layout$plot)
#> FALSE                                 # ggplot2 4.0.3
g  <- de_volcano_grid(list(A=pa,B=pa));            str(g$patches$annotation$caption)  #> NULL
g2 <- de_volcano_grid(list(A=pa,B=pa), keep_first_caption=TRUE)
g2[[1]]$labels$caption                             #> "" (nothing)
str(g2$patches$annotation$caption)                 #> NULL
```
The `&`-vs-`+` question is *not* the problem — I checked both compose the annotation:
`wrap_plots(...) & plot_annotation(caption="HELLO")` -> `chr "HELLO"`.

**Failure scenario** — a 3-panel volcano row published with dashed horizontal lines and
no text anywhere saying they mean `FDR <= 0.05 (p <= 0.0012)`. The caption is the only
place the *realised raw-p* boundary is reported, which is the whole point of the
decision-by-FDR design (R/de-volcano.R:45-49).

**Suggested fix** — read the caption off the plot object, not the build:
`first_cap <- plots[[1]]$labels$caption`. Add a test asserting the caption survives in
both `keep_first_caption` modes — `tests/testthat/test-de-volcano.R:91-112` exercises the
grid but never looks at the caption.

---

### 3. `de_pca(xlim_abs=, ylim_abs=)` deletes samples instead of zooming  [severity: medium]

**R/de-pca.R:127-128** — `xlim()`/`ylim()` set **scale** limits, so any sample outside is
converted to `NA` and dropped; the only signal is a `Removed n rows ...` warning at
*print* time, after `ggplot_build()` has already succeeded.

**Evidence** — run in container:
```r
p <- de_pca(dge, colour_by="group", xlim_abs=0.001)
ggplot_build(p)$data[[1]][,c("x","y")]
#>    x          y      <- all four x are NA
print(p)
#> Warning: Removed 4 rows containing missing values ... (`geom_point()`)
#> Warning: Removed 4 rows containing missing values ... (`geom_text()`)
```
This matches the deleted `scripts/DE/plotPCA.R:129-130` (`xlim(-x_limit_val, ...)`), so it
is inherited, not introduced — but the new roxygen ("Optional symmetric axis limits")
does not warn about it, and the sibling renderers all use `coord_cartesian()`.

**Failure scenario** — a user pins `xlim_abs = 30` to make PCA panels comparable across
organs; the one outlier sample at PC1 = 42 disappears from the figure with no visible
error, and the plot looks like it has n-1 samples.

**Suggested fix** — clip with `coord_fixed(xlim = c(-x, x), ylim = c(-y, y))` (it already
calls `coord_fixed()`, so this is free), or keep `xlim()` and document/emit a warning
naming the dropped samples.

---

### 4. `x_breaks` retunes the p-value axis in `orientation = "vertical"`, and the fold-change axis it documents gets default breaks  [severity: medium]

**R/de-volcano.R:73, 211, 272-274** — `@param x_breaks Numeric spacing between
fold-change axis ticks.` In vertical mode fold change moves to *y* (no `scale_y_*` is
set, so ggplot2 picks breaks) and `x_breaks` is applied to `-log10(p)` instead. The
limit for the FC axis also switches from the breaks-aligned `fc_tick` (horizontal, line
264) to the unaligned `fc_max` (vertical, line 274), so the two orientations do not have
the same axis contract.

**Evidence** — run in container, same data, `x_breaks = 0.5`:
```r
# horizontal: 0.5 spacing on the FOLD-CHANGE axis (as documented)
[1] -4.0 -3.5 -3.0 ... 4.0
# vertical: 0.5 spacing on the -log10(p) axis; 17 ticks
[1] 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0
```

**Failure scenario** — `de_volcano(de, orientation = "vertical", x_breaks = 0.5)` on a
dataset reaching p = 1e-40 requests ~80 ticks on the p axis (unreadable black bar) while
leaving the fold-change axis on whatever ggplot2 chose.

**Suggested fix** — apply `x_breaks` to the fold-change axis in both orientations (rename
to `fc_breaks` and add a separate `p_breaks`), and use the breaks-aligned tick maximum
for the FC limit in both.

---

### 5. An all-`NA` `adj.P.Val` column produces a `NaN` threshold line and the caption `p ≤ NaN` instead of the documented "no genes pass" path  [severity: medium]

**R/de-volcano.R:150-159** — `sig_pvals <- de_results$P.Value[sig_logic]` keeps `NA`
positions, so `length(sig_pvals) > 0` is `TRUE` even when *nothing* is significant; the
guarded `max(..., na.rm = TRUE)` then returns `-Inf` and `sig_line` becomes `NaN` with
`draw_line = TRUE`. The `draw_line = FALSE` branch (and its "No genes pass FDR" caption
and annotation) is unreachable in this case.

**Evidence** — run in container:
```r
dn <- de; dn$adj.P.Val <- NA_real_
p <- de_volcano(dn, fc_cutoff = 1)
p$labels$caption
#> "Dashed lines: horiz. – FDR ≤ 0.05 (p ≤ NaN); vert. – |log2FC| ≥ 1.0"
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: NaNs produced
```
The old code had the identical shape (`scripts/DE/plot_standard_volcano.R:183-192`), so
this is inherited, not a regression.

**Failure scenario** — a `topTable` from a fit where `eBayes` yielded `NA` adjusted
p-values (all-zero-count genes, or a hand-assembled table): the figure caption reads
`p ≤ NaN` and no threshold line is drawn, with no error.

**Suggested fix** — filter first: `sig_pvals <- de_results$P.Value[which(sig_logic)]`
plus `sig_pvals <- sig_pvals[!is.na(sig_pvals)]`, then test `length(sig_pvals) > 0`.
Same for `p_max` (see #6).

---

### 6. `P.Value == 0` makes the volcano's y limit `Inf`  [severity: low]

**R/de-volcano.R:210-211, 265** — `p_max <- ceiling(max(-log10(df$P.Value)))` is `Inf`
when any p underflows to 0, and that `Inf` goes straight into
`coord_cartesian(ylim = c(0, p_max))`.

**Evidence** — run in container:
```r
d0 <- de; d0$P.Value[5] <- 0
ggplot_build(de_volcano(d0, fc_cutoff = 1))$layout$panel_params[[1]]$y.range
#> [1]   0 Inf
```

**Failure scenario** — a DE table whose p-values were written through a format that
underflows to 0 (or a `p.adjust` of a very large fit): the panel's y range is
`c(0, Inf)`, so every point is squashed at the bottom of a blank plot.

**Suggested fix** — compute `p_max`/`p_tick` on `-log10(pmax(P.Value, .Machine$double.xmin))`,
or drop non-finite values with a message.

---

### 7. The four DE renderers duplicate the same label-selection/highlight/shade block with three different, silently divergent dedup strategies  [severity: low]

**R/de-volcano.R:184-205, R/de-md.R:87-103, R/de-bfc.R:73-89** — `get_top()`,
`label_method` switch, highlight merge, `.de_shade()` recolour and the repel layer are
re-implemented three times. The three highlight merges are *not* equivalent:
volcano subtracts (correct), MD dedups on a `gene_id` column (correct), bfc dedups on
rownames after `rbind` (finding #1). `de-utils.R` currently holds only theme / shade /
require — the shared piece that actually caused a bug is not in it.

**Evidence** — old-vs-new quote. Old `create_standard_volcano()` used
`dplyr::bind_rows(lab_df, df[rownames(df) %in% highlight_gene, ]) |> dplyr::distinct()`
(`scripts/DE/plot_standard_volcano.R:241-245`) and old `create_MD_plot()` the same
(`scripts/DE/create_MD_plot.R:120-123`) — one strategy across the module. The refactor
split it into three, and the copy that changed `distinct()` -> `!duplicated(rownames())`
without changing `bind_rows` -> `rbind` broke. No behavioural claim beyond #1 here.

**Failure scenario** — none; maintainability only. (It is the mechanism behind #1.)

**Suggested fix** — hoist one `.de_label_set(df, label_method, top_n, rank_col,
highlight_gene)` into `R/de-utils.R` returning `list(regular, highlight)`, and have all
three renderers call it.

---

### 8. `fc_cutoff` defaults differ across the DE renderers (volcano 2, bfc 2, MD 1)  [severity: low]

**R/de-volcano.R:106, R/de-bfc.R:39, R/de-md.R:45** — the same biological threshold is
`2` (4-fold) in the volcano and B/FC plots and `1` (2-fold) in the MD plot, and nothing in
the roxygen flags the difference.

**Evidence** — source quote, three defaults:
`de_volcano(... fc_cutoff = 2 ...)`, `de_bfc_plot(... fc_cutoff = 2 ...)`,
`de_md_plot(... fc_cutoff = 1 ...)`. Matches the old scripts
(`plot_standard_volcano.R:124` = 2, `create_MD_plot.R:48` = 1), so it is inherited.

**Failure scenario** — a user producing a volcano and an MD plot for the same contrast
with default arguments gets guides at |log2FC| >= 2 and >= 1 in adjacent panels of one
figure, and the volcano's `annotate_counts` up/down numbers refer to a stricter cutoff
than the MD panel's dashed lines.

**Suggested fix** — pick one default (1 is the field convention) and state in each
`@param` that the three renderers now share it, or leave the values and add an explicit
"note: the MD default differs" line.

---

## Categories with nothing to report

* **Statistical correctness of the volcano's core decision rule** — nothing found. The
  raw-p axis / FDR decision split, the inclusive `<=`, the boundary line placed at
  `max(P.Value[adj.P.Val <= cutoff])`, `fixed_p_boundary`, and the `decision_by = "p"`
  branch all match `scripts/DE/plot_standard_volcano.R:171-202` line for line, and the
  caption reports the realised raw p. Verified by execution: with `fc_cutoff = 1` on a
  7-gene table, `annotate_counts` rendered `↑ 3   ↓ 2` and
  `table(p$data$cat)["p-value & Log2FC"]` was `5` — the counts do match the points drawn,
  and they are attached to the highest-priority *populated* category as documented.
* **`de_md_plot` significance/labelling** — nothing found beyond #8. FDR-only status,
  quadrant counts, `highlight_gene` bold/black selection (`ifelse` on `gene_id`, verified:
  `black, #803800, black, #00446B` / `bold, plain, bold, plain`) and the loess trend all
  reproduce `scripts/DE/create_MD_plot.R` faithfully.
* **`de_volcano_grid` shared-limits mechanism** — nothing found. The
  `plots[[i]]$coordinates <- shared_coord` assignment does work on ggplot2 4.0.3 S7
  objects; verified both sub-panels report `x.range -0.4 8.4 / y.range -4.4 4.4` from
  inputs whose own ranges were `-0.4 8.4` and `-0.15 3.15`. Only the caption path is
  broken (#2).
* **Argument forwarding (the `.de_theme()` bug class)** — nothing found. `.de_theme()`
  now returns `$text$size == 12` (verified) and no DE renderer exposes a `base_size`,
  `base_family` or `grid` argument that it then fails to pass on. I checked every
  parameter of all six public DE functions for a use site; the only parameter that reaches
  the wrong place is `x_breaks` (#4), and the only one that is accepted and ignored is
  `fixed_p_boundary` when `decision_by = "p"` (harmless, but undocumented).
* **Layer separation** — nothing found. `grep -rn "theme_bulki" R/` shows the DE files
  touch the gene-set layer only through the deliberate, documented `.theme_bulki`
  namespace lookup in `.de_theme()`; there is no `gs_*` call in any `de-*.R`, and no
  `gs-*.R` file calls a `de_*` function.
* **`R/dge.R`** — nothing found. `build_dge()` is a faithful, better-errored port of
  `scripts/General/dge_helpers.R` (same rounding tolerance `.Machine$double.eps^0.5`,
  same order checks). `annotate_genes()` is outside my named scope; I did not audit it.

## What I did not cover

* `annotate_genes()` (lives in `R/dge.R` but is an annotation concern, not DE) — read but
  not exercised; no network in the container.
* `de_pca_3d()` beyond construction — plotly output was not visually inspected.
* No visual/PDF regression comparison against `tests/golden/` outputs; my checks are
  numeric and structural.

---

## Addendum (requested): the no-significant-genes contract

### 9. The stale docs are wrong; current behaviour matches the old CODE exactly and is already documented in `de_volcano()`'s roxygen  [severity: low — docs only]

**Verdict: (a) deliberate, faithful to the old implementation, and documented.** Nothing to
fix in `R/`. `tests/README.md:109` and `tests/VOLCANO_FIX_SUMMARY.md` describe a behaviour
that neither old implementation ever had.

**file:line** — old: `scripts/DE/plot_standard_volcano.R:186-192, 302-313`; old vertical
twin: `scripts/DE/volcano_helpers.R:85-95, 194-205`; new: `R/de-volcano.R:150-159, 286-299`
and the doc sentence at `R/de-volcano.R:48-49`.

**Evidence** — old-vs-new source quote (no execution needed; the three code paths are
textually equivalent).

Old `create_standard_volcano()`, `git show ff80de2^:scripts/DE/plot_standard_volcano.R`:
```r
      sig_pvals <- de_results$P.Value[sig_logic]
      if (length(sig_pvals) > 0) { ... } else {
        # No significant genes - don't draw a line
        horiz_line <- NA
        draw_horiz_line <- FALSE
      }
...
  if (draw_horiz_line) {
    g <- g + ggplot2::geom_hline(yintercept = horiz_line, linetype = "dashed")
  } else if (decision_by == "fdr") {
    g <- g + ggplot2::annotate("text", ...,
      label = sprintf("No genes pass FDR ≤ %.2g", p_cutoff),
      size = 4, color = "darkred", fontface = "italic")
  }
```
Old `create_vertical_volcano()`, `git show ff80de2^:scripts/DE/volcano_helpers.R:85-95` —
byte-for-byte the same logic with `vert_line`/`draw_vert_line`, same italic annotation.

New `R/de-volcano.R:154-159, 286-299` — same: `draw_line <- FALSE`, no line, italic
`"No genes pass FDR ≤ %.2g"` in `darkred` at `size = 4`, plus a caption that says
`"No genes pass FDR ≤ 0.05. Dashed lines: vert. – |log2FC| ≥ ..."`.

There is **no** code path in either old file that puts a line at `-log10(p_cutoff)` in FDR
mode. `grep -n "log10(p_cutoff)"` over both old files hits only the `decision_by == "p"`
branch (`plot_standard_volcano.R:199`, `volcano_helpers.R:101`), where it is correct. The
old `VOLCANO_FIX_SUMMARY.md:57-61` claim —

> "### Case 2: No significant genes ✅ Correct behavior — Line drawn at `-log10(0.05)` as
> reference — **This looks like "misbehavior" but is actually correct!**"

— is prose that contradicts the file it was written to describe; likewise the snippet at
`VOLCANO_FIX_SUMMARY.md:147` (`cat("No significant genes - line at:", -log10(0.05))`) is a
*diagnostic print in a scratch script*, not an assertion about the plot. Treat both as
historical noise.

`tests/testthat/test-de-volcano.R:29-38` is therefore asserting the true, inherited
contract, and it asserts all three of its parts (no `geom_hline`, the italic annotation,
the caption).

**Failure scenario** — none; documentation only. The risk is the opposite direction: if
someone "fixes" the package to match the deleted README, FDR-mode volcanoes with zero hits
would grow a dashed line at `-log10(0.05)` that marks nothing — no colour boundary exists
there, which is precisely the misalignment the original fix removed.

**Suggested fix** — delete `VOLCANO_FIX_SUMMARY.md` as planned and state the contract in
`tests/README.md` as: *"FDR mode, zero genes passing: no threshold line; an italic dark-red
`No genes pass FDR ≤ x` annotation and a matching caption instead."* No change is needed in
`R/`: the roxygen already says it at `R/de-volcano.R:48-49` —

> "When no gene passes, no line is drawn and an italic annotation says so instead."

and the `@param caption` entry adds that `NULL` "generates one documenting the thresholds
actually drawn". The only wording I would add is one clause making explicit that this
applies to `decision_by = "fdr"` **only** — in `"p"` mode the line is unconditional at
`-log10(p_cutoff)` even when nothing is significant (`R/de-volcano.R:165-166`), which is
the asymmetry that probably spawned the confused README paragraph in the first place.
