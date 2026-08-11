#' Convert a bare `gseaResult` to a minimal `gs_result` for the plot shims
#'
#' The frozen plot names (`gsea_dotplot()`, `gsea_dotplot_facet()`,
#' `gsea_barplot()`, `gsea_running_sum_plot()`) only ever took a bare
#' `gseaResult` -- no `database`/`contrast` metadata -- so they cannot go
#' through `normalize_gsea_results()`, whose old formals require both with no
#' default. This builds just enough of a `gs_result` to drive the new
#' renderers: `pathway_id`/`pathway_name` from `ID`/`Description`, `stat`/
#' `stat_type` from `NES`/`"NES"`, `p_value`/`padj` from `pvalue`/`p.adjust`,
#' `n_genes`/`n_genes_tested` both from `setSize` (the old object carries no
#' better proxy for genes-tested), and `leading_edge` from splitting
#' `core_enrichment` on `"/"`. `database` and `contrast` are filled with a
#' constant placeholder since the frozen callers never supplied either.
#'
#' @param gsea_obj A `gseaResult` object.
#' @return A `gs_result`, with `attr(., "ranks")` set to `gsea_obj@geneList`
#'   and `attr(., "gene_sets")` set to `gsea_obj@geneSets`, for
#'   `gs_plot_running()`'s fallback.
#' @keywords internal
.dep_gsea_to_gs_result <- function(gsea_obj) {
  if (!methods::is(gsea_obj, "gseaResult")) {
    stop("This deprecated function requires a `gseaResult` object.",
         call. = FALSE)
  }
  res <- as.data.frame(gsea_obj@result)
  core <- res$core_enrichment
  leading_edge <- lapply(strsplit(core, "/"), function(g) g[nzchar(g)])

  df <- data.frame(
    pathway_id = as.character(res$ID),
    pathway_name = as.character(res$Description),
    n_genes = as.integer(res$setSize),
    n_genes_tested = as.integer(res$setSize),
    stat = as.numeric(res$NES),
    p_value = as.numeric(res$pvalue),
    padj = as.numeric(res$p.adjust),
    es = as.numeric(res$enrichmentScore),
    stringsAsFactors = FALSE
  )
  df$leading_edge <- leading_edge

  out <- gs_result(
    df,
    database = "gsea",
    contrast = "gsea",
    method = "fgsea",
    stat_type = "NES"
  )
  attr(out, "ranks") <- gsea_obj@geneList
  attr(out, "gene_sets") <- gsea_obj@geneSets
  out
}

#' Enhanced GSEA dotplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_dot()] instead. This shim reproduces the old
#' formals of `gsea_dotplot()` verbatim and forwards to `gs_plot_dot()`.
#'
#' @param gsea_obj GSEA result object.
#' @param filterBy Method to sort/filter results: `"p.adjust"` (default),
#'   `"NES"`, `"NES_positive"`, `"NES_negative"`.
#' @param sortBy Secondary display sort (`"GeneRatio"` or `"p.adjust"`).
#' @param showCategory Number of top pathways to show.
#' @param padj_cutoff Significance threshold for highlighting.
#' @param title Plot title.
#' @param wrap_width Width for text wrapping.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits Numeric vector of length 2 for symmetric NES limits.
#' @param min.dotSize Minimum dot size.
#' @param max.dotSize Maximum dot size.
#' @param highlight_sig Whether to highlight significant points.
#' @param highlight_threshold FDR threshold for highlighting; `NULL` uses
#'   `padj_cutoff`.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#' @param use_gradient Logical, use continuous gradient for NES.
#' @param ... Unused; absorbs no old formals but kept for forward safety.
#'
#' @return A ggplot object, as returned by [gs_plot_dot()].
#' @note `sortBy` has no counterpart in `gs_plot_dot()` -- the new renderer
#'   always orders the y axis by `stat` magnitude within the selection -- so
#'   it is accepted and ignored.
#' @note `use_gradient = FALSE` (the old "binary colour" mode) has no
#'   counterpart: `gs_plot_dot()` always uses the continuous diverging
#'   gradient. Accepted and ignored; the plot always renders with a gradient.
#' @export
gsea_dotplot <- function(
    gsea_obj,
    filterBy = "p.adjust",
    sortBy = "GeneRatio",
    showCategory = 10,
    padj_cutoff = 0.05,
    title = "GSEA Dotplot",
    wrap_width = 50,
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    highlight_threshold = NULL,
    strip_prefix = TRUE,
    use_gradient = TRUE,
    nes_limits = NULL) {
  .Deprecated("gs_plot_dot")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  sort_by <- "padj"
  direction <- "both"
  if (identical(filterBy, "NES")) {
    sort_by <- "stat"
  } else if (identical(filterBy, "NES_positive")) {
    sort_by <- "stat"
    direction <- "up"
  } else if (identical(filterBy, "NES_negative")) {
    sort_by <- "stat"
    direction <- "down"
  }

  highlight <- if (isTRUE(highlight_sig)) {
    highlight_threshold %||% padj_cutoff
  } else {
    NULL
  }

  gs_plot_dot(
    x,
    top = showCategory,
    sort_by = sort_by,
    direction = direction,
    highlight = highlight,
    size_range = c(min.dotSize, max.dotSize),
    limits = nes_limits,
    colours = c(neg_color, mid_color, pos_color),
    wrap_width = wrap_width,
    strip_prefix = strip_prefix,
    title = title
  )
}

#' Enhanced GSEA faceted dotplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_dot()] with `facet = "direction"` instead. This
#' shim reproduces the old formals of `gsea_dotplot_facet()` verbatim.
#'
#' @param gsea_obj GSEA result object.
#' @param showCategory Number of pathways to show per direction.
#' @param padj_cutoff Adjusted p-value cutoff used for significance
#'   highlighting.
#' @param title Plot title.
#' @param wrap_width Width for text wrapping.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits NES colour scale limits.
#' @param min.dotSize Minimum dot size.
#' @param max.dotSize Maximum dot size.
#' @param highlight_sig Whether to highlight significant points.
#' @param highlight_threshold FDR threshold for highlighting; `NULL` uses
#'   `padj_cutoff`.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#'
#' @return A ggplot object, as returned by [gs_plot_dot()].
#' @export
gsea_dotplot_facet <- function(
    gsea_obj,
    showCategory = 10,
    padj_cutoff = 0.05,
    title = "GSEA Faceted Dotplot",
    wrap_width = 50,
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    nes_limits = c(-3.5, 3.5),
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    highlight_threshold = NULL,
    strip_prefix = TRUE) {
  .Deprecated("gs_plot_dot")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  highlight <- if (isTRUE(highlight_sig)) {
    highlight_threshold %||% padj_cutoff
  } else {
    NULL
  }

  gs_plot_dot(
    x,
    top = showCategory,
    sort_by = "padj",
    direction = "both",
    facet = "direction",
    highlight = highlight,
    size_range = c(min.dotSize, max.dotSize),
    limits = nes_limits,
    colours = c(neg_color, mid_color, pos_color),
    wrap_width = wrap_width,
    strip_prefix = strip_prefix,
    title = title
  )
}

#' Enhanced GSEA barplot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_bar()] instead. This shim reproduces the old
#' formals of `gsea_barplot()` verbatim.
#'
#' @param gsea_obj GSEA result object.
#' @param padj_cutoff Adjusted p-value cutoff; a hard filter, as in the old
#'   function.
#' @param top_n Number of pathways to show.
#' @param title Plot title.
#' @param neg_color Colour for negative NES.
#' @param mid_color Colour for zero NES.
#' @param pos_color Colour for positive NES.
#' @param nes_limits NES colour scale limits.
#' @param strip_prefix Logical, strip common prefixes like `"HALLMARK_"`.
#'
#' @return A ggplot object, as returned by [gs_plot_bar()].
#' @note The old function pre-filtered on `padj_cutoff` and never outlined
#'   individual bars, so `highlight` is forced to `NULL` here -- outlining an
#'   already-significance-filtered set would be a no-op difference from the
#'   old figure, but is called out explicitly rather than left implicit.
#' @export
gsea_barplot <- function(
    gsea_obj,
    padj_cutoff = 0.05,
    top_n = 30,
    title = "GSEA NES Barplot",
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    nes_limits = c(-3.5, 3.5),
    strip_prefix = TRUE) {
  .Deprecated("gs_plot_bar")

  x <- .dep_gsea_to_gs_result(gsea_obj)

  gs_plot_bar(
    x,
    top = top_n,
    sort_by = "stat",
    direction = "both",
    padj_max = padj_cutoff,
    highlight = NULL,
    limits = nes_limits,
    colours = c(neg_color, mid_color, pos_color),
    strip_prefix = strip_prefix,
    title = title
  )
}

#' Unified GSEA running-sum plot (deprecated)
#'
#' @description
#' Deprecated: use [gs_plot_running()] instead. This shim reproduces the old
#' formals of `gsea_running_sum_plot()` verbatim.
#'
#' @param gsea_obj A `gseaResult` object from clusterProfiler/fgsea.
#' @param gene_set_ids Integer vector (row indices) or character vector
#'   (pathway IDs). `NULL` (default) picks the new renderer's default top 5
#'   by \eqn{|NES|}.
#' @param palette Optional colour vector, keyed as documented in
#'   [gs_plot_running()].
#' @param labels Optional legend labels, keyed as documented in
#'   [gs_plot_running()].
#' @param legend_pos Inside-panel legend position (used when
#'   `legend_position = "inside"`).
#' @param base_size Base font size for the per-panel theme.
#' @param max_name_length Maximum character length for pathway names in the
#'   legend.
#' @param title Optional title for the ES (top) panel.
#' @param panel_heights Length-3 numeric, the ES:rug:metric height ratio.
#' @param legend_position One of `"inside"` (default), `"right"` or `"none"`.
#' @param es_ylim Optional length-2 numeric; absorbed (see `@note`).
#' @param xticks `"bottom"` (default) or `"all"`; absorbed (see `@note`).
#' @param rug_ylabels Logical; absorbed (see `@note`).
#' @param base_theme Optional complete ggplot2 theme used as the foundation
#'   for every panel.
#'
#' @return A ggplot object, as returned by [gs_plot_running()].
#' @note `es_ylim` is absorbed and ignored: the new renderer would need a
#'   per-facet zoom to honour it, which was deliberately dropped when
#'   `gs_plot_running()` was written.
#' @note `xticks` and `rug_ylabels` are absorbed and ignored: the new
#'   renderer's layout is now structurally equivalent to the old
#'   `xticks = "bottom"`, `rug_ylabels = FALSE` combination (x ticks only on
#'   the bottom panel; the rug panel's y-index labels are always blanked), so
#'   both old values of each collapse onto that one good layout.
#' @note `gsea_obj@geneList` and `gsea_obj@geneSets` are attached to the
#'   converted `gs_result` as the `ranks` and `gene_sets` attributes, which
#'   [gs_plot_running()] falls back to when its own `ranks`/`db` arguments are
#'   `NULL` -- so this shim never duplicates that lookup logic itself.
#' @export
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids = NULL,
                                  palette = NULL,
                                  labels = NULL,
                                  legend_pos = c(.98, .98),
                                  base_size = 14,
                                  max_name_length = 40,
                                  title = NULL,
                                  panel_heights = c(2.4, 0.7, 0.9),
                                  legend_position = c("inside", "right", "none"),
                                  es_ylim = NULL,
                                  xticks = c("bottom", "all"),
                                  rug_ylabels = FALSE,
                                  base_theme = NULL) {
  .Deprecated("gs_plot_running")
  legend_position <- match.arg(legend_position)
  xticks <- match.arg(xticks)

  x <- .dep_gsea_to_gs_result(gsea_obj)

  gs_plot_running(
    x,
    pathways = gene_set_ids,
    labels = labels,
    palette = palette,
    panel_heights = panel_heights,
    title = title,
    base_size = base_size,
    max_name_length = max_name_length,
    legend_position = legend_position,
    legend_pos = legend_pos,
    base_theme = base_theme
  )
}

#' Custom minimal theme with grid (deprecated)
#'
#' @description
#' Deprecated: use [theme_bulki()] instead. This shim reproduces the old
#' formals of `custom_minimal_theme_with_grid()` verbatim and forwards to
#' the grid variant of the `theme_bulki()` theme.
#'
#' @param base_size Base font size for the theme.
#' @param base_family Base font family for the theme.
#'
#' @return A ggplot2 theme object.
#' @note The old theme had no base-size floor and rendered at its documented
#'   12 pt default; `theme_bulki()` raises anything below 14 pt to 14 pt on
#'   purpose. To keep this deprecated path byte-faithful, the shim routes
#'   through the internal `.theme_bulki(floor = NULL)` rather than the public
#'   `theme_bulki()`, so `base_size = 12` really does render at 12 pt. The
#'   public theme keeps its floor. No gate could have caught the difference:
#'   the golden record for this case stores only the 144 theme element *names*,
#'   never their sizes.
#' @export
custom_minimal_theme_with_grid <- function(base_size = 12, base_family = "") {
  .Deprecated("theme_bulki")
  .theme_bulki(base_size = base_size, base_family = base_family,
               grid = TRUE, floor = NULL)
}

#' Plot all GSEA results for a given analysis (deprecated)
#'
#' @description
#' Deprecated: the new golden path is to call the [gs_plot_dot()] /
#' [gs_plot_bar()] renderer you want and save it with [gs_save()]. This shim
#' reproduces the old formals of `plot_all_gsea_results()` verbatim and
#' forwards to the internal `.gs_plot_all()`, which holds the ported body.
#'
#' @param gsea_list List of GSEA results.
#' @param analysis_name Name of the analysis.
#' @param out_root Output root directory.
#' @param n_pathways Number of pathways to show.
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param expr_data Expression data matrix for heatmaps. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param sample_annotation Sample annotation data.frame. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param sample_order Order of samples for plots. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#' @param ann_colors Colours for annotations. Not supported by
#'   `.gs_plot_all()`; accepted and ignored (see `@note`).
#'
#' @return A character vector of every written path, invisibly, as returned
#'   by `.gs_plot_all()`.
#' @note `gsea_list` (a named list of `gseaResult`/data-frame objects keyed by
#'   database) is converted to a single `gs_result` via the internal
#'   `.dep_gsea_to_gs_result()`, with each element's list name substituted in
#'   as its `database` column so `.gs_plot_all()`'s per-database loop still
#'   works.
#' @note `expr_data`, `sample_annotation`, `sample_order` and `ann_colors`
#'   drove heatmap panels in the old function; `.gs_plot_all()` (the ported
#'   body behind this shim) does not render a heatmap, so these four formals
#'   are accepted and ignored. No heatmap is written where the old function
#'   would have written one.
#' @export
plot_all_gsea_results <- function(
    gsea_list,
    analysis_name,
    out_root,
    n_pathways = 20,
    padj_cutoff = 0.05,
    expr_data = NULL,
    sample_annotation = NULL,
    sample_order = NULL,
    ann_colors = NULL) {
  .Deprecated(".gs_plot_all")

  if (is.null(gsea_list) || length(gsea_list) == 0) {
    return(invisible(character(0)))
  }

  parts <- lapply(names(gsea_list), function(nm) {
    part <- .dep_gsea_to_gs_result(gsea_list[[nm]])
    part[["database"]] <- nm
    part
  })
  x <- do.call(rbind, parts)

  .gs_plot_all(
    x,
    out_dir = file.path(out_root, analysis_name),
    name = analysis_name,
    top = n_pathways,
    padj_cutoff = padj_cutoff
  )
}

#' Save GSEA results to a text log file (deprecated)
#'
#' @description
#' Deprecated: use [gs_save()] (source-table export) or the internal
#' `.gs_write_log()` directly. This shim reproduces the old formals of
#' `save_gsea_log()` verbatim and forwards to the internal `.gs_write_log()`,
#' which holds the ported body.
#'
#' @param gsea_obj GSEA result object.
#' @param filename Output filename (with or without path).
#' @param padj_cutoff Adjusted p-value cutoff.
#' @param dir Output directory (optional).
#'
#' @return `path`, invisibly, as returned by `.gs_write_log()`.
#' @export
save_gsea_log <- function(
    gsea_obj,
    filename,
    padj_cutoff = 0.05,
    dir = NULL) {
  .Deprecated(".gs_write_log")

  if (!is.null(dir)) {
    filename <- file.path(dir, filename)
  }
  x <- .dep_gsea_to_gs_result(gsea_obj)
  .gs_write_log(x, filename, padj_cutoff = padj_cutoff)
}

#' Create a standard volcano plot (deprecated)
#'
#' @description
#' Deprecated: use [de_volcano()] instead. This shim reproduces the old
#' formals of `create_standard_volcano()` verbatim and forwards to
#' `de_volcano()`, which has the identical formals plus one new argument,
#' `orientation`, placed after `...` so no positional call to this shim can
#' shift onto it.
#'
#' @inheritParams de_volcano
#' @param de_results Data frame whose rownames are gene IDs and that contains
#'   at least `logFC`, `P.Value`, `adj.P.Val`.
#'
#' @return A `ggplot2` object.
#' @export
create_standard_volcano <- function(
    de_results,
    decision_by = c("fdr", "p"),
    p_cutoff = 0.05,
    fc_cutoff = 2,
    top_n = 5,
    highlight_gene = NULL,
    label_method = "top",
    x_breaks = 1,
    title = "Volcano plot",
    subtitle = NULL,
    caption = NULL,
    fixed_p_boundary = NULL,
    color_palette = c(
      "NS"               = "#7F7F7F",
      "Log2FC"           = "#0173B2",
      "p-value"          = "#029E73",
      "p-value & Log2FC" = "#D55E00"
    ),
    show_grid = FALSE,
    max.overlaps = 10,
    annotate_counts = FALSE,
    ...) {
  .Deprecated("de_volcano")
  de_volcano(
    de_results = de_results,
    decision_by = decision_by,
    p_cutoff = p_cutoff,
    fc_cutoff = fc_cutoff,
    top_n = top_n,
    highlight_gene = highlight_gene,
    label_method = label_method,
    x_breaks = x_breaks,
    title = title,
    subtitle = subtitle,
    caption = caption,
    fixed_p_boundary = fixed_p_boundary,
    color_palette = color_palette,
    show_grid = show_grid,
    max.overlaps = max.overlaps,
    annotate_counts = annotate_counts,
    ...
  )
}

#' Create a mean-difference (MD) plot (deprecated)
#'
#' @description
#' Deprecated: use [de_md_plot()] instead. This shim reproduces the old
#' formals of `create_MD_plot()` verbatim and forwards to `de_md_plot()`,
#' whose formals are identical (the old `fdr_cutoff` name is unchanged).
#'
#' @inheritParams de_md_plot
#'
#' @return A `ggplot2` object.
#' @export
create_MD_plot <- function(
    fit,
    coef,
    de_results = NULL,
    fc_cutoff = 1,
    fdr_cutoff = 0.05,
    top_n = 5,
    highlight_gene = NULL,
    label_method = "top",
    max.overlaps = 10,
    title = NULL,
    color_palette = c(
      Up   = "#D55E00",
      Down = "#0072B2",
      NS   = "#999999"
    ),
    show_grid = FALSE,
    show_quadrant_counts = TRUE) {
  .Deprecated("de_md_plot")
  de_md_plot(
    fit = fit,
    coef = coef,
    de_results = de_results,
    fc_cutoff = fc_cutoff,
    fdr_cutoff = fdr_cutoff,
    top_n = top_n,
    highlight_gene = highlight_gene,
    label_method = label_method,
    max.overlaps = max.overlaps,
    title = title,
    color_palette = color_palette,
    show_grid = show_grid,
    show_quadrant_counts = show_quadrant_counts
  )
}
