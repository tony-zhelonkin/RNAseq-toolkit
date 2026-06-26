#' Unified GSEA Running-Sum Plot (modular, alignment-safe)
#'
#' Builds the canonical three-panel GSEA running-sum figure:
#'   1. running enrichment-score (ES) line,
#'   2. gene-hit tick rug,
#'   3. ranked-list-metric area,
#' stacked vertically and **kept left/right aligned to a shared x-axis by
#' construction**. Panel proportions are owned by a SINGLE composer knob
#' (`panel_heights`); per-panel `aspect.ratio` is NEVER used (it desyncs panel
#' widths — that is exactly the failure mode this refactor removes).
#'
#' ## Public surface (stable API)
#' `gsea_running_sum_plot()` is the only entry point. It returns a `patchwork`
#' that ALSO carries an `attr(p, "grs_restyle")` closure — a clean EXTENSION
#' interface letting a downstream project re-skin the figure (ES y-limits,
#' single collected legend outside-right, x-ticks on the bottom panel only,
#' hidden rug y-index labels, panel_heights, a project base theme) WITHOUT ever
#' indexing into patchwork panels (`p[[i]]`). See `?attr` usage below.
#'
#' ## Internals (prefix `.grs_`, NOT part of the stable API)
#' `.grs_select_ids`, `.grs_palette`, `.grs_apply_labels`, `.grs_build_panels`
#' (the ONE place that maps gseaplot2's subplot order -> named ES/rug/metric
#' roles), `.grs_panel_theme`, `.grs_hide_x`, `.grs_style_es`, `.grs_style_rug`,
#' `.grs_style_metric`, `.grs_compose` (assembles + aligns).
#'
#' @param gsea_obj A `gseaResult` object from clusterProfiler/fgsea.
#' @param gene_set_ids Integer vector (row indices) or character vector
#'                     (pathway IDs). Default picks top 5 by |NES|.
#' @param palette Optional color vector. If NULL, uses a 9-color vibrant palette.
#'                Kept UNNAMED on purpose (see note).
#' @param labels Optional named vector of legend labels (names = gene_set_ids).
#'               If NULL, uses each set's `Description`.
#' @param legend_pos Inside-panel legend position (used when
#'                   `legend_position = "inside"`). Default top-right.
#' @param base_size Base font size for the per-panel theme (ignored where
#'                  `base_theme` is supplied).
#' @param max_name_length Maximum character length for pathway names in legend.
#' @param title Optional title for the ES (top) panel.
#' @param panel_heights Length-3 numeric, the ES:rug:metric height ratio. The
#'                      SINGLE control over panel proportions. To get a near-
#'                      square ES panel, raise its share here and pick an output
#'                      device width/height accordingly (e.g. `save_figure(...,
#'                      width = 9, height = 7)`) — never via `aspect.ratio`.
#' @param legend_position One of `"inside"` (default; single legend pinned
#'                        inside the ES panel — today's behavior), `"right"`
#'                        (single legend collected OUTSIDE on the right, panel
#'                        widths stay aligned), or `"none"`.
#' @param es_ylim Optional length-2 numeric clamping the ES panel y-range via
#'                `coord_cartesian` (zoom, never drops data). NULL = free.
#' @param xticks `"bottom"` (default) keeps x ticks/labels only on the bottom
#'               (metric) panel; `"all"` shows them on every panel.
#' @param rug_ylabels If FALSE (default) hides the rug panel's meaningless
#'                    y-index labels; TRUE keeps them.
#' @param base_theme Optional complete ggplot2 theme used as the FOUNDATION for
#'                   every panel (e.g. a project `project_theme()`). The composer
#'                   applies it FIRST, then re-asserts the per-panel chrome
#'                   (tick/label suppression, inside-legend pin) ON TOP, so a
#'                   full base theme can never clobber the layout. NULL uses the
#'                   internal `.grs_panel_theme(base_size)`.
#'
#' @return A `patchwork` of three stacked, x-aligned panels, with an attached
#'   `attr(., "grs_restyle")` closure. Re-skin without panel indexing:
#'   \preformatted{
#'   p  <- gsea_running_sum_plot(g, gene_set_ids = ids)
#'   restyle <- attr(p, "grs_restyle")
#'   p2 <- restyle(es_ylim = c(-1, 1), legend_position = "right",
#'                 panel_heights = c(2.4, 0.7, 0.9), base_theme = my_theme())
#'   }
#' @export
#'
#' @note CRITICAL: the color palette vector is NOT named. Named colors break
#'       enrichplot::gseaplot2() for custom databases (SynGO, MitoCarta) whose
#'       pathway IDs differ from the Description field. Validated in the original
#'       syngo_running_sum_plot() implementation.
#' @note Alignment guarantee: panels are composed with `patchwork::wrap_plots`
#'       (and `plot_layout(guides = "collect")` for the outside-right legend),
#'       which aligns panel REGIONS column-wise regardless of differing y-axis
#'       label widths. No panel sets `aspect.ratio`, so widths cannot desync.
#' @note Refactored 2026-06-25: modular `.grs_` panel builders + composer +
#'       restyle closure (replaces the monolithic in-body stylise()).
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids = NULL,
                                  palette = NULL,
                                  labels = NULL,
                                  legend_pos = c(.98, .98),
                                  base_size = 14,
                                  max_name_length = 40,
                                  title = NULL,
                                  panel_heights = .GRS_PANEL_HEIGHTS,
                                  legend_position = c("inside", "right", "none"),
                                  es_ylim = NULL,
                                  xticks = c("bottom", "all"),
                                  rug_ylabels = FALSE,
                                  base_theme = NULL) {

  ## ------ 0. Validation ---------------------------------------------------
  stopifnot(methods::is(gsea_obj, "gseaResult"))
  legend_position <- match.arg(legend_position)
  xticks          <- match.arg(xticks)
  res_df <- gsea_obj@result
  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("`gsea_obj` has no result rows - nothing to plot.")
  }
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }
  if (!is.null(es_ylim)) {
    es_ylim <- as.numeric(es_ylim)
    stopifnot("`es_ylim` must be length-2 finite" = length(es_ylim) == 2 && all(is.finite(es_ylim)))
  }
  stopifnot("`panel_heights` must be length-3 positive" =
              length(panel_heights) == 3 && all(panel_heights > 0))

  ## ------ 1. Resolve IDs, palette, legend labels --------------------------
  gene_set_ids <- .grs_select_ids(res_df, gene_set_ids)
  n_sets       <- length(gene_set_ids)
  palette      <- .grs_palette(palette, n_sets)
  gsea_obj     <- .grs_apply_labels(gsea_obj, gene_set_ids, labels, max_name_length)

  ## ------ 2. Build the three RAW panels (the ONE index->role mapping) ------
  panels_raw <- tryCatch(
    .grs_build_panels(gsea_obj, gene_set_ids, palette, title),
    error = function(e) stop("Error generating running sum panels: ", conditionMessage(e)))

  ## ------ 3. Composer closure — rebuild from RAW panels with override knobs.
  ## This is the extension interface: style every panel, then COMPOSE + align.
  ## Captured defaults live in `.opts0`; overrides come via `...`. The closure
  ## re-attaches itself to its output, so the result is itself re-stylable
  ## (idempotent) — and a downstream project never indexes panels.
  .opts0 <- list(panel_heights   = panel_heights,
                 legend_position = legend_position,
                 es_ylim         = es_ylim,
                 xticks          = xticks,
                 rug_ylabels     = rug_ylabels,
                 base_theme      = base_theme,
                 base_size       = base_size,
                 legend_pos      = legend_pos)

  compose <- function(...) {
    o <- utils::modifyList(.opts0, list(...))
    show_x_top <- identical(o$xticks, "all")   # ES + rug show x only in "all" mode
    es  <- .grs_style_es(panels_raw$es, base_theme = o$base_theme, base_size = o$base_size,
                         es_ylim = o$es_ylim, show_x = show_x_top,
                         legend = o$legend_position, legend_pos = o$legend_pos)
    rug <- .grs_style_rug(panels_raw$rug, base_theme = o$base_theme, base_size = o$base_size,
                          show_x = show_x_top, rug_ylabels = o$rug_ylabels)
    met <- .grs_style_metric(panels_raw$metric, base_theme = o$base_theme,
                             base_size = o$base_size, show_x = TRUE)  # bottom always keeps x
    out <- .grs_compose(list(es = es, rug = rug, metric = met),
                        panel_heights = o$panel_heights, legend = o$legend_position)
    attr(out, "grs_restyle") <- compose
    out
  }

  compose()
}

## ===========================================================================
## Internals — prefix `.grs_`. Not exported / not part of the stable API.
## ===========================================================================

## Default ES:rug:metric height ratio. SINGLE source of truth for proportions.
.GRS_PANEL_HEIGHTS <- c(2.4, 0.7, 0.9)

#' Resolve `gene_set_ids` (indices or IDs) to a clean character vector of IDs.
.grs_select_ids <- function(res_df, gene_set_ids) {
  if (is.null(gene_set_ids)) {
    ord <- order(abs(res_df$NES), decreasing = TRUE)
    gene_set_ids <- ord[seq_len(min(5L, length(ord)))]
  }
  if (is.numeric(gene_set_ids)) {
    valid_idx <- gene_set_ids[gene_set_ids >= 1 & gene_set_ids <= nrow(res_df)]
    if (length(valid_idx) == 0)
      stop("All gene_set_ids are out of bounds (max: ", nrow(res_df), ")")
    if (length(valid_idx) < length(gene_set_ids))
      warning("Some gene_set_ids were out of bounds and removed")
    gene_set_ids <- res_df$ID[valid_idx]
  } else {
    missing <- gene_set_ids[!gene_set_ids %in% res_df$ID]
    if (length(missing) > 0) {
      warning("Some pathway IDs not found: ", paste(missing, collapse = ", "))
      gene_set_ids <- gene_set_ids[gene_set_ids %in% res_df$ID]
    }
  }
  gene_set_ids <- unique(gene_set_ids[!is.na(gene_set_ids)])
  if (length(gene_set_ids) == 0) stop("No valid gene-set IDs - nothing to plot.")
  gene_set_ids
}

#' Build an UNNAMED color vector of length `n_sets`.
#' CRITICAL: never name it (breaks gseaplot2 for custom DBs — see fn @note).
.grs_palette <- function(palette, n_sets) {
  if (is.null(palette)) {
    base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
    if (n_sets > length(base_pal)) grDevices::colorRampPalette(base_pal)(n_sets)
    else base_pal[seq_len(n_sets)]
  } else if (length(palette) < n_sets) {
    rep(palette, length.out = n_sets)
  } else {
    palette
  }
}

#' Set the plotted set Descriptions so the legend reads clean.
#' gseaplot2 draws the legend from `@result$Description`; we mutate a copy here so
#' `labels`/`max_name_length` take effect (rownames preserved for lookup). Long
#' names are WRAPPED onto multiple lines at `max_name_length`, never truncated with
#' an ellipsis — the full pathway name always stays readable (craft: don't truncate
#' labels). The legend grows taller, not clipped.
.grs_apply_labels <- function(gsea_obj, gene_set_ids, labels, max_name_length) {
  res <- gsea_obj@result
  idx <- match(gene_set_ids, res$ID)
  cur <- if (is.null(labels)) res$Description[idx] else unname(labels[gene_set_ids])
  res$Description[idx] <- vapply(cur, .grs_wrap_label,
                                 FUN.VALUE = character(1), width = max_name_length)
  gsea_obj@result <- res
  gsea_obj
}

#' Soft-wrap a label to <= `width` chars per line (break on whitespace; full text
#' preserved across `\n`-joined lines). Never inserts an ellipsis.
.grs_wrap_label <- function(x, width) {
  if (is.na(x)) return("(Unknown)")
  if (nchar(x) <= width) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}

#' Build the three RAW panels. The ONLY place that knows gseaplot2's subplot
#' order (1 = ES line, 2 = gene-hit rug, 3 = ranked metric) and binds it to
#' named roles, so nothing downstream needs positional `[[i]]` indexing.
.grs_build_panels <- function(gsea_obj, gene_set_ids, palette, title) {
  raw <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID    = gene_set_ids,
    title        = title %||% "",
    subplots     = c(1, 2, 3),
    pvalue_table = FALSE,
    rel_heights  = c(1.5, 0.5, 0.5),
    color        = palette)            # pass UNNAMED color vector
  list(es = raw[[1]], rug = raw[[2]], metric = raw[[3]])
}

#' Per-panel base theme (used when no project `base_theme` is supplied).
.grs_panel_theme <- function(base_size) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.title     = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(colour = "black", linewidth = 0.4),
      plot.margin      = ggplot2::margin(5, 10, 5, 5))
}

#' Suppress the x axis (ticks + text + title) on a panel.
.grs_hide_x <- function() {
  ggplot2::theme(axis.text.x  = ggplot2::element_blank(),
                 axis.title.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())
}

#' Apply the foundation theme (project `base_theme` or the internal default),
#' then return — caller layers per-panel chrome ON TOP so it cannot be clobbered.
.grs_apply_base <- function(p, base_theme, base_size) {
  p + (base_theme %||% .grs_panel_theme(base_size))
}

#' ES (top) panel — running enrichment score. Owns the SINGLE legend.
#' Legend keys are colored by gseaplot2's own `scale_color_manual(values = palette)`;
#' we do NOT add `override.aes` (it must match the key count exactly, and collides
#' when two truncated Descriptions map to the same legend entry — a custom-DB crash).
.grs_style_es <- function(p, base_theme, base_size, es_ylim, show_x,
                          legend, legend_pos) {
  p <- .grs_apply_base(p, base_theme, base_size) + ggplot2::labs(color = NULL)
  if (!is.null(es_ylim)) p <- p + ggplot2::coord_cartesian(ylim = es_ylim)
  p <- p + switch(legend,
    inside = ggplot2::theme(
      legend.position        = "inside",
      legend.position.inside = legend_pos,
      legend.justification   = c(1, 1),
      legend.background      = ggplot2::element_rect(fill = "white", colour = "grey90"),
      legend.key.size        = ggplot2::unit(0.8, "lines")),
    right  = ggplot2::theme(legend.position = "right"),  # collected outside by composer
    none   = ggplot2::theme(legend.position = "none"))
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' Rug (middle) panel — gene-hit ticks. Never carries its own legend (its colour
#' guide is dropped so `guides = "collect"` yields ONE legend from the ES panel).
.grs_style_rug <- function(p, base_theme, base_size, show_x, rug_ylabels) {
  p <- .grs_apply_base(p, base_theme, base_size) +
    ggplot2::labs(color = NULL) +
    ggplot2::guides(colour = "none", fill = "none") +
    ggplot2::theme(legend.position = "none")
  if (!isTRUE(rug_ylabels))
    p <- p + ggplot2::theme(axis.text.y  = ggplot2::element_blank(),
                            axis.title.y = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank())
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' Metric (bottom) panel — ranked-list statistic. Keeps the shared x-axis.
.grs_style_metric <- function(p, base_theme, base_size, show_x) {
  p <- .grs_apply_base(p, base_theme, base_size) +
    ggplot2::scale_y_continuous(n.breaks = 4) +
    ggplot2::theme(legend.position = "none")
  if (!isTRUE(show_x)) p <- p + .grs_hide_x()
  p
}

#' Compose + ALIGN. patchwork aligns panel regions column-wise regardless of
#' y-label widths; `guides = "collect"` keeps a right-side legend OUTSIDE the
#' stack so panel widths stay equal. No `aspect.ratio` anywhere by design.
.grs_compose <- function(panels, panel_heights, legend) {
  pw <- patchwork::wrap_plots(panels$es, panels$rug, panels$metric,
                              ncol = 1, heights = panel_heights)
  if (identical(legend, "right")) {
    pw <- pw + patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "right")
  }
  pw
}

## Portable null-coalesce (toolkit scripts are sourced, not packaged).
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}
