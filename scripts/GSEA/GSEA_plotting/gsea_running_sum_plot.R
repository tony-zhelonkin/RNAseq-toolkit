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
#' hidden rug y-index labels, panel_heights, colours, labels, a project base
#' theme) WITHOUT ever indexing into patchwork panels (`p[[i]]`). See `?attr`
#' usage below.
#'
#' ## Internals (prefix `.grs_`, NOT part of the stable API)
#' `.grs_select_ids`, `.grs_resolve_labels`, `.grs_palette` (re-keys any accepted
#' palette keying onto the final plotted label), `.grs_apply_labels`,
#' `.grs_raw_panels` / `.grs_build_panels` (the ONE place that maps gseaplot2's
#' subplot order -> named ES/rug/metric roles), `.grs_set_color_scale` (the ONE
#' place that owns the colour scale), `.grs_panel_theme`, `.grs_hide_x`,
#' `.grs_style_es`, `.grs_style_rug`, `.grs_style_metric`, `.grs_compose`
#' (assembles + aligns).
#'
#' @param gsea_obj A `gseaResult` object from clusterProfiler/fgsea.
#' @param gene_set_ids Integer vector (row indices) or character vector
#'                     (pathway IDs). Default picks top 5 by |NES|.
#' @param palette Optional color vector. If NULL, uses a 9-color vibrant palette.
#'                Three keyings are accepted and resolved in this precedence
#'                order: (1) named by the FINAL plotted label (the value the
#'                colour aesthetic carries — i.e. `labels` after wrapping, else
#'                `Description`), (2) named by pathway ID, (3) unnamed, zipped to
#'                `gene_set_ids` **in the order the caller supplied** (recycled or
#'                truncated to length). Whichever you pass, the palette is
#'                re-keyed internally onto the final label, so colours can never
#'                be permuted by ggplot's alphabetical level order (see note).
#'                Names matching neither labels nor IDs warn and fall back to (3).
#' @param labels Optional legend labels for the plotted sets. Named by pathway ID
#'               (keys not present keep that set's `Description`, so a partial map
#'               degrades gracefully), or unnamed and zipped to `gene_set_ids` in
#'               the caller's declared order. NULL uses each set's `Description`.
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
#'                 panel_heights = c(2.4, 0.7, 0.9), base_theme = my_theme(),
#'                 palette = c("#A6611A", "#2166AC"))
#'   }
#'   Restyle-overridable knobs: `panel_heights`, `legend_position`, `es_ylim`,
#'   `xticks`, `rug_ylabels`, `base_theme`, `base_size`, `legend_pos`, `palette`
#'   (re-resolved and re-applied at STYLE time, so no panel rebuild is needed),
#'   plus the three build-time knobs `labels`, `max_name_length` and `title`,
#'   which DO force a raw-panel rebuild — they mutate `Description` inside the
#'   gseaResult before `gseaplot2()` ever sees it — and are handled by re-running
#'   the builder rather than being silently ignored. Overrides are per-call, not
#'   cumulative: every `restyle()` starts from the arguments the original call was
#'   built with.
#' @export
#'
#' @note CRITICAL colour-keying invariant: every per-set argument (`palette`,
#'       `labels`) is keyed by **pathway ID or final plotted label**, and this
#'       function re-keys it onto whatever the plot actually maps on — which for
#'       `enrichplot::gseaplot2()` is `Description`. NEVER rely on ggplot's
#'       discrete level order: ggplot sorts levels ALPHABETICALLY, so an unnamed
#'       palette forwarded to `gseaplot2(color = )` (which does
#'       `scale_color_manual(values = color)`, i.e. positionally) permutes the
#'       colours whenever alphabetical order differs from the caller's declared
#'       order — e.g. `c("… up", "… down")` renders SWAPPED because "… down" sorts
#'       first, and the name-keyed `labels` keep the legend TEXT correct, so the
#'       swap reads as a deliberate colour choice rather than a bug. This function
#'       therefore OWNS the colour scale (`.grs_set_color_scale()`: `values` +
#'       `breaks = names()`, matching toolkit house style in
#'       `plot_standard_volcano.R`, `volcano_helpers.R`, `create_MD_plot.R`,
#'       `create_fc_b_plot.R`) instead of delegating to gseaplot2's positional
#'       `color=`.
#' @note WHY the old "never name the palette" rule existed (do not reintroduce
#'       it): the Feb-2026 breakage was a palette named by PATHWAY ID, which never
#'       matches the mapped aesthetic for custom databases (SynGO, MitoCarta)
#'       whose IDs differ from `Description` — giving grey50 keys and a
#'       "No shared levels found between `names(values)` …" warning. That
#'       observation was correct, but banning names was the wrong remedy: it
#'       forbade the only keying that works and traded a loud failure for a silent
#'       one. ID-named palettes are now re-keyed, not forbidden.
#' @note Alignment guarantee: panels are composed with `patchwork::wrap_plots`
#'       (and `plot_layout(guides = "collect")` for the outside-right legend),
#'       which aligns panel REGIONS column-wise regardless of differing y-axis
#'       label widths. No panel sets `aspect.ratio`, so widths cannot desync.
#' @note Refactored 2026-06-25: modular `.grs_` panel builders + composer +
#'       restyle closure (replaces the monolithic in-body stylise()).
#' @note Fixed 2026-07-26: palette resolution re-keyed onto the final plotted
#'       label and the colour scale moved into the styling layer, so the caller's
#'       declared order — not ggplot's alphabetical order — decides both colour
#'       assignment and legend key order.
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

  ## ------ 1. Resolve IDs, then build the three RAW panels ------------------
  ## Label resolution comes FIRST and is shared: the final (wrapped) label is
  ## both what gseaplot2 draws and the key the colour scale is built on, so it is
  ## computed once and handed to the palette resolver and the panel builder.
  gene_set_ids <- .grs_select_ids(res_df, gene_set_ids)
  gsea_raw     <- gsea_obj                    # pristine; label overrides re-mutate it
  built        <- .grs_raw_panels(gsea_raw, res_df, gene_set_ids,
                                  labels, max_name_length, title)

  ## ------ 2. Composer closure — restyle from RAW panels with override knobs.
  ## This is the extension interface: style every panel, then COMPOSE + align.
  ## Captured defaults live in `.opts0`; overrides come via `...`. The closure
  ## re-attaches itself to its output, so the result is itself re-stylable
  ## (idempotent) — and a downstream project never indexes panels.
  ##
  ## `palette` is captured UNRESOLVED and re-keyed inside compose(), because its
  ## resolution depends on the labels; the resolved (label-named) vector is then
  ## applied at STYLE time, so overriding colours needs no panel rebuild.
  ## `labels`/`max_name_length`/`title` are baked into the gseaResult before
  ## gseaplot2() runs, so overriding those re-runs the builder (one gseaplot2
  ## call) instead of being silently ignored.
  .opts0 <- list(panel_heights   = panel_heights,
                 legend_position = legend_position,
                 es_ylim         = es_ylim,
                 xticks          = xticks,
                 rug_ylabels     = rug_ylabels,
                 base_theme      = base_theme,
                 base_size       = base_size,
                 legend_pos      = legend_pos,
                 palette         = palette,
                 labels          = labels,
                 max_name_length = max_name_length,
                 title           = title)

  compose <- function(...) {
    o <- utils::modifyList(.opts0, list(...))
    b <- if (identical(o$labels,          .opts0$labels) &&
             identical(o$max_name_length, .opts0$max_name_length) &&
             identical(o$title,           .opts0$title)) {
      built                                   # nothing build-time changed
    } else {
      .grs_raw_panels(gsea_raw, res_df, gene_set_ids,
                      o$labels, o$max_name_length, o$title)
    }
    pal        <- .grs_palette(o$palette, gene_set_ids, b$set_labels)
    show_x_top <- identical(o$xticks, "all")  # ES + rug show x only in "all" mode
    es  <- .grs_style_es(b$panels$es, base_theme = o$base_theme, base_size = o$base_size,
                         es_ylim = o$es_ylim, show_x = show_x_top,
                         legend = o$legend_position, legend_pos = o$legend_pos,
                         palette = pal)
    rug <- .grs_style_rug(b$panels$rug, base_theme = o$base_theme, base_size = o$base_size,
                          show_x = show_x_top, rug_ylabels = o$rug_ylabels,
                          palette = pal)
    met <- .grs_style_metric(b$panels$metric, base_theme = o$base_theme,
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

#' Resolve the FINAL plotted label for every selected set.
#' Returns a character vector NAMED BY PATHWAY ID whose values are the exact
#' strings that end up in `@result$Description` — i.e. what the colour aesthetic
#' will carry, and therefore the only correct key for a manual colour scale.
#' `labels` may be named by pathway ID (unmatched keys keep that set's
#' `Description`, so a partial map degrades gracefully instead of silently
#' becoming "(Unknown)") or unnamed, in which case it is zipped to `gene_set_ids`
#' in the caller's declared order.
.grs_resolve_labels <- function(res_df, gene_set_ids, labels, max_name_length) {
  raw <- res_df$Description[match(gene_set_ids, res_df$ID)]
  if (!is.null(labels) && length(labels) > 0) {
    nm <- names(labels)
    if (is.null(nm) || !any(nzchar(nm))) {
      raw <- rep(unname(as.character(labels)), length.out = length(gene_set_ids))
    } else {
      hit <- match(gene_set_ids, nm)
      raw[!is.na(hit)] <- as.character(labels)[hit[!is.na(hit)]]
    }
  }
  stats::setNames(vapply(raw, .grs_wrap_label, FUN.VALUE = character(1),
                         width = max_name_length, USE.NAMES = FALSE),
                  gene_set_ids)
}

#' Whitespace-normalise a label for MATCHING only (never for display).
#' `.grs_wrap_label()` soft-wraps long labels by turning spaces into newlines, so a
#' caller who keys a palette by the label they passed to `labels` would otherwise
#' fail to match the wrapped string the plot actually carries. Collapsing runs of
#' whitespace makes the wrapped and unwrapped forms compare equal.
.grs_norm_label <- function(x) trimws(gsub("[[:space:]]+", " ", x))

#' Resolve `palette` to a color vector NAMED BY THE FINAL PLOTTED LABEL — the
#' value the colour aesthetic actually carries — so ggplot matches BY NAME and
#' never by its alphabetically sorted level order (see fn @note).
#' Accepted keyings, in precedence order:
#'   1. already named by final label -> reordered to the caller's declared order
#'      (matched whitespace-insensitively, so pre-wrap labels still match)
#'   2. named by pathway ID          -> re-keyed to the corresponding label
#'   3. unnamed                      -> zipped to `gene_set_ids` in the order the
#'                                      caller supplied (recycled / truncated)
#' NULL keeps the internal 9-colour ramp, just label-keyed.
#' @param set_labels Output of `.grs_resolve_labels()` (named by pathway ID).
.grs_palette <- function(palette, gene_set_ids, set_labels) {
  n_sets <- length(gene_set_ids)
  labs   <- unname(set_labels[gene_set_ids])
  vals <- if (is.null(palette) || length(palette) == 0) {
    base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
    if (n_sets > length(base_pal)) grDevices::colorRampPalette(base_pal)(n_sets)
    else base_pal[seq_len(n_sets)]
  } else {
    nm <- names(palette)
    if (!is.null(nm) && all(.grs_norm_label(labs) %in% .grs_norm_label(nm))) {
      unname(palette[match(.grs_norm_label(labs),
                           .grs_norm_label(nm))])  # (1) keyed by final label
    } else if (!is.null(nm) && all(gene_set_ids %in% nm)) {
      unname(palette[gene_set_ids])               # (2) keyed by pathway ID
    } else {
      if (!is.null(nm) && any(nzchar(nm)))
        warning("`palette` names match neither the plotted labels nor the pathway IDs; ",
                "zipping colors to `gene_set_ids` order instead.")
      rep(unname(palette), length.out = n_sets)   # (3) unnamed, declared order
    }
  }
  pal <- stats::setNames(as.character(vals), labs)
  pal[!duplicated(names(pal))]   # two sets can wrap to one label -> one legend key
}

#' Set the plotted set Descriptions so the legend reads clean.
#' gseaplot2 draws the legend from `@result$Description`; we mutate a copy here so
#' `labels`/`max_name_length` take effect (rownames preserved for lookup). Long
#' names are WRAPPED onto multiple lines at `max_name_length`, never truncated with
#' an ellipsis — the full pathway name always stays readable (craft: don't truncate
#' labels). The legend grows taller, not clipped.
#' @param set_labels Output of `.grs_resolve_labels()` (named by pathway ID).
.grs_apply_labels <- function(gsea_obj, set_labels) {
  res <- gsea_obj@result
  res$Description[match(names(set_labels), res$ID)] <- unname(set_labels)
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

#' Resolve labels -> mutate the gseaResult -> build the raw panels, as ONE unit.
#' These three steps share the final label vector, which is returned alongside the
#' panels so the palette resolver (and any later restyle) keys on exactly the
#' strings the built panels carry.
.grs_raw_panels <- function(gsea_obj, res_df, gene_set_ids, labels, max_name_length, title) {
  set_labels <- .grs_resolve_labels(res_df, gene_set_ids, labels, max_name_length)
  obj        <- .grs_apply_labels(gsea_obj, set_labels)
  panels <- tryCatch(
    .grs_build_panels(obj, gene_set_ids, title),
    error = function(e) stop("Error generating running sum panels: ", conditionMessage(e)))
  list(panels = panels, set_labels = set_labels)
}

#' Build the three RAW panels. The ONLY place that knows gseaplot2's subplot
#' order (1 = ES line, 2 = gene-hit rug, 3 = ranked metric) and binds it to
#' named roles, so nothing downstream needs positional `[[i]]` indexing.
#' Deliberately does NOT forward a palette: gseaplot2's `color=` is applied
#' POSITIONALLY against alphabetically sorted `Description` levels, and is
#' silently dropped altogether unless `length(color) == length(geneSetID)`. The
#' colour scale is owned by `.grs_set_color_scale()` at style time instead.
.grs_build_panels <- function(gsea_obj, gene_set_ids, title) {
  raw <- enrichplot::gseaplot2(
    gsea_obj,
    geneSetID    = gene_set_ids,
    title        = title %||% "",
    subplots     = c(1, 2, 3),
    pvalue_table = FALSE,
    rel_heights  = c(1.5, 0.5, 0.5))
  list(es = raw[[1]], rug = raw[[2]], metric = raw[[3]])
}

#' Apply the label-named colour scale — the ONE place colours are bound.
#' `values` maps BY NAME, and `breaks = names(values)` pins LEGEND order to the
#' caller's declared order (ggplot would otherwise sort levels alphabetically).
#' This mirrors toolkit house style (`plot_standard_volcano.R`,
#' `volcano_helpers.R`, `create_MD_plot.R`, `create_fc_b_plot.R`). Replacing a
#' colour scale gseaplot2 may already have added emits ggplot2's "Scale for
#' colour is already present" message; that one message is muffled, nothing else.
.grs_set_color_scale <- function(p, palette) {
  if (is.null(palette) || length(palette) == 0 || is.null(names(palette))) return(p)
  withCallingHandlers(
    p + ggplot2::scale_color_manual(values = palette, breaks = names(palette)),
    message = function(m) {
      if (grepl("already present", conditionMessage(m), fixed = TRUE))
        invokeRestart("muffleMessage")
    })
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

#' ES (top) panel — running enrichment score. Owns the SINGLE legend, and with it
#' the label-named colour scale (values + breaks, so key colours AND key order
#' both follow the caller's declared order). We do NOT add `override.aes` (it must
#' match the key count exactly, and collides when two wrapped Descriptions map to
#' the same legend entry — a custom-DB crash).
.grs_style_es <- function(p, base_theme, base_size, es_ylim, show_x,
                          legend, legend_pos, palette = NULL) {
  p <- .grs_set_color_scale(p, palette)
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
#' guide is dropped so `guides = "collect"` yields ONE legend from the ES panel),
#' but it still needs the SAME label-named scale so each tick block matches the ES
#' curve it belongs to. Single-set exception preserved from gseaplot2: with one
#' set the rug is a plain black tick rug, not a coloured one.
.grs_style_rug <- function(p, base_theme, base_size, show_x, rug_ylabels,
                           palette = NULL) {
  if (length(palette) == 1L) palette <- stats::setNames("black", names(palette))
  p <- .grs_set_color_scale(p, palette)
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

#' Metric (bottom) panel — ranked-list statistic. Keeps the shared x-axis. Carries
#' no colour aesthetic (grey segments), so no colour scale is applied here.
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
