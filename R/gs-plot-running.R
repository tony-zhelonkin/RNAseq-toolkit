#' GSEA running-sum (enrichment curve) plot
#'
#' The canonical three-panel GSEA figure, rebuilt from scratch on
#' [fgsea::plotEnrichmentData()]:
#'
#' 1. the running enrichment score (ES) curve,
#' 2. one lane of gene-hit ticks per pathway,
#' 3. the ranked metric itself.
#'
#' # Why this is one faceted ggplot, not three composed plots
#'
#' The old `gsea_running_sum_plot()` built three separate ggplots and glued
#' them with `patchwork`, then spent 16 helpers and three `@note` blocks
#' defending the alignment and colour invariants that composition kept
#' breaking. Here the three panels are **rows of a single
#' [ggplot2::facet_grid()]**, so they share one x scale *by construction* --
#' alignment is not something a composer has to preserve, and there is no
#' `patchwork` dependency. The return value is a plain `ggplot`, which is what
#' every other `gs_plot_*` renderer returns and what `gs_save()` expects.
#'
#' # Why colours cannot be permuted
#'
#' The colour aesthetic is mapped to **`pathway_id`**, never to a display
#' label, and the scale is an explicit
#' `scale_colour_manual(values = <named vector>, breaks = names(values))` keyed
#' by that same id. Legend *text* is applied through the scale's `labels`, so
#' renaming a pathway can never move its colour. The old failure mode -- an
#' unnamed palette handed to `enrichplot::gseaplot2(color = )`, which assigns
#' positionally while ggplot orders levels alphabetically, silently swapping
#' two colours while the name-keyed legend text stayed correct -- is
#' structurally impossible here.
#'
#' # Panel heights
#'
#' `panel_heights` is the single knob for panel proportions, and it is honoured
#' exactly: each panel's y range is padded (never clipped) so that the ranges
#' stand in the requested ratio, and `facet_grid(space = "free_y")` sizes the
#' panels from those ranges. No `aspect.ratio` is set anywhere, so panel widths
#' cannot desynchronise.
#'
#' @param x A [gs_result][gs_result-class] (the usual case), a `gs_db`, or a
#'   named list of character vectors. A `gs_result` supplies pathway ids,
#'   display names and the `stat` used to pick the default pathways; it does
#'   **not** carry gene-set membership, so `db` is required with it.
#' @param ranks Named numeric vector of gene-level statistics, decreasing, as
#'   returned by [gs_ranks()]. If `NULL`, `attr(x, "ranks")` is used -- the hook
#'   a deprecation shim can fill in.
#' @param db A `gs_db` or named list of character vectors giving set
#'   membership. If `NULL`, `attr(x, "gene_sets")` is used, or `x` itself when
#'   `x` is already a `gs_db` / named list.
#' @param pathways Character vector of pathway ids, or (for a `gs_result`)
#'   integer row indices. `NULL` picks the top `top` pathways by `abs(stat)`
#'   for a `gs_result`, or the first `top` sets otherwise. The order given is
#'   the order colours and legend keys are assigned in.
#' @param top Integer. How many pathways to pick when `pathways` is `NULL`.
#' @param labels Optional character vector of legend labels, named by pathway
#'   id (a partial map degrades gracefully -- unnamed ids keep their
#'   `pathway_name`), or unnamed and zipped to `pathways` in the given order.
#' @param palette Optional colours, named by pathway id or unnamed and zipped
#'   to `pathways` in the given order (recycled). `NULL` uses a nine-colour
#'   vibrant ramp.
#' @param panel_heights Length-3 numeric, the ES : ticks : metric height ratio.
#' @param gsea_param Numeric exponent passed to
#'   [fgsea::plotEnrichmentData()]'s `gseaParam`.
#' @param metric_label Character. Strip label for the ranked-metric panel; name
#'   the statistic you ranked by (e.g. `"t statistic"`, `"log2 FC"`).
#' @param title Optional plot title.
#' @param base_size Base font size, floored at 14 to stay legible at
#'   journal-column width.
#' @param max_name_length Integer. Legend labels longer than this are *wrapped*
#'   onto several lines, never truncated.
#' @param legend_position One of `"right"` (default), `"bottom"`, `"inside"`
#'   or `"none"`.
#' @param legend_pos Length-2 numeric, the inside-legend position in npc units,
#'   used only when `legend_position = "inside"`.
#' @param base_theme Optional complete [ggplot2::theme()] used as the
#'   foundation; the panel chrome is re-applied on top so a project theme
#'   cannot clobber the layout. `NULL` uses `theme_bulki()` when available and
#'   [ggplot2::theme_minimal()] otherwise.
#' @return A `ggplot`.
#' @seealso [gs_ranks()], [gs_test()], [gs_leading_edge()]
#' @importFrom rlang .data
#' @importFrom fgsea plotEnrichmentData
#' @examples
#' sets <- list(SET_A = c("A", "B", "C"), SET_B = c("D", "E", "F"))
#' ranks <- stats::setNames(seq(3, -3, length.out = 8), LETTERS[1:8])
#' gs_plot_running(sets, ranks = ranks)
#' @export
gs_plot_running <- function(x,
                            ranks = NULL,
                            db = NULL,
                            pathways = NULL,
                            top = 5L,
                            labels = NULL,
                            palette = NULL,
                            panel_heights = c(2.4, 0.7, 0.9),
                            gsea_param = 1,
                            metric_label = "Ranked metric",
                            title = NULL,
                            base_size = 14,
                            max_name_length = 40,
                            legend_position = c("right", "bottom", "inside",
                                                "none"),
                            legend_pos = c(0.98, 0.98),
                            base_theme = NULL) {
  legend_position <- match.arg(legend_position)
  if (!is.numeric(panel_heights) || length(panel_heights) != 3L ||
        any(!is.finite(panel_heights)) || any(panel_heights <= 0)) {
    stop("`panel_heights` must be three positive finite numbers ",
         "(ES : ticks : metric).", call. = FALSE)
  }

  ranks <- .grs_ranks(x, ranks)
  sets <- .grs_sets(x, db)
  ids <- .grs_select(x, sets, pathways, top)
  set_labels <- .grs_labels(x, ids, labels, max_name_length)
  pal <- .grs_palette(palette, ids)

  panels <- .grs_panel_levels(metric_label)
  curves <- .grs_curves(sets[ids], ranks, gsea_param)

  es_df <- do.call(rbind, lapply(ids, function(id) {
    d <- curves[[id]]$curve
    data.frame(rank = d$rank, y = d$ES, pathway_id = id,
               stringsAsFactors = FALSE)
  }))
  tick_df <- do.call(rbind, lapply(seq_along(ids), function(i) {
    d <- curves[[ids[i]]]$ticks
    n <- length(ids)
    data.frame(rank = d$rank, ymin = n - i + 0.2, ymax = n - i + 0.8,
               pathway_id = ids[i], stringsAsFactors = FALSE)
  }))
  stat_df <- curves[[ids[1L]]]$stats
  stat_df <- data.frame(rank = stat_df$rank, y = stat_df$stat,
                        stringsAsFactors = FALSE)

  win <- .grs_windows(
    list(es = range(c(0, es_df$y)),
         ticks = c(0, length(ids)),
         stats = range(c(0, stat_df$y))),
    panel_heights
  )
  pad_df <- do.call(rbind, lapply(names(win), function(k) {
    data.frame(panel = k, rank = 1, y = c(0, win[[k]]$h),
               stringsAsFactors = FALSE)
  }))
  pad_df$panel <- factor(pad_df$panel, levels = names(panels))

  es_df$y <- .grs_map(es_df$y, win$es)
  tick_df$ymin <- .grs_map(tick_df$ymin, win$ticks)
  tick_df$ymax <- .grs_map(tick_df$ymax, win$ticks)
  stat_df$y <- .grs_map(stat_df$y, win$stats)
  es_df$panel <- factor("es", levels = names(panels))
  tick_df$panel <- factor("ticks", levels = names(panels))
  stat_df$panel <- factor("stats", levels = names(panels))
  stat_df$base <- .grs_map(0, win$stats)
  zero_df <- data.frame(
    panel = factor(c("es", "stats"), levels = names(panels)),
    y = c(.grs_map(0, win$es), .grs_map(0, win$stats))
  )

  p <- ggplot() +
    geom_blank(
      data = pad_df,
      mapping = aes(x = .data$rank, y = .data$y)
    ) +
    geom_hline(
      data = zero_df,
      mapping = aes(yintercept = .data$y),
      colour = "grey60", linewidth = 0.3
    ) +
    geom_ribbon(
      data = stat_df,
      mapping = aes(x = .data$rank, ymin = .data$base, ymax = .data$y),
      fill = "grey70", colour = "transparent"
    ) +
    geom_segment(
      data = tick_df,
      mapping = aes(x = .data$rank, xend = .data$rank,
                    y = .data$ymin, yend = .data$ymax,
                    colour = .data$pathway_id),
      linewidth = 0.35
    ) +
    geom_line(
      data = es_df,
      mapping = aes(x = .data$rank, y = .data$y, colour = .data$pathway_id),
      linewidth = 0.7
    ) +
    scale_colour_manual(
      values = pal,
      breaks = names(pal),
      labels = unname(set_labels[names(pal)]),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0)),
      breaks = .grs_break_fun(win)
    ) +
    facet_grid(rows = vars(.data$panel), scales = "free_y",
               space = "free_y", switch = "y",
               labeller = as_labeller(panels)) +
    labs(x = "Rank in ranked gene list", y = NULL, title = title)

  p + .grs_theme(base_theme, base_size, legend_position, legend_pos)
}

# ---- internals (prefix `.grs_`) ---------------------------------------------

#' Resolve the rank vector for [gs_plot_running()]
#'
#' @param x The renderer's `x`.
#' @param ranks The renderer's `ranks`, possibly `NULL`.
#' @return A named numeric vector, sorted decreasing.
#' @keywords internal
.grs_ranks <- function(x, ranks) {
  ranks <- ranks %||% attr(x, "ranks")
  if (is.null(ranks)) {
    stop("`ranks` is required: supply the named numeric vector from ",
         "`gs_ranks()`. A `gs_result` does not carry the ranked list.",
         call. = FALSE)
  }
  if (!is.numeric(ranks) || is.null(names(ranks))) {
    stop("`ranks` must be a *named* numeric vector; see `gs_ranks()`.",
         call. = FALSE)
  }
  sort(ranks, decreasing = TRUE)
}

#' Resolve gene-set membership for [gs_plot_running()]
#'
#' @param x The renderer's `x`.
#' @param db The renderer's `db`, possibly `NULL`.
#' @return A named list of character vectors.
#' @keywords internal
.grs_sets <- function(x, db) {
  cand <- db %||% attr(x, "gene_sets")
  if (is.null(cand) && !inherits(x, "gs_result") && .grs_is_sets(x)) {
    cand <- x
  }
  if (is.null(cand)) {
    stop("`db` is required: the running curve needs full gene-set membership, ",
         "which a `gs_result` does not carry. Pass the `gs_db` you tested ",
         "against.", call. = FALSE)
  }
  if (!.grs_is_sets(cand)) {
    stop("`db` must be a `gs_db` or a named list of character vectors.",
         call. = FALSE)
  }
  sets <- unclass(cand)
  attributes(sets) <- list(names = names(sets))
  lapply(sets, as.character)
}

#' Is `x` shaped like a named list of gene sets?
#'
#' @param x Object to test.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
.grs_is_sets <- function(x) {
  is.list(x) && length(x) > 0L && !is.null(names(x)) &&
    all(nzchar(names(x))) &&
    all(vapply(x, function(z) is.character(z) || is.factor(z), logical(1L)))
}

#' Choose which pathways to draw
#'
#' @param x The renderer's `x`.
#' @param sets Named list of gene sets.
#' @param pathways Ids, row indices, or `NULL`.
#' @param top Integer default count.
#' @return A character vector of pathway ids, in plotting order.
#' @keywords internal
.grs_select <- function(x, sets, pathways, top) {
  if (is.null(pathways)) {
    if (inherits(x, "gs_result")) {
      if (!nrow(x)) {
        stop("`x` has no rows, so there is nothing to plot.", call. = FALSE)
      }
      ord <- order(abs(x$stat), decreasing = TRUE)
      ids <- x$pathway_id[ord][seq_len(min(top, length(ord)))]
    } else {
      ids <- names(sets)[seq_len(min(top, length(sets)))]
    }
  } else if (is.numeric(pathways)) {
    if (!inherits(x, "gs_result")) {
      stop("Integer `pathways` index the rows of a `gs_result`; pass pathway ",
           "ids instead.", call. = FALSE)
    }
    bad <- pathways < 1 | pathways > nrow(x)
    if (any(bad)) {
      stop("`pathways` index outside the ", nrow(x), " rows of `x`: ",
           paste(pathways[bad], collapse = ", "), ".", call. = FALSE)
    }
    ids <- x$pathway_id[as.integer(pathways)]
  } else {
    ids <- as.character(pathways)
  }
  ids <- unique(ids)
  missing_ids <- setdiff(ids, names(sets))
  if (length(missing_ids)) {
    stop("Pathway(s) not found in `db`: ",
         paste(sQuote(missing_ids), collapse = ", "), ".", call. = FALSE)
  }
  ids
}

#' Legend labels, keyed by pathway id
#'
#' @param x The renderer's `x`.
#' @param ids Character vector of pathway ids.
#' @param labels The renderer's `labels`.
#' @param max_name_length Wrap width.
#' @return A character vector of labels, named by pathway id.
#' @keywords internal
.grs_labels <- function(x, ids, labels, max_name_length) {
  base <- stats::setNames(ids, ids)
  if (inherits(x, "gs_result")) {
    hit <- match(ids, x$pathway_id)
    ok <- !is.na(hit)
    base[ok] <- x$pathway_name[hit[ok]]
  } else {
    pn <- attr(x, "pathway_names")
    if (!is.null(pn)) {
      hit <- pn[ids]
      base[!is.na(hit)] <- unname(hit[!is.na(hit)])
    }
  }
  if (!is.null(labels)) {
    if (!is.null(names(labels))) {
      hit <- intersect(ids, names(labels))
      base[hit] <- as.character(labels[hit])
    } else {
      labels <- rep(as.character(labels), length.out = length(ids))
      base[] <- labels
    }
  }
  vapply(base, .grs_wrap, character(1L), width = max_name_length,
         USE.NAMES = TRUE)
}

#' Wrap a label onto several lines rather than truncating it
#'
#' @param s Character scalar.
#' @param width Maximum line width.
#' @return A character scalar, possibly containing newlines.
#' @keywords internal
.grs_wrap <- function(s, width) {
  if (is.na(s) || !nzchar(s)) return("")
  paste(strwrap(s, width = max(width, 8L)), collapse = "\n")
}

#' Resolve a palette keyed by pathway id
#'
#' Named palettes are matched by id; unnamed ones are zipped to `ids` in the
#' caller's declared order. The result is always named by id, which is what the
#' colour aesthetic maps, so ggplot's alphabetical level order cannot permute
#' the colours.
#'
#' @param palette The renderer's `palette`.
#' @param ids Character vector of pathway ids, in plotting order.
#' @return A character vector of colours, named by pathway id.
#' @keywords internal
.grs_palette <- function(palette, ids) {
  n <- length(ids)
  if (is.null(palette) || length(palette) == 0L) {
    base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                  "#FFFF33", "#A65628", "#F781BF", "#999999")
    vals <- if (n > length(base_pal)) {
      grDevices::colorRampPalette(base_pal)(n)
    } else {
      base_pal[seq_len(n)]
    }
    return(stats::setNames(vals, ids))
  }
  nm <- names(palette)
  if (!is.null(nm) && all(ids %in% nm)) {
    return(stats::setNames(as.character(palette[ids]), ids))
  }
  if (!is.null(nm) && any(nzchar(nm))) {
    warning("`palette` names do not cover every plotted pathway id; zipping ",
            "colours to `pathways` order instead.", call. = FALSE)
  }
  stats::setNames(rep(as.character(palette), length.out = n), ids)
}

#' Panel keys and their strip labels
#'
#' @param metric_label Strip label for the ranked-metric panel.
#' @return A named character vector, `panel key -> strip label`.
#' @keywords internal
.grs_panel_levels <- function(metric_label = "Ranked metric") {
  if (!is.character(metric_label) || length(metric_label) != 1L) {
    stop("`metric_label` must be a single string.", call. = FALSE)
  }
  c(es = "Enrichment score", ticks = "Genes", stats = metric_label)
}

#' Running-curve data for each pathway
#'
#' Delegates to [fgsea::plotEnrichmentData()], which returns `curve`, `ticks`
#' and `stats` -- exactly the three panels. The cumulative sum is never
#' recomputed here.
#'
#' @param sets Named list of gene sets to draw.
#' @param ranks Named numeric vector, decreasing.
#' @param gsea_param `gseaParam` for fgsea.
#' @return A named list of `plotEnrichmentData()` outputs.
#' @keywords internal
.grs_curves <- function(sets, ranks, gsea_param) {
  out <- lapply(names(sets), function(id) {
    genes <- intersect(sets[[id]], names(ranks))
    if (!length(genes)) {
      stop("Pathway ", sQuote(id), " has no genes in `ranks`, so it has no ",
           "running curve.", call. = FALSE)
    }
    fgsea::plotEnrichmentData(pathway = genes, stats = ranks,
                              gseaParam = gsea_param)
  })
  stats::setNames(out, names(sets))
}

#' Panel windows: exact panel heights without clipping any data
#'
#' `facet_grid(space = "free_y")` sizes a panel in proportion to its y range,
#' and the natural ranges of an ES curve, a tick lane and a ranked metric have
#' nothing to do with the proportions a reader wants. So each panel's data is
#' mapped onto its own window `[0, panel_heights[k]]` -- a linear, invertible
#' rescaling that clips nothing -- and the panel heights then come out exactly in
#' the requested ratio. The axis still reads in original units, because
#' [.grs_break_fun()] inverts the map and returns *named* breaks, which ggplot2
#' uses as the labels.
#'
#' A 10% margin is added to the ES and metric ranges here rather than left to
#' scale expansion, so the plotted window limits are exactly `c(0, h)` and each
#' panel stays identifiable by them.
#'
#' @param ranges Named list of length-2 numeric natural ranges, in panel order.
#' @param heights Length-3 numeric ratio.
#' @return A named list, one entry per panel, each `list(lo =, hi =, h =)`.
#' @keywords internal
.grs_windows <- function(ranges, heights) {
  heights <- as.numeric(heights) / sum(heights)
  # Panels are told apart by their window height, so those must be distinct.
  for (i in seq_along(heights)) {
    while (any(abs(heights[-i] - heights[i]) < 1e-9)) {
      heights[i] <- heights[i] + 1e-6
    }
  }
  out <- lapply(seq_along(ranges), function(i) {
    r <- as.numeric(ranges[[i]])
    if (names(ranges)[i] != "ticks") {
      pad <- max(diff(r), .Machine$double.eps) * 0.05
      r <- c(r[1L] - pad, r[2L] + pad)
    }
    if (diff(r) <= 0) r <- r + c(-0.5, 0.5)
    list(lo = r[1L], hi = r[2L], h = heights[[i]])
  })
  stats::setNames(out, names(ranges))
}

#' Map data values onto a panel window
#'
#' @param y Numeric vector in original units.
#' @param w A window from [.grs_windows()].
#' @return `y` rescaled onto `[0, w$h]`.
#' @keywords internal
.grs_map <- function(y, w) {
  (y - w$lo) / (w$hi - w$lo) * w$h
}

#' Break function that restores original units panel by panel
#'
#' With `scales = "free_y"` every panel resolves its own breaks and a break
#' function gets no panel context -- but each panel's limits are exactly
#' `c(0, h)` (pinned by `geom_blank()` with zero expansion, heights made
#' distinct in [.grs_windows()]), so the panel is identifiable from its limits
#' alone. Breaks are chosen in original units, mapped into window space, and
#' returned *named* with their original-unit labels. The tick panel's y index
#' is meaningless, so it gets no breaks at all.
#'
#' @param win The output of [.grs_windows()].
#' @return A function suitable for `scale_y_continuous(breaks = )`.
#' @keywords internal
.grs_break_fun <- function(win) {
  function(limits) {
    h <- max(as.numeric(limits))
    key <- names(win)[which.min(vapply(win, function(w) abs(w$h - h),
                                       numeric(1L)))]
    if (identical(key, "ticks")) return(numeric(0))
    w <- win[[key]]
    b <- scales::extended_breaks()(c(w$lo, w$hi))
    b <- b[b >= w$lo & b <= w$hi]
    stats::setNames(.grs_map(b, w), format(b, trim = TRUE))
  }
}

#' Panel theme and legend placement
#'
#' @param base_theme Optional complete theme used as the foundation.
#' @param base_size Base font size, floored at 14.
#' @param legend_position One of `"right"`, `"bottom"`, `"inside"`, `"none"`.
#' @param legend_pos Length-2 numeric inside-legend position.
#' @return A list of ggplot2 theme objects, applied in order.
#' @keywords internal
.grs_theme <- function(base_theme, base_size, legend_position, legend_pos) {
  base_size <- max(base_size, 14)
  base <- base_theme %||% .grs_default_theme(base_size)
  chrome <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90),
    panel.spacing.y = unit(2, "pt"),
    legend.position = if (legend_position == "inside") {
      "inside"
    } else {
      legend_position
    },
    legend.position.inside = legend_pos,
    legend.justification = if (legend_position == "inside") {
      c(1, 1)
    } else {
      "center"
    },
    legend.key.spacing.y = unit(2, "pt")
  )
  list(base, chrome)
}

#' Default foundation theme
#'
#' Uses `theme_bulki()` when the package provides it (B3's renderer theme) and
#' falls back to [ggplot2::theme_minimal()] so this file has no copy of it.
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme.
#' @keywords internal
.grs_default_theme <- function(base_size) {
  tb <- get0("theme_bulki", envir = asNamespace("bulkiRNA"),
             mode = "function", ifnotfound = NULL)
  if (!is.null(tb)) tb(base_size = base_size) else theme_minimal(base_size)
}
