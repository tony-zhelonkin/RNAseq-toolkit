#' Significant-gene up/down counts tied to the highest-priority populated legend
#' category (internal helper; not part of the stable API).
#'
#' Walks the colour categories in priority order and returns the up/down split of
#' the FIRST one that has any members, so the rendered counts always belong to
#' the strongest claim the legend actually makes:
#'   1. `"both"` — combined *FDR & |log2FC|* (`sig_dec & sig_fc`). If ≥1 gene
#'      qualifies, counts are those genes split by `logFC` sign.
#'   2. `"fdr"`  — *FDR*-only (`sig_dec`). Used when priority-1 is empty; counts
#'      are the `sig_dec` genes split by sign.
#'   3. `"none"` — nothing crosses FDR; counts are `0 / 0`.
#' `legend_key` names the colour category whose legend label should receive the
#' appended count line (`"p-value & Log2FC"` for `"both"`, `"p-value"` for
#' `"fdr"`, `NA` for `"none"`). `↑` counts `logFC > 0`, `↓` counts `logFC < 0`.
#'
#' @param df Data frame carrying `logFC`, `sig_dec`, `sig_fc` (as built inside
#'           `create_standard_volcano()`).
#' @return list(category, n_up, n_down, legend_key)
#' @keywords internal
.volcano_sig_counts <- function(df) {
  if (any(df$sig_dec & df$sig_fc, na.rm = TRUE)) {
    sel <- df$sig_dec & df$sig_fc
    category   <- "both"
    legend_key <- "p-value & Log2FC"
  } else if (any(df$sig_dec, na.rm = TRUE)) {
    sel <- df$sig_dec
    category   <- "fdr"
    legend_key <- "p-value"
  } else {
    return(list(category = "none", n_up = 0L, n_down = 0L,
                legend_key = NA_character_))
  }
  list(category   = category,
       n_up       = sum(sel & df$logFC > 0, na.rm = TRUE),
       n_down     = sum(sel & df$logFC < 0, na.rm = TRUE),
       legend_key = legend_key)
}

#' Create a Standard Volcano Plot for Differential Expression Results
#'

#' ## Why raw *p* on the *y*‑axis and FDR for the decision?
#' * Raw *p*-values are per‑gene test statistics; FDR is an average over the
#'   rejected set. Plotting the basic data gives maximal resolution and avoids
#'   the stair‑step artefacts created when many *p*s map to one adjusted value.
#' * A raw‑*p* cut‑off that realises a desired FDR always exists, so we can draw
#'   a horizontal line at that *p* while keeping the more informative axis.
#'
#' The implementation therefore uses **−log10(p‑value)** on the *y*-axis but
#' decides significance with `adj.P.Val`.
#
#' ## Two alternative decision rules
#' | `decision_by` | Colour/label decision uses | Horizontal dashed line shows |
#' |---------------|----------------------------|-----------------------------|
#' | "fdr" (default)| `adj.P.Val ≤ p_cutoff`      | −log10 of the **largest raw* p*** whose *adj.P.Val* ≤ *p_cutoff* |
#' | "p"            | `P.Value ≤ p_cutoff`        | −log10(*p_cutoff*) |
#'
#' *Why keep raw‑p on the Y axis even in FDR mode?*  Raw *p* is a per‑gene test
#' statistic; FDR is an average over the rejected set. Keeping *p* preserves
#' resolution, then using FDR for the decision provides the multiple‑testing
#' guarantee while avoiding the staircase artefact of −log10(FDR).
#'
#' ### Helper utilities preserved from the original version
#' * **shade()** – darkens a colour (used for label text)
#' * **text_col()** – black/white contrast helper (currently unused but kept)
#' * **custom_minimal_theme_with_grid()** – loads user theme or falls back to
#'   `theme_minimal()`
#'
#' @param de_results  Data frame whose **rownames are gene IDs** and that
#'                    contains at least the columns `logFC`, `P.Value`,
#'                    `adj.P.Val` (limma‑style).
#' @param decision_by Character. Either "fdr" (default) or "p". Controls which
#'                    column is used for the significance call and therefore the
#'                    colour of the dots.
#' @param p_cutoff    Numeric. Interpreted as *FDR* threshold when
#'                    `decision_by = "fdr"` and as raw *p* threshold when
#'                    `decision_by = "p"`.
#' @param fc_cutoff   Numeric. Absolute log2 fold‑change threshold (default 2).
#' @param top_n       Integer. Number of top genes (ranked by the same statistic
#'                    used for the decision) to label on each side (default 5).
#' @param highlight_gene Character vector. Extra genes to label regardless of
#'                    rank (default NULL).
#' @param label_method Character. Legacy labelling behaviour retained for full
#'                    backward compatibility. "top" delegates to the new logic.
#' @param x_breaks    Numeric. Spacing between x‑axis ticks (default 1).
#' @param title       Character. Plot title (default "Volcano plot").
#' @param subtitle    Character. Optional subtitle to clarify the plot approach
#'                    (e.g., "Highlighting by FDR"). Default NULL (no subtitle).
#' @param caption     Character. Optional caption; if `NULL` a caption that
#'                    documents the thresholds will be generated.
#' @param color_palette Named character vector of four colours for the point
#'                    categories *NS*, *Log2FC*, *p‑value*, *Both*.
#' @param show_grid   Logical. If `TRUE` keep panel grid; otherwise drop it.
#' @param max.overlaps Passed to **ggrepel**.
#' @param annotate_counts Logical (default `FALSE`). When `TRUE`, append a
#'                    second line of significant‑gene up/down counts
#'                    (`↑ n   ↓ m`) directly beneath the legend label of the
#'                    highest‑priority POPULATED significance category. Priority:
#'                    the combined *FDR & |log2FC|* line first; if empty, the
#'                    *FDR*‑only line; otherwise `0 / 0` under the combined line.
#'                    `↑` counts `logFC > 0`, `↓` counts `logFC < 0`. The colour
#'                    key dot stays aligned with the threshold text (line 1) and
#'                    the counts sit on line 2, so the annotation survives any
#'                    later legend repositioning. Opt‑in; default keeps the
#'                    legend (and all existing consumers/tests) unchanged.
#' @param ...         Soft‑absorbs deprecated args such as `use_fdr`.
#'
#' @return A `ggplot2` object.
#' @export
#'
#' @examples
#' # Default FDR decision rule
#' g1 <- create_standard_volcano(de_results)
#'
#' # Classic raw‑p decision rule at p ≤ 1e‑3
#' g2 <- create_standard_volcano(de_results,
#'                               decision_by = "p",
#'                               p_cutoff    = 1e-3)
create_standard_volcano <- function(
    de_results,
    decision_by   = c("fdr", "p"),
    p_cutoff      = 0.05,
    fc_cutoff     = 2,
    top_n         = 5,
    highlight_gene= NULL,
    label_method  = "top",
    x_breaks      = 1,
    title         = "Volcano plot",
    subtitle      = NULL,
    caption       = NULL,
    fixed_p_boundary = NULL,
    color_palette = c(
      "NS"               = "#7F7F7F",   # grey
      "Log2FC"           = "#0173B2",   # blue
      "p-value"          = "#029E73",   # green
      "p-value & Log2FC" = "#D55E00"    # orange
    ),
    show_grid     = FALSE,
    max.overlaps  = 10,
    annotate_counts = FALSE,    # opt-in: counts under the top sig legend line
    ...                         # absorb deprecated args (e.g. use_fdr)
) {
  decision_by <- match.arg(decision_by)

  # ───────────────────────────────── helpers ──────────────────────────
  shade <- function(hex, factor = .6) {
    rgb <- grDevices::col2rgb(hex)/255 * factor
    grDevices::rgb(pmax(pmin(rgb,1),0)[1],
                   pmax(pmin(rgb,1),0)[2],
                   pmax(pmin(rgb,1),0)[3])
  }
  text_col <- function(hex) {
    lum <- sum(grDevices::col2rgb(hex) * c(0.299,0.587,0.114))/255
    ifelse(lum > .55, "black", "white")
  }
  custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
    source("scripts/custom_minimal_theme.R", local = TRUE)
    custom_minimal_theme_with_grid
  } else {
    function() ggplot2::theme_minimal()
  }

  # ──────────────────────── sanity checks ────────────────────────────
  stopifnot(all(c("logFC","P.Value","adj.P.Val") %in% colnames(de_results)))

  if (isTRUE(list(...)[["use_fdr"]])) {
    warning("`use_fdr` is deprecated: please use `decision_by = \"fdr\"` instead.")
  }

  # ────────────────── 1. annotate significance ───────────────────────
  if (decision_by == "fdr") {
    sig_stat   <- de_results$adj.P.Val
    stat_name  <- "FDR"
    sig_logic  <- sig_stat <= p_cutoff                # inclusive

    # Find the boundary p-value: the largest raw p among significant genes
    # This ensures the line aligns with the actual color boundary.
    # If fixed_p_boundary is provided (e.g. from a pre-filtered dataset), use it.
    if (!is.null(fixed_p_boundary)) {
      horiz_line <- -log10(fixed_p_boundary)
      draw_horiz_line <- TRUE
    } else {
      sig_pvals <- de_results$P.Value[sig_logic]
      if (length(sig_pvals) > 0) {
        p_thresh <- max(sig_pvals, na.rm = TRUE)
        horiz_line <- -log10(p_thresh)
        draw_horiz_line <- TRUE
      } else {
        # No significant genes - don't draw a line
        horiz_line <- NA
        draw_horiz_line <- FALSE
      }
    }
    legend_sig <- sprintf("FDR ≤ %.2g", p_cutoff)
  } else {  # decision_by == "p"
    sig_stat   <- de_results$P.Value
    stat_name  <- "p-value"
    sig_logic  <- sig_stat <= p_cutoff
    horiz_line <- -log10(p_cutoff)
    draw_horiz_line <- TRUE
    legend_sig <- sprintf("p ≤ %.2g", p_cutoff)
  }

  df <- dplyr::mutate(de_results,
    sig_fc  = abs(logFC) >= fc_cutoff,
    sig_dec = sig_logic,
    significance_value = sig_stat,  # Add this line
    cat = dplyr::case_when(sig_fc & sig_dec ~ "p-value & Log2FC",
                         sig_fc          ~ "Log2FC",
                         sig_dec         ~ "p-value",
                         TRUE            ~ "NS"))

  # ────────────────── 2. label selection ─────────────────────────────
    # helper: select top_n genes by the relevant significance statistic on a given
  get_top <- function(side) {
  if (side == "up") {
    df |>
      dplyr::filter(.data$logFC > 0) |>
      dplyr::arrange(.data$significance_value) |>
      dplyr::slice_head(n = top_n)
  } else {
    df |>
      dplyr::filter(.data$logFC < 0) |>
      dplyr::arrange(.data$significance_value) |>
      dplyr::slice_head(n = top_n)
  }
}

  if (label_method == "top") {
    lab_df <- dplyr::bind_rows(get_top("up"), get_top("down"))
  } else if (label_method == "sig") {
    lab_df <- df[df$cat == "p-value & Log2FC", ]
  } else if (label_method == "p") {
    lab_df <- df[df$sig_dec, ]
  } else if (label_method == "log2fc") {
    lab_df <- df[df$sig_fc, ]
  } else {
    lab_df <- df[0, ]
  }

  if (!is.null(highlight_gene)) {
    lab_df <- dplyr::bind_rows(lab_df,
                               df[rownames(df) %in% highlight_gene, ]) |>
              dplyr::distinct()
  }

  # ────────────────── 3. axis limits & colours ───────────────────────
  xmax <- ceiling(max(abs(df$logFC))/x_breaks)*x_breaks
  ymax <- ceiling(max(-log10(df$P.Value)))
  dark_pal <- vapply(color_palette, shade, character(1))

  # Legend labels (one per colour category). When annotate_counts is on, append
  # a second line of up/down counts under the highest-priority populated sig
  # line so the colour key dot stays on line 1 (the threshold) and the counts
  # sit on line 2 — anchored to that exact legend entry, not floating.
  legend_labels <- c(
    "p-value & Log2FC" = sprintf("%s & |log2FC| ≥ %.1f", legend_sig, fc_cutoff),
    "Log2FC"           = sprintf("|log2FC| ≥ %.1f", fc_cutoff),
    "p-value"          = legend_sig,
    "NS"               = "NS")

  if (isTRUE(annotate_counts)) {
    sc  <- .volcano_sig_counts(df)
    # "none" → no populated sig line; pin the 0/0 under the combined (top) line.
    key <- if (is.na(sc$legend_key)) "p-value & Log2FC" else sc$legend_key
    legend_labels[key] <- paste0(legend_labels[key],
                                 sprintf("\n↑ %d   ↓ %d", sc$n_up, sc$n_down))
  }

  # ────────────────── 4. build ggplot ────────────────────────────────
  g <- ggplot2::ggplot(df, ggplot2::aes(logFC, -log10(P.Value), colour = cat)) +
       ggplot2::geom_point(size = 2, alpha = .65) +
       ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       ggplot2::scale_colour_manual(name = NULL,
         values = color_palette,
         breaks = names(color_palette),
         labels = legend_labels) +
       ggplot2::scale_x_continuous(breaks = seq(-xmax, xmax, by = x_breaks),
                                   limits = c(-xmax, xmax)) +
       ggplot2::coord_cartesian(ylim = c(0, ymax)) +
       ggplot2::labs(x = "log2(FC)",
                     y = expression(-log[10](p-value)),
                     title = title,
                     subtitle = subtitle,
                     caption = if (is.null(caption)) {
                       if (decision_by == "fdr") {
                         if (draw_horiz_line) {
                           sprintf("Dashed lines: horiz. – FDR ≤ %.2g (p ≤ %.2g); vert. – |log2FC| ≥ %.1f",
                                   p_cutoff, signif(10^(-horiz_line),2), fc_cutoff)
                         } else {
                           sprintf("No genes pass FDR ≤ %.2g. Dashed lines: vert. – |log2FC| ≥ %.1f",
                                   p_cutoff, fc_cutoff)
                         }
                       } else {
                         sprintf("Dashed lines: horiz. – p ≤ %.2g; vert. – |log2FC| ≥ %.1f",
                                 p_cutoff, fc_cutoff)
                       }
                     } else caption) +
       custom_minimal_theme_with_grid()

  # Add horizontal line only if there are significant genes
  if (draw_horiz_line) {
    g <- g + ggplot2::geom_hline(yintercept = horiz_line, linetype = "dashed")
  } else if (decision_by == "fdr") {
    # Add annotation when no genes pass FDR threshold - positioned on the right
    g <- g + ggplot2::annotate("text",
                               x = xmax * 0.5,
                               y = ymax * 0.95,
                               label = sprintf("No genes pass FDR ≤ %.2g", p_cutoff),
                               size = 4,
                               color = "darkred",
                               fontface = "italic")
  }

  if (!show_grid) {
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                             panel.grid.minor = ggplot2::element_blank())
  }

  # ────────────────── 5. labels ───────────────────────────────────────
  if (nrow(lab_df)) {
    # Separate highlighted genes (e.g., calcium genes) from regular labels
    # This ensures priority genes are never suppressed by ggrepel overlap detection
    if (!is.null(highlight_gene) && any(rownames(lab_df) %in% highlight_gene)) {
      calcium_labs <- lab_df[rownames(lab_df) %in% highlight_gene, , drop = FALSE]
      regular_labs <- lab_df[!rownames(lab_df) %in% highlight_gene, , drop = FALSE]
    } else {
      calcium_labs <- lab_df[0, , drop = FALSE]  # Empty data frame
      regular_labs <- lab_df
    }

    # Layer 1: Regular labels (standard max.overlaps)
    if (nrow(regular_labs) > 0) {
      g <- g + ggrepel::geom_text_repel(
        data            = regular_labs,
        ggplot2::aes(label = rownames(regular_labs)),
        colour          = dark_pal[regular_labs$cat],
        fontface        = "plain",
        size            = 3.5,
        box.padding     = .4,
        point.padding   = .3,
        max.overlaps    = max.overlaps,
        min.segment.length = 0,
        show.legend     = FALSE)
    }

    # Layer 2: Highlighted genes (priority - never suppressed)
    if (nrow(calcium_labs) > 0) {
      g <- g + ggrepel::geom_text_repel(
        data            = calcium_labs,
        ggplot2::aes(label = rownames(calcium_labs)),
        colour          = "black",
        fontface        = "bold",
        size            = 3.5,
        box.padding     = .5,
        point.padding   = .3,
        max.overlaps    = Inf,     # NEVER suppress highlighted genes
        force           = 5,        # Push harder to find space
        min.segment.length = 0,
        show.legend     = FALSE)
    }
  }

  return(g)
}
