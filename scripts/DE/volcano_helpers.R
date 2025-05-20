# R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# --------------------------------------------------
#  Vertical volcano logic with smarter caption control
# --------------------------------------------------
#' Create a *vertical* volcano plot (−log10 p on *x*, log2FC on *y*).
#'
#' ### Why vertical?
#' Turning the volcano 90° is useful when you need several panels stacked in a
#' column or grid: the fold-change axis then runs top-to-bottom, matching the
#' reader’s eye movement across rows.
#'
#' ### Caption logic
#' `caption` accepts three flavours and is **vectorised** so you can decide per
#' panel:
#' | value            | behaviour                                   |
#' |------------------|---------------------------------------------|
#' | `TRUE` (default) | auto-generate the dashed-line explanation    |
#' | `FALSE`          | no caption                                  |
#' | character string | use that string verbatim                    |
#'
#' In a `patchwork` grid keep `caption = TRUE` for the first panel and
#' `caption = FALSE` for the rest, then—if you prefer a global caption—add it
#' once with `plot_annotation(caption = …)`.
#'
#' All colour categories, labelling rules and helper utilities mirror those of
#' `create_standard_volcano()` so the figures stay visually consistent.
#'
#' @inheritParams create_standard_volcano
#' @param caption Logical/character.  See *Caption logic* table.
#' @return A `ggplot2` object.
#' @export
create_vertical_volcano <- function(
    de_results,
    decision_by   = c("fdr", "p"),
    p_cutoff      = 0.05,
    fc_cutoff     = 2,
    top_n         = 5,
    highlight_gene= NULL,
    label_method  = "top",
    x_breaks      = 1,
    title         = "Vertical Volcano Plot",
    color_palette = c(
      "NS"               = "#7F7F7F",   # grey
      "Log2FC"           = "#0173B2",   # blue
      "p-value"          = "#029E73",   # green
      "p-value & Log2FC" = "#D55E00"    # orange
    ),
    show_grid     = FALSE,
    max.overlaps  = 10
) {
  decision_by <- match.arg(decision_by)

  # ───────────────────────────────── helpers ──────────────────────────
  shade <- function(hex, factor = .6) {
    rgb <- grDevices::col2rgb(hex)/255 * factor
    grDevices::rgb(pmax(pmin(rgb,1),0)[1],
                   pmax(pmin(rgb,1),0)[2],
                   pmax(pmin(rgb,1),0)[3])
  }
  
  custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
    source("scripts/custom_minimal_theme.R", local = TRUE)
    custom_minimal_theme_with_grid
  } else {
    function() ggplot2::theme_minimal()
  }

  # ──────────────────────── sanity checks ────────────────────────────
  stopifnot(all(c("logFC","P.Value","adj.P.Val") %in% colnames(de_results)))

  # ────────────────── 1. annotate significance ───────────────────────
  if (decision_by == "fdr") {
    sig_stat   <- de_results$adj.P.Val
    stat_name  <- "FDR"
    sig_logic  <- sig_stat <= p_cutoff                # inclusive
    p_thresh   <- max(de_results$P.Value[sig_logic], p_cutoff, na.rm = TRUE)
    vert_line <- -log10(p_thresh)
    legend_sig <- sprintf("FDR ≤ %.2g", p_cutoff)
  } else {  # decision_by == "p"
    sig_stat   <- de_results$P.Value
    stat_name  <- "p-value"
    sig_logic  <- sig_stat <= p_cutoff
    vert_line <- -log10(p_cutoff)
    legend_sig <- sprintf("p ≤ %.2g", p_cutoff)
  }

  df <- dplyr::mutate(de_results,
      sig_fc  = abs(logFC) >= fc_cutoff,
      sig_dec = sig_logic,
      significance_value = sig_stat,  # Add this for sorting
      cat = dplyr::case_when(sig_fc & sig_dec ~ "p-value & Log2FC",
                             sig_fc          ~ "Log2FC",
                             sig_dec         ~ "p-value",
                             TRUE            ~ "NS"))

  # ────────────────── 2. label selection ─────────────────────────────
  # Helper: select top_n genes by the relevant significance statistic on a given side
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

  # IMPORTANT CHANGE: Create a separate dataframe for highlight genes
  if (!is.null(highlight_gene)) {
    # Get the data for highlight genes
    highlight_df <- df[rownames(df) %in% highlight_gene, ]
    
    # Remove any highlight genes that might be in the regular label dataframe
    # to avoid duplicates
    lab_df <- lab_df[!rownames(lab_df) %in% highlight_gene, ]
  } else {
    highlight_df <- NULL
  }

  # ────────────────── 3. axis limits & colours ───────────────────────
  xmax <- ceiling(max(-log10(df$P.Value), na.rm = TRUE)/x_breaks)*x_breaks
  ymax <- ceiling(max(abs(df$logFC), na.rm = TRUE))
  dark_pal <- vapply(color_palette, shade, character(1))

  # ────────────────── 4. build ggplot ────────────────────────────────
  # Key difference: x and y axes are swapped compared to standard volcano plot
  g <- ggplot2::ggplot(df, ggplot2::aes(-log10(P.Value), logFC, colour = cat)) +
       ggplot2::geom_point(size = 2, alpha = .65) +
       ggplot2::geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       ggplot2::geom_vline(xintercept = vert_line, linetype = "dashed") +
       ggplot2::scale_colour_manual(name = NULL,
         values = color_palette,
         breaks = names(color_palette),
         labels = c(
           "p-value & Log2FC" = sprintf("%s & |log2FC| ≥ %.1f", legend_sig, fc_cutoff),
           "Log2FC"           = sprintf("|log2FC| ≥ %.1f", fc_cutoff),
           "p-value"          = legend_sig,
           "NS"               = "NS")) +
       ggplot2::scale_x_continuous(breaks = seq(0, xmax, by = x_breaks),
                                  limits = c(0, xmax)) +
       ggplot2::coord_cartesian(ylim = c(-ymax, ymax)) +
       ggplot2::labs(y = "log2(FC)",
                    x = expression(-log[10](p-value)),
                    title = title,
                    caption = if (decision_by == "fdr") {
                      sprintf("Dashed lines: vert. – FDR ≤ %.2g (p ≤ %.2g); horiz. – |log2FC| ≥ %.1f",
                             p_cutoff, signif(10^(-vert_line),2), fc_cutoff)
                    } else {
                      sprintf("Dashed lines: vert. – p ≤ %.2g; horiz. – |log2FC| ≥ %.1f",
                             p_cutoff, fc_cutoff)
                    }) +
       custom_minimal_theme_with_grid()

  if (!show_grid) {
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank())
  }

  # ────────────────── 5. labels ───────────────────────────────────────
  # Regular labels with max.overlaps constraint
  if (nrow(lab_df) > 0) {
    g <- g + ggrepel::geom_text_repel(
      data            = lab_df,
      ggplot2::aes(label = rownames(lab_df)),
      colour          = dark_pal[lab_df$cat],
      size            = 3.5,
      box.padding     = .4,
      point.padding   = .3,
      max.overlaps    = max.overlaps,  # Using user-provided max.overlaps
      min.segment.length = 0,
      show.legend     = FALSE)
  }
  
  # Highlight gene labels in a separate layer with no overlap constraint
  if (!is.null(highlight_df) && nrow(highlight_df) > 0) {
    g <- g + ggrepel::geom_text_repel(
      data            = highlight_df,
      ggplot2::aes(label = rownames(highlight_df)),
      colour          = "black",  # Always black for highlight genes
      fontface        = "bold",   # Always bold for highlight genes
      size            = 3.5,
      box.padding     = .6,       # Slightly larger padding for emphasis
      point.padding   = .4,
      max.overlaps    = Inf,      # No limit on overlaps for highlight genes
      min.segment.length = 0,
      force           = 3,        # Stronger repulsion force
      show.legend     = FALSE)
  }

  return(g)
}



#' Combine Multiple Vertical Volcano Plots into a Row
#'
#' Creates a multi-panel display of vertical volcano plots with consistent scales
#' and a single unified legend.
#' 
#' @param volcano_list List of ggplot objects from create_vertical_volcano().
#' @param labels Character vector of labels for each panel (default: names of list).
#' @param guide_position Character, position for the combined legend (default: "bottom").
#'
#' @return A combined ggplot object with multiple panels and a single legend.
#' @export
#' @import ggplot2 patchwork
combine_volcano_row <- function(
    volcano_list,
    labels = names(volcano_list),
    guide_position = "bottom",
    keep_first_caption = FALSE          # set TRUE if you really want caption under panel 1
) {
  ## 1. global limits ---------------------------------------------------
  global_y <- max(sapply(volcano_list, function(p)
    max(abs(ggplot2::layer_scales(p)$y$range$range), na.rm = TRUE)))
  global_x <- max(sapply(volcano_list, function(p)
    max(ggplot2::layer_scales(p)$x$range$range, na.rm = TRUE)))

  ## 2. capture caption from first plot ---------------------------------
  # Fixed: using ggplot_build instead of plot_build
  first_cap <- ggplot2::ggplot_build(volcano_list[[1]])$layout$plot$labels$caption
  if (is.null(first_cap) || !nzchar(first_cap)) first_cap <- NULL

  ## 3. strip captions, unify titles & limits ---------------------------
  for (i in seq_along(volcano_list)) {
    volcano_list[[i]] <- volcano_list[[i]] +
      ggplot2::coord_cartesian(xlim = c(0, global_x),
                               ylim = c(-global_y, global_y),
                               clip = "off") +
      ggplot2::ggtitle(labels[i]) +
      ggplot2::labs(caption = if (keep_first_caption && i == 1) first_cap else NULL)
  }

  ## 4. assemble, collect guides, add one caption -----------------------
  combined <- patchwork::wrap_plots(volcano_list, nrow = 1) +
              patchwork::plot_layout(guides = "collect") &
              ggplot2::theme(legend.position = guide_position)

  if (!is.null(first_cap) && !keep_first_caption) {
    combined <- combined & patchwork::plot_annotation(caption = first_cap)
  }

  combined
}
