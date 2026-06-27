#' Enhanced GSEA Faceted Dotplot with Continuous NES Gradient
#'
#' Creates a faceted dotplot showing pathways split by up/down regulation.
#' This function now shares the same rendering contract as `gsea_dotplot()`:
#'   1. the base dots use a filled-circle glyph (`shape = 21`) with NES mapped
#'      to `fill` via the same continuous diverging scale,
#'   2. significance is indicated only by a black outline overlay,
#'   3. pathway selection is driven by significance (`p.adjust`/`qvalue`), with
#'      the top `showCategory` pathways chosen WITHIN each direction facet.
#'
#' @param gsea_obj GSEA result object
#' @param showCategory Number of pathways to show per direction
#' @param padj_cutoff Adjusted p-value cutoff used for significance highlighting
#' @param title Plot title
#' @param wrap_width Width for text wrapping
#' @param neg_color Color for negative NES (default: colorblind-safe blue #2166AC)
#' @param mid_color Color for zero NES (default: white #F7F7F7)
#' @param pos_color Color for positive NES (default: colorblind-safe orange #B35806)
#' @param nes_limits NES color scale limits (default: c(-3.5, 3.5))
#' @param min.dotSize Minimum dot size
#' @param max.dotSize Maximum dot size
#' @param highlight_sig Whether to highlight significant points with outline
#' @param highlight_threshold FDR threshold for highlighting significant points.
#'        If NULL (default), uses `padj_cutoff`. Set explicitly to override.
#' @param strip_prefix Logical, whether to strip common prefixes like "HALLMARK_"
#'
#' @return A ggplot2 object
#' @export
#'
#' @note Requires format_pathway_name() function to be available in environment.
#'       This is typically sourced by run_gsea_analysis() before calling this function.
#'
#' @note Rendering/selection contract updated 2026-06-26 to match
#'       `gsea_dotplot()`: fill-based NES gradient, outline-only significance,
#'       and top-by-significance selection within each Up/Down facet.

facet_grid_with_left_border <- function(...) {
  facet <- ggplot2::facet_grid(...)

  facet$params$strip.background.y <- list(
    ggplot2::element_rect(color = "black", fill = NA, size = 1.5, linewidth = 1.5,
                          linetype = "solid", inherit.blank = FALSE)
  )

  facet
}


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
  strip_prefix = TRUE
) {
  gsea_data <- as.data.frame(gsea_obj@result)

  if (is.null(gsea_data) || nrow(gsea_data) == 0) {
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No pathways)")))
  }

  # Use qvalue if present, otherwise p.adjust
  sig_col <- if ("qvalue" %in% colnames(gsea_data)) "qvalue" else "p.adjust"

  gsea_data <- gsea_data %>%
    dplyr::mutate(
      Direction = ifelse(.data$NES > 0, "Up", "Down"),
      count = stringr::str_count(.data$core_enrichment, "/") + ifelse(nchar(.data$core_enrichment) > 0, 1, 0),
      GeneRatio = .data$count / .data$setSize,
      negLogPval = -log10(.data[[sig_col]])
    )

  # Format pathway names using smart capitalization with biological exceptions
  gsea_data$Description <- format_pathway_name(
    gsea_data$Description,
    use_formatting = TRUE,
    strip_prefix = strip_prefix
  )

  # Apply text wrapping
  gsea_data$Description <- sapply(gsea_data$Description, smart_wrap, width = wrap_width)

  # Select top-N by significance WITHIN each direction facet (same display contract as gsea_dotplot)
  plot_data <- gsea_data %>%
    dplyr::group_by(.data$Direction) %>%
    dplyr::arrange(.data[[sig_col]], .by_group = TRUE) %>%
    dplyr::slice_head(n = showCategory) %>%
    dplyr::ungroup()

  if (nrow(plot_data) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste(title, "(No pathways after selection)")))
  }

  # Reorder within each facet for y-axis display (same GeneRatio-first ordering as gsea_dotplot)
  plot_data <- plot_data %>%
    dplyr::group_by(.data$Direction) %>%
    dplyr::arrange(dplyr::desc(.data$GeneRatio), .by_group = TRUE) %>%
    dplyr::mutate(Description = factor(.data$Description,
                                       levels = rev(unique(.data$Description)))) %>%
    dplyr::ungroup()

  # Base layer: filled-circle dots with NES mapped to fill (same as gsea_dotplot)
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = .data$Description)) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$negLogPval,
        fill = .data$NES
      ),
      shape = 21,
      stroke = 0,
      color = "transparent"
    )

  # Add outline for significant points if requested
  if (highlight_sig) {
    if (!is.null(highlight_threshold)) {
      high_sig_threshold <- as.numeric(highlight_threshold)
    } else {
      padj_cutoff_num <- as.numeric(padj_cutoff)
      if (is.na(padj_cutoff_num)) {
        warning("padj_cutoff is not numeric, using default value of 0.05")
        padj_cutoff_num <- 0.05
      }
      high_sig_threshold <- padj_cutoff_num
    }

    highlight_data <- plot_data %>%
      dplyr::filter(.data[[sig_col]] < high_sig_threshold)

    if (nrow(highlight_data) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = highlight_data,
          ggplot2::aes(size = .data$negLogPval),
          shape = 21, color = "black", fill = NA, stroke = 2
        )
    }
  }

  y_font_size <- ifelse(nrow(plot_data) > 20, 8, 9)

  pval_label <- if ("qvalue" %in% colnames(gsea_data)) {
    bquote(-log[10](q-value))
  } else {
    bquote(-log[10](p-value))
  }

  p +
    ggplot2::scale_fill_gradient2(
      low = neg_color,
      mid = mid_color,
      high = pos_color,
      midpoint = 0,
      limits = nes_limits,
      oob = scales::squish,
      name = "NES"
    ) +
    ggplot2::scale_size_continuous(
      name = pval_label,
      range = c(min.dotSize, max.dotSize)
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        override.aes = list(
          shape = 16,
          fill = "black",
          color = "black"
        )
      )
    ) +
    facet_grid_with_left_border(Direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = title,
      x = "Gene Ratio",
      y = NULL
    ) +
    custom_minimal_theme_with_grid() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "plain"),
      strip.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 5, l = 5)),
      panel.spacing.y = ggplot2::unit(1, "lines"),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      axis.text.y = ggplot2::element_text(size = y_font_size)
    )
}
