#' Enhanced GSEA Faceted Dotplot
#' 
#' Creates a faceted dotplot showing pathways split by up/down regulation
#' with improved sizing and highlighting.
#'
#' @param gsea_obj GSEA result object
#' @param showCategory Number of pathways to show per direction
#' @param padj_cutoff Adjusted p-value cutoff
#' @param title Plot title
#' @param wrap_width Width for text wrapping
#' @param pos_color Color for positive NES
#' @param neg_color Color for negative NES
#' @param min.dotSize Minimum dot size
#' @param max.dotSize Maximum dot size
#' @param highlight_sig Whether to highlight significant points with outline
#' @param strip_prefix Logical, whether to strip common prefixes like "HALLMARK_"
#'
#' @return A ggplot2 object
#' @export
gsea_dotplot_facet <- function(
  gsea_obj,
  showCategory = 10,
  padj_cutoff = 0.05,
  title = "GSEA Faceted Dotplot",
  wrap_width = 50,
  pos_color = "#fc8d59", 
  neg_color = "#91bfdb",
  min.dotSize = 2,
  max.dotSize = 10,
  highlight_sig = TRUE,
  strip_prefix = TRUE
) {
  # Extract and filter data
  gsea_data <- as.data.frame(gsea_obj@result)
  
  # Use qvalue if present, otherwise p.adjust
  sig_col <- if("qvalue" %in% colnames(gsea_data)) "qvalue" else "p.adjust"
  
  gsea_data <- gsea_data %>%
    dplyr::filter(.data[[sig_col]] < padj_cutoff) %>%
    dplyr::mutate(
      Direction = ifelse(.data$NES > 0, "Up", "Down"),
      count = stringr::str_count(.data$core_enrichment, "/") + ifelse(nchar(.data$core_enrichment) > 0, 1, 0),
      GeneRatio = .data$count / .data$setSize,
      negLogPval = -log10(.data[[sig_col]])
    )
  
  if (nrow(gsea_data) == 0) {
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No significant pathways)")))
  }
  
  # Clean up description text
  gsea_data$Description <- stringr::str_replace_all(gsea_data$Description, "_", " ")
  
  # Strip common prefixes if requested
  if (strip_prefix) {
    common_prefixes <- c(
      "HALLMARK ", "KEGG ", "REACTOME ", "BIOCARTA ", "GOBP ", "GOCC ", "GOMF ",
      "PID ", "WIKIPATHWAY ", "^GO "
    )
    
    for (prefix in common_prefixes) {
      gsea_data$Description <- stringr::str_replace(gsea_data$Description, paste0("^", prefix), "")
    }
  }
  
  gsea_data$Description <- stringr::str_to_title(gsea_data$Description)
  gsea_data$Description <- sapply(gsea_data$Description, smart_wrap, width = wrap_width)
  
  # Get top N for each direction
  plot_data <- gsea_data %>%
    dplyr::group_by(.data$Direction) %>%
    dplyr::slice_max(order_by = abs(.data$NES), n = showCategory) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$NES)
  
  # For ordering within each facet
  plot_data <- plot_data %>%
    dplyr::group_by(.data$Direction) %>%
    dplyr::mutate(Description = factor(.data$Description, 
                                      levels = unique(.data$Description[order(.data$GeneRatio, decreasing = FALSE)]))) %>%
    dplyr::ungroup()
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = .data$Description)) +
    ggplot2::geom_point(
      ggplot2::aes(
        size = .data$negLogPval,
        color = .data$Direction
      )
    )
  
# Add outline for significant points if requested
  if (highlight_sig) {
    high_sig_threshold <- padj_cutoff / 10  # More stringent threshold for highlighting
    highlight_data <- plot_data %>% 
      dplyr::filter(.data[[sig_col]] < high_sig_threshold)
    
    if (nrow(highlight_data) > 0) {
      p <- p + 
        ggplot2::geom_point(
          data = highlight_data,
          ggplot2::aes(size = .data$negLogPval),
          shape = 21, color = "black", fill = NA, stroke = 1
        )
    }
  }
  
  # Complete the plot with scales, facets and theme
    # Complete the plot with scales, facets and theme
  y_font_size <- ifelse(nrow(plot_data) > 20, 8, 9)

  # Create pval_label before building the plot
  pval_label <- if("qvalue" %in% colnames(gsea_data)) {
    bquote(-log[10](q-value))
  } else {
    bquote(-log[10](p-value))
  }

  p <- p +
    ggplot2::scale_color_manual(values = c("Up" = pos_color, "Down" = neg_color)) +
    ggplot2::scale_size_continuous(
      name = pval_label,
      range = c(min.dotSize, max.dotSize)
    ) +
    ggplot2::facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      color = "Direction"
    ) +
    custom_minimal_theme_with_grid() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey90", color = "grey70"),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      axis.text.y = ggplot2::element_text(size = y_font_size)
    )

  
  return(p)
}
