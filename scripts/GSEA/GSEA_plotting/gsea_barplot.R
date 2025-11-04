#' Enhanced GSEA Barplot
#' 
#' Creates a horizontal barplot of NES values with cleaner pathway names.
#'
#' @param gsea_obj GSEA result object
#' @param padj_cutoff Adjusted p-value cutoff
#' @param top_n Number of pathways to show
#' @param title Plot title
#' @param pos_color Color for positive NES
#' @param neg_color Color for negative NES
#' @param strip_prefix Logical, whether to strip common prefixes like "HALLMARK_"
#'
#' @return A ggplot2 object
#' @export
#'
#' @note Requires format_pathway_name() function to be available in environment.
#'       This is typically sourced by run_gsea_analysis() before calling this function.

gsea_barplot <- function(
  gsea_obj,
  padj_cutoff = 0.05,
  top_n = 30,
  title = "GSEA NES Barplot",
  pos_color = "#fc8d59",
  neg_color = "#91bfdb",
  strip_prefix = TRUE
) {
  # Extract and filter data
  # Ensure padj_cutoff is numeric
  padj_cutoff_num <- as.numeric(padj_cutoff)
  if (is.na(padj_cutoff_num)) {
    warning("padj_cutoff is not numeric, using default value of 0.05")
    padj_cutoff_num <- 0.05
  }
  
  gsea_data <- as.data.frame(gsea_obj@result) %>%
    dplyr::filter(.data$p.adjust < padj_cutoff_num)
  
  if (nrow(gsea_data) == 0) {
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No significant pathways)")))
  }
  
  # Sort by absolute NES and take top N
  gsea_data <- gsea_data %>%
    dplyr::arrange(dplyr::desc(abs(.data$NES))) %>%
    utils::head(top_n)
  
  # Format pathway names using smart capitalization with biological exceptions
  gsea_data$Description <- format_pathway_name(
    gsea_data$Description,
    use_formatting = TRUE,
    strip_prefix = strip_prefix
  )
  
  # Sort by NES for display
  gsea_data <- gsea_data %>% dplyr::arrange(.data$NES)
  
  # Create plot
  p <- ggplot2::ggplot(gsea_data, ggplot2::aes(x = stats::reorder(.data$Description, .data$NES), y = .data$NES)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = .data$NES > 0)) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = neg_color, `TRUE` = pos_color)) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = "NES"
    ) +
    custom_minimal_theme_with_grid() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank()
    )
  
  return(p)
}
