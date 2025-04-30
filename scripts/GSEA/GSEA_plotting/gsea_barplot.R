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
  gsea_data <- as.data.frame(gsea_obj@result) %>%
    dplyr::filter(.data$p.adjust < padj_cutoff)
  
  if (nrow(gsea_data) == 0) {
    return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No significant pathways)")))
  }
  
  # Sort by absolute NES and take top N
  gsea_data <- gsea_data %>%
    dplyr::arrange(dplyr::desc(abs(.data$NES))) %>%
    utils::head(top_n)
  
  # Clean up description text
  gsea_data$Description <- stringr::str_replace_all(gsea_data$Description, "_", " ")
  
  # Strip common prefixes like "HALLMARK ", "GOBP ", etc.
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
