#' Plot a Heatmap for a Single GSEA Pathway
#'
#' This function creates a heatmap visualization for a specific pathway from GSEA results,
#' showing the expression patterns of core enrichment genes across samples.
#'
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param pathway_name Character, name or partial name of the pathway to visualize.
#'        The function will search both Description and ID fields.
#' @param expression_data Matrix or data frame of normalized expression values
#'        (e.g., log-CPM, log2 TPM). Rows = genes, columns = samples.
#' @param sample_order Character vector specifying the order of samples in the heatmap.
#' @param annotation_col Data frame of sample annotations for column annotation bars (default: NULL).
#' @param annotation_colors List of colors for annotation bars (default: NULL).
#' @param output_prefix Character, file path prefix for saving the PDF (default: "pathway_heatmap").
#' @param scale_expr Character, scaling parameter passed to pheatmap: "row", "column", or "none" (default: "row").
#' @param gaps_col Numeric vector of column indices where to place gaps (default: NULL).
#'
#' @return The pheatmap object (invisibly).
#' @export
#'
#' @examples
#' # Assuming gsea_obj is a GSEA result object, expression_data is a normalized expression matrix,
#' # and sample_order is a vector of sample IDs in the desired order
#' plot_single_pathway_heatmap(
#'   gsea_obj = gsea_kegg_results,
#'   pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
#'   expression_data = norm_counts,
#'   sample_order = sample_order
#' )
library(pheatmap)
library(dplyr)
library(stringr)  # for string manipulation

plot_single_pathway_heatmap <- function(
  gsea_obj,
  pathway_name,
  expression_data,
  sample_order,
  annotation_col = NULL,
  annotation_colors = NULL,
  output_prefix = "pathway_heatmap",
  scale_expr = "row",
  gaps_col = NULL  # New parameter for column gaps, NULL means no gaps
) {
  # Convert GSEA results to a data frame
  gsea_df <- as.data.frame(gsea_obj@result)
  
  # Find matching pathway
  pathway_row <- gsea_df %>%
    dplyr::filter(
      grepl(pathway_name, Description, ignore.case = TRUE) | 
      grepl(pathway_name, ID, ignore.case = TRUE)
    )
  
  if (nrow(pathway_row) == 0) {
    message("No matching pathway found for: ", pathway_name)
    return(NULL)
  }
  
  # Extract the first matching row
  pathway_info <- pathway_row[1, ]
  
  # Extract core genes
  core_genes <- unlist(strsplit(pathway_info$core_enrichment, "/"))
  core_genes <- core_genes[core_genes %in% rownames(expression_data)]
  
  if (length(core_genes) == 0) {
    message("No genes in ", pathway_name, " found in expression data.")
    return(NULL)
  }
  
  # Subset expression matrix
  expr_sub <- expression_data[core_genes, sample_order, drop = FALSE]
  
  # Clean up the pathway name for the heatmap title
  # e.g. remove prefixes like "GOBP_", "KEGG_", etc., then wrap text
  cleaned_pathway_name <- pathway_info$Description %>%
    str_replace_all(c("^GOBP_" = "", "^KEGG_" = "", "^REACTOME_" = "", "^HALLMARK_" = "")) %>%
    str_replace_all("_", " ") %>% 
    str_wrap(width = 40)  # wrap text to 40 characters per line
  
  heatmap_title <- paste("Heatmap of", cleaned_pathway_name)
  
  # Create output filename
  out_filename <- paste0(
    output_prefix, "_", 
    gsub("[^A-Za-z0-9_]+", "_", pathway_info$Description),
    ".pdf"
  )
  
  # Use pdf() and dev.off() for consistent output
  pdf(out_filename, width = 10, height = 8)
  
  # Create pheatmap parameters list
  pheatmap_params <- list(
    mat = expr_sub,
    scale = scale_expr,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("navy", "white", "red"))(50),
    main = heatmap_title,
    fontsize = 9,
    fontsize_row = 7,
    fontsize_col = 7
  )
  
  # Add gaps_col parameter only if it's not NULL
  if (!is.null(gaps_col)) {
    pheatmap_params$gaps_col <- gaps_col
  }
  
  # Plot with pheatmap
  p <- do.call(pheatmap::pheatmap, pheatmap_params)
  
  dev.off()
  
  return(invisible(p))
}
