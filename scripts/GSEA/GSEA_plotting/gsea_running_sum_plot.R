#' Create and Save a GSEA Running Sum Enrichment Plot
#' 
#' This function creates and saves a GSEA enrichment plot showing the running sum
#' statistic for specified gene sets, using the enrichplot package.
#' 
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param gene_set_ids Vector of gene set IDs to plot (numeric indices or character IDs).
#' @param title Character, plot title and base filename for saving (required).
#' @param subplots Numeric vector, indicating which subplots to include (default: c(1,2,3)):
#'        1 = running enrichment score, 2 = positions of gene set members, 3 = ranking metric scores.
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 10).
#' @param height Numeric, height of the saved plot in inches (default: 5).
#' @param dpi Numeric, resolution of the saved plot (default: 300).
#'
#' @return The generated plot object (invisibly).
#' @export
#'
#' @examples
#' # Basic usage
#' gsea_running_sum_plot(gsea_obj, gene_set_ids = 1:5, 
#'                      title = "HALLMARK running sum for top 5 pathways")
#'
#' # Save the plot
#' gsea_running_sum_plot(gsea_obj, gene_set_ids = 1:3, 
#'                      title = "Top pathways", save_plot = TRUE)
library(enrichplot)
library(ggplot2)

gsea_running_sum_plot <- function(gsea_obj, 
                                 gene_set_ids,
                                 title,
                                 subplots = c(1, 2, 3),
                                 save_plot = FALSE,
                                 output_dir = "plots/",
                                 width = 10,
                                 height = 5,
                                 dpi = 300) {
    
    # Input validation
    if (missing(title) || is.null(title) || title == "") {
        stop("Title must be provided")
    }
    
    # Create the GSEA plot
    p <- enrichplot::gseaplot2(
        gsea_obj,
        geneSetID = gene_set_ids,
        title = title,
        subplots = subplots
    )
    
    # Save the plot if requested
    if (save_plot) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Create filename from title
        filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_runsum.pdf"))
        
        ggsave(
            filename = filename,
            plot = p,
            width = width,
            height = height,
            dpi = dpi
        )
        
        message("Plot saved to: ", filename)
    }
    
    # Return the plot object
    return(p)
}
