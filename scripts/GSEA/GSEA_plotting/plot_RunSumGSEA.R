#' Create and Save a GSEA Running Sum Enrichment Plot
#' 
#' This function creates and saves a GSEA enrichment plot showing the running sum
#' statistic for specified gene sets, using the enrichplot package.
#' 
#' @param gse_object A GSEA result object from clusterProfiler.
#' @param gene_set_ids Vector of gene set IDs to plot (numeric indices or character IDs).
#' @param title Character, plot title and base filename for saving (required).
#' @param output_dir Character, directory path for saving plots (default: "3_Results/imgs/RunSum/").
#' @param width Numeric, plot width in inches (default: 10).
#' @param height Numeric, plot height in inches (default: 5).
#' @param dpi Numeric, plot resolution in dots per inch (default: 300).
#' @param subplots Numeric vector, indicating which subplots to include (default: c(1,2,3)):
#'        1 = running enrichment score, 2 = positions of gene set members, 3 = ranking metric scores.
#'
#' @return The generated plot object (invisibly).
#' @export
#'
#' @examples
#' # Assuming gse_object is a GSEA result object from clusterProfiler
#' runSumGSEAplot(gse_object, gene_set_ids = 1:5, 
#'               title = "HALLMARK running sum for top 5 pathways")
runSumGSEAplot <- function(gse_object, 
                            gene_set_ids,
                            title,
                            output_dir = "3_Results/imgs/RunSum/",
                            width = 10,
                            height = 5,
                            dpi = 300,
                            subplots = c(1, 2, 3)) {
    
    # Input validation
    if (missing(title) || is.null(title) || title == "") {
        stop("Title must be provided")
    }
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Create the GSEA plot
    p <- enrichplot::gseaplot2(
        gse_object,
        geneSetID = gene_set_ids,
        title = title,
        subplots = subplots
    )
    
    # Create filename from title (replace spaces with underscores)
    filename <- file.path(output_dir, paste0(gsub(" ", "_", title), ".pdf"))
    
    # Save the plot
    ggsave(
        filename = filename,
        plot = p,
        width = width,
        height = height,
        dpi = dpi
    )
    
    # Print confirmation message
    message("Plot saved to: ", filename)
    
    # Return the plot object
    return(p)
}
