#' Create and save a GSEA enrichment plot
#' 
#' @param gse_object A GSE object containing enrichment analysis results
#' @param gene_set_ids Vector of gene set IDs to plot (numeric or character)
#' @param title Character string for plot title and filename
#' @param output_dir Directory path for saving plots (default: "../imgs/RunSum/")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 5)
#' @param dpi Plot resolution (default: 300)
#' @param subplots Vector indicating which subplots to include (default: c(1,2,3))
#' @return The generated plot object
#' 
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