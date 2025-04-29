#' Get Core Enrichment Genes for Significant Pathways
#'
#' This function extracts core enrichment genes for significant pathways from a GSEA result object,
#' based on adjusted p-value cutoff. It includes robust error handling for various edge cases.
#'
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param q_cutoff Numeric, q-value (adjusted p-value) cutoff for significance (default: 0.01).
#' @param top Integer, number of top pathways to return, ranked by adjusted p-value (default: NULL, returns all).
#'
#' @return A named list where names are pathway IDs and values are vectors of core enrichment genes.
#' @export
#'
#' @examples
#' # Assuming gsea_obj is a GSEA result object from clusterProfiler
#' pathway_genes <- get_pathway_genes(gsea_obj, q_cutoff = 0.05, top = 10)
get_pathway_genes <- function(gsea_obj, q_cutoff = 0.01, top = NULL) {
    # Check if valid GSEA result
    if (is.null(gsea_obj) || length(gsea_obj) == 0 || is.character(gsea_obj) || nrow(gsea_obj@result) == 0) {
        return(list())  # Return empty list if no results
    }
    
    # Extract results data frame
    results_df <- gsea_obj@result
    
    # Filter for significant pathways
    sig_results <- results_df[results_df$p.adjust < q_cutoff, ]
    
    # Check if any significant results
    if (nrow(sig_results) == 0) {
        return(list())  # Return empty list if no significant results
    }
    
    # Sort by adjusted p-value
    sig_results <- sig_results[order(sig_results$p.adjust), ]
    
    # Take top N pathways if specified
    if (!is.null(top) && top < nrow(sig_results)) {
        sig_results <- head(sig_results, top)
    }
    
    # Create list of gene sets
    pathway_genes <- lapply(sig_results$core_enrichment, function(x) {
        if (is.character(x) && !is.na(x) && nchar(x) > 0) {
            strsplit(x, "/")[[1]]
        } else {
            character(0)  # Return empty character vector if no genes
        }
    })
    names(pathway_genes) <- sig_results$ID
    
    return(pathway_genes)
}
