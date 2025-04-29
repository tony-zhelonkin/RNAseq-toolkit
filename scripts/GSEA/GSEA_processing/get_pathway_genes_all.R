#' Get Core Enrichment Genes for Significant Pathways Across All Contrasts
#'
#' This function extracts core enrichment genes for significant pathways from GSEA results
#' across all contrasts, and ranks pathways by their minimum adjusted p-value.
#'
#' @param gsea_results_list A nested list of GSEA result objects, organized by contrast and database.
#' @param database Character, name of the database to extract pathways from (e.g., "hallmark", "kegg").
#' @param top Integer, number of top pathways to return, ranked by minimum adjusted p-value (default: NULL, returns all).
#'
#' @return A named list where names are pathway IDs and values are vectors of core enrichment genes.
#' @export
#'
#' @examples
#' # Assuming gsea_results_list is a nested list of GSEA results by contrast and database
#' hallmark_genes <- get_pathway_genes_all(gsea_results_list, "hallmark", top = 20)
get_pathway_genes_all <- function(gsea_results_list, database, top = NULL) {
    # Combine results from all contrasts
    all_pathway_genes <- list()
    
    for(contrast in names(gsea_results_list)) {
        if(!is.null(gsea_results_list[[contrast]][[database]])) {
            results_df <- gsea_results_list[[contrast]][[database]]@result
            sig_results <- results_df[results_df$p.adjust < 0.01, ]
            pathway_genes <- strsplit(sig_results$core_enrichment, "/")
            names(pathway_genes) <- sig_results$ID
            all_pathway_genes <- c(all_pathway_genes, pathway_genes)
        }
    }
    
    # Remove duplicates
    unique_pathways <- unique(names(all_pathway_genes))
    unique_pathway_genes <- all_pathway_genes[unique_pathways]
    
    # Get minimum p-value for each unique pathway across all contrasts
    min_pvals <- sapply(unique_pathways, function(pathway) {
        min_pval <- Inf
        for(contrast in names(gsea_results_list)) {
            if(!is.null(gsea_results_list[[contrast]][[database]])) {
                results_df <- gsea_results_list[[contrast]][[database]]@result
                if(pathway %in% results_df$ID) {
                    min_pval <- min(min_pval, results_df$p.adjust[results_df$ID == pathway])
                }
            }
        }
        return(min_pval)
    })
    
    # Sort by minimum p-value and take top N if specified
    if(!is.null(top)) {
        ordered_pathways <- names(sort(min_pvals))[1:min(top, length(min_pvals))]
        unique_pathway_genes <- unique_pathway_genes[ordered_pathways]
    }
    
    return(unique_pathway_genes)
}
