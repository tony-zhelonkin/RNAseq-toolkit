#' Get Core Enrichment Genes for Significant Pathways
#'
#' This function extracts core enrichment genes for significant pathways from a GSEA result object,
#' based on adjusted p-value cutoff.
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
  # Extract results data frame
  results_df <- gsea_obj@result
  
  # Filter for significant pathways
  sig_results <- results_df[results_df$p.adjust < q_cutoff, ]
  
  # Sort by adjusted p-value if needed
  if (!is.null(top)) {
    sig_results <- sig_results[order(sig_results$p.adjust), ]
    sig_results <- sig_results[1:min(top, nrow(sig_results)), ]
  }
  
  # Extract core enrichment genes
  pathway_genes <- strsplit(sig_results$core_enrichment, "/")
  names(pathway_genes) <- sig_results$ID
  
  return(pathway_genes)
}
