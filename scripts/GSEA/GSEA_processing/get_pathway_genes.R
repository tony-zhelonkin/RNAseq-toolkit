# Define function for single GSEA result
get_pathway_genes <- function(gsea_results, top_n = 30) {
    # Extract results data frame
    results_df <- gsea_results@result
    
    # Get significant pathways and sort by adjusted p-value
    sig_results <- results_df[results_df$p.adjust < 0.01, ]
    sig_results <- sig_results[order(sig_results$p.adjust), ]
    
    # Take top N pathways
    sig_results <- head(sig_results, top_n)
    
    # Create list of gene sets
    pathway_genes <- strsplit(sig_results$core_enrichment, "/")
    names(pathway_genes) <- sig_results$ID
    
    return(pathway_genes)
}