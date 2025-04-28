calculate_pathway_scores <- function(expression_data, pathway_genes, method = "mean") {
    scores <- matrix(NA, 
                    nrow = ncol(expression_data),
                    ncol = length(pathway_genes),
                    dimnames = list(colnames(expression_data), 
                                  names(pathway_genes)))
    
    for(pathway in names(pathway_genes)) {
        genes <- pathway_genes[[pathway]]
        # Get genes that are present in expression data
        genes <- genes[genes %in% rownames(expression_data)]
        pathway_exp <- expression_data[genes, ]
        scores[, pathway] <- colMeans(pathway_exp, na.rm = TRUE)
    }
    return(scores)
}