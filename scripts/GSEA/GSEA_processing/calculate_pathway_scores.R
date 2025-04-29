#' Calculate Pathway Scores from Expression Data
#'
#' This function calculates pathway scores by averaging the expression values of genes
#' in each pathway across all samples. It includes robust error handling for various edge cases.
#'
#' @param expression_data A matrix or data frame of normalized expression data,
#'        with genes as rows and samples as columns.
#' @param pathway_genes A named list where names are pathway IDs and values are
#'        vectors of gene identifiers.
#' @param method Character, method to calculate pathway scores (default: "mean").
#'        Currently only "mean" is implemented.
#' @param verbose Logical, whether to print debug messages (default: FALSE).
#'
#' @return A matrix of pathway scores, with samples as rows and pathways as columns.
#' @export
#'
#' @examples
#' # Assuming expression_data is a normalized expression matrix and
#' # pathway_genes is a list of genes for each pathway
#' pathway_scores <- calculate_pathway_scores(expression_data, pathway_genes)
calculate_pathway_scores <- function(expression_data, pathway_genes, method = "mean", verbose = FALSE) {
    # Debug information
    if (verbose) {
        message(sprintf("calculate_pathway_scores: expression_data has %d genes and %d samples", 
                       nrow(expression_data), ncol(expression_data)))
        message(sprintf("calculate_pathway_scores: processing %d pathways", length(pathway_genes)))
    }
    
    # If no pathways, return empty matrix with all samples
    if (length(pathway_genes) == 0) {
        warning("No pathways provided to calculate_pathway_scores, returning empty matrix")
        return(matrix(NA, 
                     nrow = ncol(expression_data), 
                     ncol = 0,
                     dimnames = list(colnames(expression_data), character(0))))
    }
    
    # Initialize scores matrix with all samples
    scores <- matrix(NA, 
                    nrow = ncol(expression_data),
                    ncol = length(pathway_genes),
                    dimnames = list(colnames(expression_data), 
                                  names(pathway_genes)))
    
    for (pathway in names(pathway_genes)) {
        genes <- pathway_genes[[pathway]]
        # Get genes that are present in expression data
        genes <- genes[genes %in% rownames(expression_data)]
        
        if (length(genes) > 0) {
            # Extract expression for these genes across ALL samples
            pathway_exp <- expression_data[genes, , drop = FALSE]
            
            # Calculate mean expression for each sample
            scores[, pathway] <- colMeans(pathway_exp, na.rm = TRUE)
        } else {
            if (verbose) {
                warning(sprintf("No genes found for pathway '%s'", pathway))
            }
            # Set to NA for this pathway
            scores[, pathway] <- NA
        }
    }
    
    # Debug information
    if (verbose) {
        message(sprintf("calculate_pathway_scores: returning scores matrix with %d samples and %d pathways", 
                       nrow(scores), ncol(scores)))
    }
    
    return(scores)
}
