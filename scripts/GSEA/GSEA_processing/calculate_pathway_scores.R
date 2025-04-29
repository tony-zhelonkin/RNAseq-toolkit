#' Calculate Pathway Scores from Expression Data
#'
#' This function calculates pathway scores by averaging the expression values of genes
#' in each pathway across all samples.
#'
#' @param expression_data A matrix or data frame of normalized expression data,
#'        with genes as rows and samples as columns.
#' @param pathway_genes A named list where names are pathway IDs and values are
#'        vectors of gene identifiers.
#' @param method Character, method to calculate pathway scores (default: "mean").
#'        Currently only "mean" is implemented.
#'
#' @return A matrix of pathway scores, with samples as rows and pathways as columns.
#' @export
#'
#' @examples
#' # Assuming expression_data is a normalized expression matrix and
#' # pathway_genes is a list of genes for each pathway
#' pathway_scores <- calculate_pathway_scores(expression_data, pathway_genes)
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
