#' Calculate Pathway Scores from Expression Data
#'
#' Calculates pathway scores by averaging the expression of genes within each pathway
#' across all samples. Handles cases where pathways have no genes found in the
#' expression data.
#'
#' @param expression_data A numeric matrix or data frame of normalized expression data
#'        (genes as rows, samples as columns). Rownames must be gene identifiers,
#'        and colnames must be sample identifiers.
#' @param pathway_genes A named list where each element is a character vector of
#'        gene identifiers belonging to a pathway. The names of the list elements
#'        are the pathway IDs.
#' @param method Character, the method for calculating scores. Currently, only
#'        "mean" (default) is implemented.
#' @param verbose Logical, if TRUE, print status messages (default: FALSE).
#'
#' @return A numeric matrix of pathway scores, with samples as rows and pathways
#'         as columns. Returns an empty matrix with correct sample names if
#'         `pathway_genes` is empty or NULL.
#' @export
#'
#' @examples
#' # Create dummy data
#' expr_mat <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10,
#'                    dimnames = list(paste0("Gene", 1:100), paste0("Sample", 1:10)))
#' path_list <- list(PathwayA = paste0("Gene", 1:10),
#'                   PathwayB = paste0("Gene", 11:25),
#'                   PathwayC = paste0("Gene", 95:105)) # PathwayC has some missing genes
#'
#' # Calculate scores
#' pathway_scores <- calculate_pathway_scores(expr_mat, path_list, verbose = TRUE)
#' print(head(pathway_scores))
#'
calculate_pathway_scores <- function(expression_data, pathway_genes, method = "mean", verbose = FALSE) {

    # --- Input Validation ---
    if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
        stop("`expression_data` must be a matrix or data frame.")
    }
    if (is.null(rownames(expression_data))) {
        stop("`expression_data` must have rownames (gene identifiers).")
    }
    if (is.null(colnames(expression_data))) {
        stop("`expression_data` must have colnames (sample identifiers).")
    }
    if (!is.list(pathway_genes) && !is.null(pathway_genes)) {
        stop("`pathway_genes` must be a named list or NULL.")
    }
     if (!is.null(pathway_genes) && (is.null(names(pathway_genes)) || any(names(pathway_genes) == ""))) {
        stop("`pathway_genes` list must be named (pathway IDs).")
    }
    if (method != "mean") {
        stop("Invalid `method`. Currently only 'mean' is supported.")
    }
    # ------------------------

    if (verbose) {
        message(sprintf("Input: %d genes, %d samples.", nrow(expression_data), ncol(expression_data)))
        message(sprintf("Processing %d pathways using '%s' method.", length(pathway_genes), method))
    }

    # Handle empty or NULL pathway list
    if (is.null(pathway_genes) || length(pathway_genes) == 0) {
        warning("Empty or NULL `pathway_genes` provided. Returning an empty score matrix.")
        return(matrix(NA_real_,
                     nrow = ncol(expression_data),
                     ncol = 0,
                     dimnames = list(colnames(expression_data), character(0))))
    }

    # Initialize scores matrix (samples x pathways)
    sample_names <- colnames(expression_data)
    pathway_names <- names(pathway_genes)
    scores <- matrix(NA_real_,
                    nrow = length(sample_names),
                    ncol = length(pathway_names),
                    dimnames = list(sample_names, pathway_names))

    # Get gene universe from expression data once
    expression_genes <- rownames(expression_data)

    # Loop through each pathway
    for (pathway_id in pathway_names) {
        pathway_gene_set <- pathway_genes[[pathway_id]]

        # Find intersection of pathway genes and expression data genes
        genes_in_pathway_and_expr <- intersect(pathway_gene_set, expression_genes)

        if (length(genes_in_pathway_and_expr) > 0) {
            # Extract expression for these genes across all samples
            # Ensure it remains a matrix even if only one gene is found
            pathway_exp_subset <- expression_data[genes_in_pathway_and_expr, , drop = FALSE]

            # Calculate score for each sample based on the method
            if (method == "mean") {
                # Use colMeans for efficiency and NA handling
                pathway_sample_scores <- colMeans(pathway_exp_subset, na.rm = TRUE)
                scores[, pathway_id] <- pathway_sample_scores
            }
            # Add other methods here with 'else if (method == "...")' in the future

            if (verbose && any(is.na(pathway_sample_scores))) {
                 message(sprintf("Pathway '%s': Some samples have NA scores (possibly due to all genes having NA).", pathway_id))
            }

        } else {
            # No genes for this pathway found in the expression data
            if (verbose) {
                warning(sprintf("Pathway '%s': No overlapping genes found in expression data. Assigning NA scores.", pathway_id))
            }
            # Scores are already NA by initialization, but can be explicit:
            # scores[, pathway_id] <- NA_real_
        }
    }

    if (verbose) {
        message(sprintf("Finished calculation. Returning scores matrix: %d samples x %d pathways.",
                       nrow(scores), ncol(scores)))
    }

    return(scores)
}
