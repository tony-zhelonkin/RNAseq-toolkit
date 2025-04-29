#' Get Significant Pathways from GSEA Results
#'
#' This function extracts significant pathways from GSEA results based on a q-value cutoff.
#' It includes robust error handling for various edge cases.
#'
#' @param gsea_results A list of GSEA result objects.
#' @param q_cutoff Numeric, q-value (adjusted p-value) cutoff for significance (default: 0.01).
#'
#' @return A character vector of unique pathway IDs that pass the significance threshold.
#' @export
#'
#' @examples
#' # Assuming gsea_results is a list of GSEA result objects
#' significant_pathways <- get_significant_pathways(gsea_results, q_cutoff = 0.05)
get_significant_pathways <- function(gsea_results, q_cutoff = 0.01) {
    significant_paths <- lapply(gsea_results, function(x) {
        if(!is.null(x) && length(x) > 0 && !is.character(x) && nrow(x@result) > 0) {
            x$ID[x$p.adjust < q_cutoff]
        } else {
            character(0)  # Return empty character vector if no results
        }
    })
    unique(unlist(significant_paths))
}
