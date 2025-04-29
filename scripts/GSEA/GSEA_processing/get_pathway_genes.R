#' Get Core Enrichment Genes for Significant Pathways from a Single GSEA Result
#'
#' Extracts core enrichment genes for significant pathways from a single `gseaResult`
#' object (from `clusterProfiler`), based on an adjusted p-value cutoff.
#' Pathways are ranked by adjusted p-value.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param padj_cutoff Numeric, adjusted p-value cutoff for significance (default: 0.05).
#' @param top Integer or NULL, the number of top pathways to return, ranked by `p.adjust`.
#'        If NULL (default), returns all pathways meeting the `padj_cutoff`.
#'
#' @return A named list where names are pathway IDs (from the `ID` column of the GSEA result)
#'         and values are character vectors of core enrichment genes. Returns an empty list
#'         if the input is invalid or no significant pathways are found.
#' @export
#' @importFrom methods is slot
#' @importFrom utils head
#'
#' @examples
#' # Create a dummy gseaResult-like object (replace with actual results)
#' make_dummy_gsea <- function(ids, padjs, cores) {
#'   res <- data.frame(ID = ids, p.adjust = padjs, core_enrichment = cores)
#'   obj <- new("gseaResult", result = res) # Simplified object
#'   return(obj)
#' }
#' gsea_res <- make_dummy_gsea(c("P1", "P2", "P3"), c(0.01, 0.06, 0.04), c("G1/G2", "G3/G4", "G5/G6"))
#'
#' # Get top 2 pathways with padj < 0.05
#' top_genes <- get_pathway_genes(gsea_res, padj_cutoff = 0.05, top = 2)
#' print(top_genes) # Should be P1 and P3 genes
#'
get_pathway_genes <- function(gsea_obj, padj_cutoff = 0.05, top = NULL) {

    # --- Input Validation ---
    if (!methods::is(gsea_obj, "gseaResult")) {
         warning("Input `gsea_obj` is not a gseaResult object. Returning empty list.")
         return(list())
    }
     # Check if the result slot is accessible and is a data frame with rows
    if (!methods::.hasSlot(gsea_obj, "result") || !is.data.frame(gsea_obj@result) || nrow(gsea_obj@result) == 0) {
        warning("Input `gsea_obj` has an invalid or empty result slot. Returning empty list.")
        return(list())
    }
     if (!is.numeric(padj_cutoff) || length(padj_cutoff) != 1 || padj_cutoff < 0 || padj_cutoff > 1) {
        stop("`padj_cutoff` must be a single numeric value between 0 and 1.")
    }
    if (!is.null(top) && (!is.numeric(top) || length(top) != 1 || top <= 0 || top != round(top))) {
        stop("`top` must be NULL or a single positive integer.")
    }
    # Check required columns exist
    required_cols <- c("ID", "p.adjust", "core_enrichment")
    if (!all(required_cols %in% colnames(gsea_obj@result))) {
        stop(sprintf("Input `gsea_obj@result` is missing required columns: %s",
                     paste(setdiff(required_cols, colnames(gsea_obj@result)), collapse=", ")))
    }
    # ------------------------

    results_df <- gsea_obj@result

    # Filter for significant pathways
    sig_results <- results_df[results_df$p.adjust < padj_cutoff, , drop = FALSE]

    # Check if any significant results
    if (nrow(sig_results) == 0) {
        # warning(sprintf("No significant pathways found with p.adjust < %f.", padj_cutoff)) # Optional warning
        return(list())
    }

    # Sort by adjusted p-value
    sig_results <- sig_results[order(sig_results$p.adjust), , drop = FALSE]

    # Take top N pathways if specified
    if (!is.null(top)) {
        n_to_keep <- min(top, nrow(sig_results))
        sig_results <- utils::head(sig_results, n_to_keep)
    }

    # Create list of gene sets, handling potential NA/empty core_enrichment
    pathway_genes_list <- lapply(sig_results$core_enrichment, function(core_str) {
        if (is.character(core_str) && !is.na(core_str) && nchar(core_str) > 0) {
            # Split the string by '/'
            genes <- strsplit(core_str, "/")[[1]]
            # Remove any empty strings resulting from splitting (e.g., trailing '/')
            genes <- genes[nzchar(genes)]
            return(genes)
        } else {
            # Return empty character vector if core_enrichment is NA, empty, or not character
            return(character(0))
        }
    })

    # Set names of the list using the pathway IDs
    names(pathway_genes_list) <- sig_results$ID

    return(pathway_genes_list)
}
