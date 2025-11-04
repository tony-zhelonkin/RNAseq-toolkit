#' Get Unique Significant Pathway IDs from a List of GSEA Results
#'
#' Extracts unique pathway IDs that are significant (below `padj_cutoff`) in at least
#' one `gseaResult` object within a provided list.
#' Assumes input is a list where each element is a `gseaResult` object from `clusterProfiler`.
#'
#' @param gsea_results_list A list where each element is a `gseaResult` object.
#' @param padj_cutoff Numeric, adjusted p-value cutoff for significance (default: 0.05).
#' @param verbose Logical, whether to print detailed logging (default: FALSE).
#'
#' @return A character vector of unique significant pathway IDs. Returns an empty
#'         character vector if the input list is empty or no significant pathways are found.
#' @export
#' @importFrom methods is slot
#'
#' @examples
#' # Create dummy gseaResult-like objects (replace with actual results)
#' make_dummy_gsea <- function(ids, padjs, cores) {
#'   res <- data.frame(ID = ids, p.adjust = padjs, core_enrichment = cores)
#'   obj <- new("gseaResult", result = res) # Simplified object
#'   return(obj)
#' }
#' gsea1 <- make_dummy_gsea(c("P1", "P2"), c(0.01, 0.06), c("G1/G2", "G3/G4"))
#' gsea2 <- make_dummy_gsea(c("P1", "P3"), c(0.07, 0.04), c("G1/G2/G5", "G6/G7"))
#' gsea_list <- list(result1 = gsea1, result2 = gsea2)
#'
#' # Get unique significant pathways with padj < 0.05
#' sig_paths <- get_significant_pathways(gsea_list, padj_cutoff = 0.05)
#' print(sig_paths) # Should be c("P1", "P3")
#'
get_significant_pathways <- function(gsea_results_list, padj_cutoff = 0.05, verbose = FALSE) {

    # --- Input Validation ---
    if (!is.list(gsea_results_list)) {
        stop("`gsea_results_list` must be a list.")
    }
    if (!is.numeric(padj_cutoff) || length(padj_cutoff) != 1 || padj_cutoff < 0 || padj_cutoff > 1) {
        stop("`padj_cutoff` must be a single numeric value between 0 and 1.")
    }
    # ------------------------

    if (verbose) {
        message("[get_significant_pathways] Processing ", length(gsea_results_list), " GSEA result(s)")
        message("[get_significant_pathways] Using padj cutoff: ", padj_cutoff)
    }

    # Use lapply to process each element in the list
    significant_path_ids_list <- lapply(names(gsea_results_list), function(result_name) {
        gsea_obj <- gsea_results_list[[result_name]]
        # Validate individual element
        if (methods::is(gsea_obj, "gseaResult") &&
            methods::.hasSlot(gsea_obj, "result") &&
            is.data.frame(gsea_obj@result) &&
            nrow(gsea_obj@result) > 0 &&
            all(c("ID", "p.adjust") %in% colnames(gsea_obj@result))) {

            # Filter results and extract IDs
            results_df <- gsea_obj@result
            n_total <- nrow(results_df)
            sig_ids <- results_df$ID[results_df$p.adjust < padj_cutoff]
            n_sig <- length(sig_ids)
            n_nonsig <- n_total - n_sig

            if (verbose) {
                message(sprintf("[get_significant_pathways]   %s: %d total pathways, %d significant (padj < %.3f), %d filtered out",
                               result_name, n_total, n_sig, padj_cutoff, n_nonsig))
            }

            # Ensure it's not NULL if no rows meet criteria
            if (n_sig == 0) {
                return(character(0))
            } else {
                return(sig_ids)
            }
        } else {
            # Return empty character vector for invalid/empty objects
            if (verbose) {
                message(sprintf("[get_significant_pathways]   %s: Invalid/empty gseaResult object, skipping", result_name))
            } else {
                warning("An element in `gsea_results_list` was not a valid/non-empty gseaResult object. Skipping.")
            }
            return(character(0))
        }
    })

    # Combine all significant IDs and get unique ones
    all_sig_ids <- unlist(significant_path_ids_list)
    n_before_unique <- length(all_sig_ids)
    unique_significant_ids <- unique(all_sig_ids)
    n_after_unique <- length(unique_significant_ids)

    if (verbose) {
        message(sprintf("[get_significant_pathways] Total significant pathways before deduplication: %d", n_before_unique))
        message(sprintf("[get_significant_pathways] Unique significant pathways: %d", n_after_unique))
        if (n_after_unique > 0) {
            message(sprintf("[get_significant_pathways] Removed %d duplicates", n_before_unique - n_after_unique))
        }
    }

    # Return empty vector if NULL (can happen if input list was empty)
    if (is.null(unique_significant_ids)) {
        return(character(0))
    } else {
        return(unique_significant_ids)
    }
}
