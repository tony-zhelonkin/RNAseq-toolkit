#' Get Core Enrichment Genes for Top Significant Pathways Across Contrasts
#'
#' Extracts core enrichment genes for pathways found significant (below `padj_cutoff`)
#' in at least one contrast for a given database. Ranks pathways by their minimum
#' adjusted p-value across all contrasts and returns the gene lists for the top N.
#' Assumes input list contains `gseaResult` objects from `clusterProfiler`.
#'
#' @param gsea_results_list A nested list: `list(ContrastName = list(DatabaseName = gseaResult))`.
#' @param database Character, the name of the database (e.g., "HALLMARK", "KEGG") to extract from.
#'        Must match a key in the inner lists of `gsea_results_list`.
#' @param padj_cutoff Numeric, the adjusted p-value threshold for significance (default: 0.05).
#' @param top Integer or NULL, the number of top pathways to return, ranked by minimum
#'        `p.adjust` across contrasts. If NULL (default), returns all significant pathways found.
#'
#' @return A named list where names are pathway IDs (Description or ID from GSEA results)
#'         and values are character vectors of core enrichment genes. Returns an empty list
#'         if no significant pathways are found for the specified database.
#' @export
#' @importFrom methods is slot
#' @importFrom stats setNames
#'
#' @examples
#' # Create dummy gseaResult-like objects (replace with actual results)
#' make_dummy_gsea <- function(ids, padjs, cores) {
#'   res <- data.frame(ID = ids, p.adjust = padjs, core_enrichment = cores)
#'   obj <- new("gseaResult", result = res) # Simplified object
#'   return(obj)
#' }
#' gsea1_h <- make_dummy_gsea(c("P1", "P2"), c(0.01, 0.06), c("G1/G2", "G3/G4"))
#' gsea2_h <- make_dummy_gsea(c("P1", "P3"), c(0.04, 0.005), c("G1/G2/G5", "G6/G7"))
#' gsea_list <- list(Contrast1 = list(HALLMARK = gsea1_h), Contrast2 = list(HALLMARK = gsea2_h))
#'
#' # Get top 1 pathway genes based on min p.adjust < 0.05
#' top_genes <- get_pathway_genes_all(gsea_list, "HALLMARK", padj_cutoff = 0.05, top = 1)
#' print(top_genes) # Should be P3 genes
#'
#' # Get all significant pathway genes
#' all_sig_genes <- get_pathway_genes_all(gsea_list, "HALLMARK", padj_cutoff = 0.05)
#' print(all_sig_genes) # Should be P1 and P3 genes
#'
get_pathway_genes_all <- function(gsea_results_list, database, padj_cutoff = 0.05, top = NULL) {

    # --- Input Validation ---
    if (!is.list(gsea_results_list) || length(gsea_results_list) == 0) {
        warning("`gsea_results_list` is empty or not a list. Returning empty list.")
        return(list())
    }
    if (!is.character(database) || length(database) != 1) {
        stop("`database` must be a single character string.")
    }
    if (!is.numeric(padj_cutoff) || length(padj_cutoff) != 1 || padj_cutoff < 0 || padj_cutoff > 1) {
        stop("`padj_cutoff` must be a single numeric value between 0 and 1.")
    }
    if (!is.null(top) && (!is.numeric(top) || length(top) != 1 || top <= 0 || top != round(top))) {
        stop("`top` must be NULL or a single positive integer.")
    }
    # ------------------------

    # Store results: pathway_id -> list(min_padj = Inf, genes = character(0))
    pathway_data <- list()

    # Single pass through the results list
    for (contrast_name in names(gsea_results_list)) {
        contrast_results <- gsea_results_list[[contrast_name]]

        if (database %in% names(contrast_results)) {
            gsea_obj <- contrast_results[[database]]

            # Check if it's a valid GSEA result object (basic check)
            if (methods::is(gsea_obj, "gseaResult") && nrow(gsea_obj@result) > 0) {
                results_df <- gsea_obj@result

                # Filter for significant pathways in this contrast
                sig_results_contrast <- results_df[results_df$p.adjust < padj_cutoff, , drop = FALSE]

                if (nrow(sig_results_contrast) > 0) {
                    for (i in 1:nrow(sig_results_contrast)) {
                        pathway_id <- sig_results_contrast$ID[i] # Or Description? Use ID for consistency
                        current_padj <- sig_results_contrast$p.adjust[i]
                        core_enrichment_str <- sig_results_contrast$core_enrichment[i]

                        # Extract genes
                        current_genes <- character(0)
                        if (is.character(core_enrichment_str) && !is.na(core_enrichment_str) && nchar(core_enrichment_str) > 0) {
                            current_genes <- strsplit(core_enrichment_str, "/")[[1]]
                        }

                        # Update pathway data
                        if (pathway_id %in% names(pathway_data)) {
                            # Update p-value if current is lower
                            pathway_data[[pathway_id]]$min_padj <- min(pathway_data[[pathway_id]]$min_padj, current_padj)
                            # Optionally merge genes? For now, keep genes from the first encounter or the one with min p-value?
                            # Let's keep the genes from the result with the minimum p-value encountered so far.
                            if (current_padj < pathway_data[[pathway_id]]$min_padj) {
                                 pathway_data[[pathway_id]]$genes <- current_genes
                            }
                        } else {
                            # Add new pathway
                            pathway_data[[pathway_id]] <- list(min_padj = current_padj, genes = current_genes)
                        }
                    }
                }
            } else {
                 warning(sprintf("Item '%s' in contrast '%s' is not a valid or non-empty gseaResult object. Skipping.",
                                 database, contrast_name))
            }
        } else {
             # warning(sprintf("Database '%s' not found in contrast '%s'. Skipping.", database, contrast_name))
             # This might be too verbose if many contrasts lack the database
        }
    }

    # If no significant pathways found across all contrasts
    if (length(pathway_data) == 0) {
        warning(sprintf("No pathways found for database '%s' with p.adjust < %f.", database, padj_cutoff))
        return(list())
    }

    # Extract minimum p-values and sort
    min_pvals <- sapply(pathway_data, function(x) x$min_padj)
    ordered_pathway_ids <- names(sort(min_pvals))

    # Select top N if specified
    if (!is.null(top)) {
        n_to_keep <- min(top, length(ordered_pathway_ids))
        selected_pathway_ids <- ordered_pathway_ids[1:n_to_keep]
    } else {
        selected_pathway_ids <- ordered_pathway_ids
    }

    # Extract the gene lists for the selected pathways
    final_pathway_genes <- lapply(pathway_data[selected_pathway_ids], function(x) x$genes)
    # Ensure names are preserved
    names(final_pathway_genes) <- selected_pathway_ids

    return(final_pathway_genes)
}
