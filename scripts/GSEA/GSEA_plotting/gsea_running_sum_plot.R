#' Create GSEA Running Sum Enrichment Plot(s)
#'
#' Generates GSEA enrichment plot(s) showing the running sum statistic for one
#' or more specified gene sets, using `enrichplot::gseaplot2`.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param gene_set_ids Character or numeric vector. Specifies which gene set(s) to plot,
#'        either by their ID (character) or their index (numeric) in the `gsea_obj` result table.
#' @param title Character, the main title for the plot(s). If multiple `gene_set_ids` are
#'        provided, this title is often reused or adapted by `gseaplot2`.
#' @param subplots Numeric vector, indicating which subplots to include (default: c(1, 2, 3)):
#'        1 = running enrichment score, 2 = positions of gene set members along the ranked list,
#'        3 = ranking metric scores for all genes.
#'
#' @return A ggplot object, or potentially a list of ggplot objects if multiple `gene_set_ids`
#'         are plotted individually by `gseaplot2` (behavior depends on `enrichplot` version).
#'         The plot is NOT saved automatically.
#' @export
#' @importFrom enrichplot gseaplot2
#' @importFrom methods is slot
#'
#' @examples
#' # Assuming gsea_res is a valid gseaResult object
#'
#' # Plot top pathway (index 1)
#' # p1 <- gsea_running_sum_plot(gsea_res, gene_set_ids = 1, title = gsea_res$Description[1])
#' # print(p1)
#'
#' # Plot pathways with specific IDs
#' # pathway_ids_to_plot <- c("HALLMARK_APOPTOSIS", "HALLMARK_DNA_REPAIR")
#' # p_list <- gsea_running_sum_plot(gsea_res, gene_set_ids = pathway_ids_to_plot,
#' #                                title = "Selected Hallmark Pathways")
#' # print(p_list[[1]]) # Print the first plot if multiple are returned
#' # To save: ggsave("apoptosis_running_sum.png", p_list[[1]]) # Save individual plots if needed

gsea_running_sum_plot <- function(gsea_obj,
                                 gene_set_ids,
                                 title = NULL, # Made title optional, gseaplot2 uses pathway name
                                 subplots = c(1, 2, 3)) {

    # --- Input Validation ---
    if (!methods::is(gsea_obj, "gseaResult")) {
         stop("Input `gsea_obj` must be a gseaResult object from clusterProfiler.")
    }
     if (!methods::.hasSlot(gsea_obj, "result") || !is.data.frame(gsea_obj@result) || nrow(gsea_obj@result) == 0) {
        stop("Input `gsea_obj` has an invalid or empty result slot.")
    }
    if (missing(gene_set_ids) || length(gene_set_ids) == 0) {
        stop("`gene_set_ids` must be provided (numeric index or character ID).")
    }
    # Check if IDs/indices are valid
    if (is.numeric(gene_set_ids)) {
        if (any(gene_set_ids > nrow(gsea_obj@result)) || any(gene_set_ids <= 0)) {
            stop("Numeric `gene_set_ids` are out of bounds.")
        }
    } else if (is.character(gene_set_ids)) {
        if (!all(gene_set_ids %in% gsea_obj@result$ID)) {
             missing_ids <- setdiff(gene_set_ids, gsea_obj@result$ID)
             warning("Some character `gene_set_ids` not found in gsea_obj: ", paste(missing_ids, collapse=", "))
             gene_set_ids <- intersect(gene_set_ids, gsea_obj@result$ID) # Keep only valid ones
             if (length(gene_set_ids) == 0) stop("No valid character `gene_set_ids` found.")
        }
    } else {
        stop("`gene_set_ids` must be a numeric or character vector.")
    }
    # ------------------------


    # Create the GSEA plot using enrichplot::gseaplot2
    # Note: gseaplot2 might return a single plot or a list depending on input/version
    p <- enrichplot::gseaplot2(
        x = gsea_obj,
        geneSetID = gene_set_ids,
        title = title, # Pass title, gseaplot2 might override/use pathway name
        subplots = subplots
        # Add other gseaplot2 parameters if needed (e.g., pvalue_table = TRUE)
    )

    # Return the plot object(s)
    return(p)
}
