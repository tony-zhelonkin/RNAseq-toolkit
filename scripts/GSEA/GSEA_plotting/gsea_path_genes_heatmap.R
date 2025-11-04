#' Create Heatmap for Core Enrichment Genes of a GSEA Pathway
#'
#' Generates a heatmap using `pheatmap` showing the expression patterns of core
#' enrichment genes for a specified pathway identified from a `gseaResult` object.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param pathway_name Character, the exact pathway ID or Description to visualize.
#'        The function first searches for an exact match in the `ID` column, then
#'        in the `Description` column of the `gsea_obj@result`.
#' @param expression_data A numeric matrix or data frame of normalized expression values
#'        (e.g., log-CPM). Rownames must be gene identifiers matching those in GSEA results,
#'        and colnames must be sample identifiers.
#' @param sample_order Character vector specifying the exact order (and subset) of samples
#'        (colnames in `expression_data`) to include in the heatmap columns.
#' @param annotation_col Data frame or NULL. Sample annotations for column annotation bars.
#'        Rownames must match the sample names specified in `sample_order`. Default: NULL.
#' @param annotation_colors List or NULL. A list specifying colors for annotation tracks,
#'        matching the structure required by `pheatmap`. Default: NULL.
#' @param title Character or NULL. Optional title for the heatmap. If NULL (default),
#'        a title is generated from the pathway Description (cleaned and wrapped).
#' @param scale Character, scaling parameter passed to `pheatmap`: "row" (default),
#'        "column", or "none".
#' @param cluster_rows Logical, whether to cluster rows (genes) (default: TRUE).
#' @param show_rownames Logical, whether to show gene names as rownames (default: TRUE).
#' @param gaps_col Numeric vector or NULL. Column indices where to draw vertical gaps
#'        in the heatmap (default: NULL). Indices refer to the order in `sample_order`.
#' @param base_fontsize Numeric, base font size for the heatmap (default: 9).
#' @param ... Additional arguments passed directly to `pheatmap::pheatmap`.
#'
#' @return The `pheatmap` object, returned invisibly (useful for saving with `pheatmap::save_pheatmap_pdf` etc.).
#' @export
#' @import pheatmap
#' @importFrom dplyr %>% filter slice
#' @importFrom stringr str_replace_all str_wrap
#' @importFrom methods is slot
#' @importFrom grDevices colorRampPalette
#'
#' @examples
#' # Assuming gsea_res is a valid gseaResult object,
#' # expr_mat is a normalized expression matrix,
#' # sample_ordering is a character vector of sample names,
#' # and sample_annot is an annotation data frame.
#'
#' # Basic usage with pathway ID:
#' # p_obj <- gsea_pathway_heatmap(gsea_res, "GO:0006955", expr_mat, sample_ordering)
#'
#' # With annotations and custom title:
#' # p_obj_annot <- gsea_pathway_heatmap(gsea_res, "HALLMARK_INFLAMMATORY_RESPONSE", expr_mat,
#' #                                     sample_ordering, annotation_col = sample_annot,
#' #                                     title = "Inflammatory Response Genes")
#' # To save: pheatmap::save_pheatmap_pdf(p_obj_annot, "inflammation_heatmap.pdf")

library(pheatmap) # Ensure pheatmap is loaded for colorRampPalette if not explicitly imported

gsea_pathway_heatmap <- function(gsea_obj,
                                pathway_name,
                                expression_data,
                                sample_order,
                                annotation_col = NULL,
                                annotation_colors = NULL,
                                title = NULL,
                                scale = "row", # Renamed scale_expr
                                cluster_rows = TRUE, # Added
                                show_rownames = TRUE, # Added
                                gaps_col = NULL,
                                base_fontsize = 9, # Added
                                ...) {

    # --- Input Validation ---
    if (!methods::is(gsea_obj, "gseaResult")) stop("`gsea_obj` must be a gseaResult object.")
    if (!methods::.hasSlot(gsea_obj, "result") || !is.data.frame(gsea_obj@result) || nrow(gsea_obj@result) == 0) stop("`gsea_obj` has an invalid or empty result slot.")
    required_cols <- c("ID", "Description", "core_enrichment")
    if (!all(required_cols %in% colnames(gsea_obj@result))) stop(sprintf("`gsea_obj@result` missing required columns: %s", paste(setdiff(required_cols, colnames(gsea_obj@result)), collapse=", ")))
    if (!is.matrix(expression_data) && !is.data.frame(expression_data)) stop("`expression_data` must be a matrix or data frame.")
    if (is.null(rownames(expression_data))) stop("`expression_data` must have rownames (gene identifiers).")
    if (!all(sample_order %in% colnames(expression_data))) stop("Not all `sample_order` names found in `expression_data` colnames.")
    if (!is.null(annotation_col) && !all(sample_order %in% rownames(annotation_col))) stop("Rownames of `annotation_col` must match `sample_order`.")
    # ------------------------

    gsea_df <- as.data.frame(gsea_obj@result)

    # Find matching pathway (exact ID first, then exact Description)
    pathway_row <- gsea_df %>% dplyr::filter(.data$ID == pathway_name)
    if (nrow(pathway_row) == 0) {
        pathway_row <- gsea_df %>% dplyr::filter(.data$Description == pathway_name)
    }

    if (nrow(pathway_row) == 0) {
        warning("No pathway found matching ID or Description: ", pathway_name)
        return(invisible(NULL))
    }

    if (nrow(pathway_row) > 1) {
        warning("Multiple pathways matched '", pathway_name, "'. Using the first match (ID: ", pathway_row$ID[1], ").")
        pathway_info <- pathway_row %>% dplyr::slice(1)
    } else {
        pathway_info <- pathway_row
    }

    # Extract core genes
    core_genes <- unlist(strsplit(pathway_info$core_enrichment, "/"))
    core_genes <- core_genes[nzchar(core_genes)] # Remove empty strings if any

    # Find which core genes are present in the expression data
    genes_present <- intersect(core_genes, rownames(expression_data))

    if (length(genes_present) == 0) {
        warning("No core enrichment genes for pathway '", pathway_name, "' found in the provided expression data.")
        return(invisible(NULL))
    }
     if (length(genes_present) < 2 && cluster_rows) {
        warning("Only one gene found for pathway '", pathway_name, "'. Setting cluster_rows = FALSE.")
        cluster_rows <- FALSE
    }


    # Subset expression matrix using the provided sample order
    # Ensure it remains a matrix even with one gene
    expr_sub <- expression_data[genes_present, sample_order, drop = FALSE]

    # Clean up the pathway name for the default title
    cleaned_pathway_name <- pathway_info$Description %>%
        stringr::str_replace_all(c("^GOBP_" = "", "^KEGG_" = "", "^REACTOME_" = "", "^HALLMARK_" = "")) %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_wrap(width = 50) # Wrap title text

    # Use provided title or generate default
    heatmap_title <- if (!is.null(title)) title else cleaned_pathway_name # Simpler default title

    # Prepare annotation data frame if provided
    if (!is.null(annotation_col)) {
        annotation_col_ordered <- annotation_col[sample_order, , drop = FALSE]
    } else {
        annotation_col_ordered <- NA # Use NA for pheatmap when no annotation
    }

    # Create heatmap using pheatmap directly
    p <- pheatmap::pheatmap(
        mat = expr_sub,
        scale = scale,
        cluster_rows = cluster_rows,
        cluster_cols = FALSE, # Keep user-defined sample order
        show_rownames = show_rownames,
        show_colnames = FALSE, # Usually too cluttered
        annotation_col = annotation_col_ordered,
        annotation_colors = annotation_colors,
        color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(50), # Slightly different red
        main = heatmap_title,
        fontsize = base_fontsize,
        fontsize_row = max(5, base_fontsize - 2), # Adjust row font size
        fontsize_col = max(5, base_fontsize - 2), # Col names not shown, but set anyway
        gaps_col = gaps_col,
        silent = TRUE, # Suppress plotting to console within function
        ... # Pass additional arguments
    )

    # Return the pheatmap object invisibly
    return(invisible(p))
}
