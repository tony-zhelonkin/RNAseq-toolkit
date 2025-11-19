#' Plot Per-Contrast Cross-Database Pooled Dotplot
#'
#' Generates a dotplot showing the top pathways from ALL databases for a single contrast.
#' This visualization helps identify which pathways are most enriched across different
#' gene set collections (e.g., Hallmark, KEGG, GO) for a given comparison.
#'
#' @param gsea_results_list Named list of gseaResult objects from different databases
#'                          for a single contrast
#' @param contrast_name Name of the contrast being visualized
#' @param top_n Number of top pathways to show per database
#' @param padj_cutoff Adjusted p-value cutoff for significance
#' @param output_file Path to save the plot (optional)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param show_qvalue Whether to use q-value instead of p.adjust (default: FALSE)
#'
#' @return A ggplot2 object (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have GSEA results for one contrast across multiple databases
#' plot_pooled_contrast_dotplot(
#'   gsea_results_list = all_gsea_results[["KO_vs_WT_LPS"]],
#'   contrast_name = "KO_vs_WT_LPS",
#'   top_n = 10,
#'   padj_cutoff = 0.05,
#'   output_file = "pooled_KO_vs_WT_LPS_dotplot.pdf"
#' )
#' }
#'
plot_pooled_contrast_dotplot <- function(
    gsea_results_list,
    contrast_name,
    top_n = 20,
    padj_cutoff = 0.05,
    output_file = NULL,
    width = 12,
    height = 10,
    show_qvalue = FALSE
) {

    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required for this function.")
    }
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        stop("Package 'dplyr' is required for this function.")
    }

    # Input validation
    if (is.null(gsea_results_list) || length(gsea_results_list) == 0) {
        stop("gsea_results_list is empty or NULL")
    }

    # Choose significance column
    sig_col <- if (show_qvalue) "qvalue" else "p.adjust"

    # Collect top pathways from each database
    all_pathways <- list()

    for (db_name in names(gsea_results_list)) {
        gsea_obj <- gsea_results_list[[db_name]]

        # Validate object
        if (is.null(gsea_obj) || !methods::is(gsea_obj, "gseaResult")) {
            message(sprintf("  Skipping %s (NULL or invalid object)", db_name))
            next
        }

        # Extract results
        res_df <- as.data.frame(gsea_obj@result)

        if (nrow(res_df) == 0) {
            message(sprintf("  Skipping %s (no results)", db_name))
            next
        }

        # Check if significance column exists
        if (!sig_col %in% colnames(res_df)) {
            if (sig_col == "qvalue") {
                message(sprintf("  %s: qvalue not found, falling back to p.adjust", db_name))
                sig_col <- "p.adjust"
            } else {
                message(sprintf("  Skipping %s (no p.adjust column)", db_name))
                next
            }
        }

        # Filter for significance
        sig_pathways <- res_df[res_df[[sig_col]] < padj_cutoff, ]

        if (nrow(sig_pathways) == 0) {
            message(sprintf("  Skipping %s (no significant pathways)", db_name))
            next
        }

        # Sort by absolute NES and take top N
        sig_pathways <- sig_pathways[order(abs(sig_pathways$NES), decreasing = TRUE), ]
        sig_pathways <- head(sig_pathways, top_n)

        # Add database column
        sig_pathways$Database <- db_name

        # Select relevant columns
        pathway_data <- data.frame(
            Pathway = sig_pathways$Description,
            Database = sig_pathways$Database,
            NES = sig_pathways$NES,
            pvalue = sig_pathways$pvalue,
            padj = sig_pathways[[sig_col]],
            setSize = sig_pathways$setSize,
            stringsAsFactors = FALSE
        )

        all_pathways[[db_name]] <- pathway_data
    }

    # Combine all pathways
    if (length(all_pathways) == 0) {
        stop("No significant pathways found in any database for contrast: ", contrast_name)
    }

    combined_pathways <- do.call(rbind, all_pathways)
    rownames(combined_pathways) <- NULL

    message(sprintf("  Found %d pathways across %d databases",
                   nrow(combined_pathways), length(unique(combined_pathways$Database))))

    # Create -log10(padj) for visualization
    combined_pathways$negLogPadj <- -log10(combined_pathways$padj)

    # Create a unique pathway ID combining pathway name and database
    combined_pathways$PathwayID <- paste0(combined_pathways$Database, ": ", combined_pathways$Pathway)

    # Sort pathways by NES for plotting
    combined_pathways <- combined_pathways[order(combined_pathways$NES, decreasing = TRUE), ]
    combined_pathways$PathwayID <- factor(combined_pathways$PathwayID,
                                          levels = rev(combined_pathways$PathwayID))

    # Create the plot
    p <- ggplot2::ggplot(combined_pathways,
                         ggplot2::aes(x = Database, y = PathwayID,
                                     fill = NES, size = negLogPadj)) +
        ggplot2::geom_point(shape = 21, color = "black") +
        ggplot2::scale_fill_gradient2(
            low = "#91bfdb",
            mid = "white",
            high = "#fc8d59",
            midpoint = 0,
            name = "NES"
        ) +
        ggplot2::scale_size_continuous(
            name = bquote(-log[10](.(if(show_qvalue) "q" else "p"))),
            range = c(2, 10)
        ) +
        ggplot2::labs(
            title = sprintf("Top Pathways Across Databases: %s", contrast_name),
            subtitle = sprintf("Showing top %d pathways per database (padj < %.3f)",
                             top_n, padj_cutoff),
            x = "Database",
            y = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = ggplot2::element_text(size = 8),
            panel.grid.major.y = ggplot2::element_line(color = "grey90"),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "right",
            plot.title = ggplot2::element_text(face = "bold"),
            plot.subtitle = ggplot2::element_text(size = 9, color = "grey30")
        )

    # Save if output file specified
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = p, width = width, height = height)
        message(sprintf("  Saved plot to: %s", output_file))
    }

    return(invisible(p))
}
