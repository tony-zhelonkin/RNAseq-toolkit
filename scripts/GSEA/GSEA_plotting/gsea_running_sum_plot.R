#' Unified GSEA Running-Sum Plot
#'
#' Creates a three-panel running sum plot for GSEA results. Works with all
#' databases including MSigDB collections, SynGO, and MitoCarta.
#'
#' @param gsea_obj A `gseaResult` object from clusterProfiler/fgsea
#' @param gene_set_ids Integer vector (row indices) or character vector
#'                     (pathway IDs). Default picks top 5 by |NES|.
#' @param palette Optional color vector. If NULL, uses a 9-color vibrant palette.
#' @param labels Optional named vector of legend labels. If NULL, uses Description.
#' @param legend_pos Legend position inside the first panel (default: top-right).
#' @param base_size Base font size for the theme.
#' @param max_name_length Maximum character length for pathway names in legend.
#' @param title Optional title for the plot.
#'
#' @return A patchwork object combining three stacked panels.
#' @export
#'
#' @note CRITICAL: This function does NOT name the color palette vector.
#'       Named colors cause issues with enrichplot::gseaplot2() for custom
#'       databases like SynGO and MitoCarta. This design was validated in
#'       the original syngo_running_sum_plot() implementation.
#'
#' @note Updated 2025-12-02 to unify MSigDB and custom database support.
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids = NULL,
                                  palette = NULL,
                                  labels = NULL,
                                  legend_pos = c(.98, .98),
                                  base_size = 14,
                                  max_name_length = 40,
                                  title = NULL) {

  ## ------ 0. Package and object validation --------------------------------
  stopifnot(methods::is(gsea_obj, "gseaResult"))
  res_df <- gsea_obj@result

  if (is.null(res_df) || nrow(res_df) == 0) {
    stop("`gsea_obj` has no result rows - nothing to plot.")
  }

  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' is required.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required.")
  }

  ## ------ 1. Gene set ID selection (robust handling) ----------------------
  if (is.null(gene_set_ids)) {
    # Default: top 5 pathways by absolute NES
    ord <- order(abs(res_df$NES), decreasing = TRUE)
    gene_set_ids <- ord[seq_len(min(5L, length(ord)))]
  }

  # Convert numeric indices to pathway IDs
  if (is.numeric(gene_set_ids)) {
    # Validate indices are within bounds
    valid_idx <- gene_set_ids[gene_set_ids >= 1 & gene_set_ids <= nrow(res_df)]
    if (length(valid_idx) == 0) {
      stop("All gene_set_ids are out of bounds (max: ", nrow(res_df), ")")
    }
    if (length(valid_idx) < length(gene_set_ids)) {
      warning("Some gene_set_ids were out of bounds and removed")
    }
    gene_set_ids <- res_df$ID[valid_idx]
  } else {
    # Validate character IDs exist in results
    missing <- gene_set_ids[!gene_set_ids %in% res_df$ID]
    if (length(missing) > 0) {
      warning("Some pathway IDs not found: ", paste(missing, collapse = ", "))
      gene_set_ids <- gene_set_ids[gene_set_ids %in% res_df$ID]
    }
  }

  # Final validation
  gene_set_ids <- unique(gene_set_ids[!is.na(gene_set_ids)])
  if (length(gene_set_ids) == 0) {
    stop("No valid gene-set IDs - nothing to plot.")
  }

  n_sets <- length(gene_set_ids)

  ## ------ 2. Color palette (CRITICAL: do NOT name the colors) -------------
  # Named color vectors cause issues with enrichplot::gseaplot2() for custom

  # databases like SynGO and MitoCarta. Keep colors as unnamed vector.
  if (is.null(palette)) {
    base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                  "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
    if (n_sets > length(base_pal)) {
      # Interpolate additional colors for large pathway sets
      palette <- grDevices::colorRampPalette(base_pal)(n_sets)
    } else {
      palette <- base_pal[seq_len(n_sets)]
    }
  } else if (length(palette) < n_sets) {
    palette <- rep(palette, length.out = n_sets)
  }
  # NOTE: Explicitly NOT naming the palette - this is intentional!
  # names(palette) <- gene_set_ids  # DO NOT DO THIS

  ## ------ 3. Legend labels (truncated for readability) --------------------
  if (is.null(labels)) {
    labels <- res_df$Description[match(gene_set_ids, res_df$ID)]
  }
  # Truncate long labels
  labels <- vapply(labels, function(x) {
    if (is.na(x)) return("(Unknown)")
    if (nchar(x) > max_name_length) {
      paste0(substr(x, 1, max_name_length - 3), "...")
    } else {
      x
    }
  }, FUN.VALUE = character(1))
  names(labels) <- gene_set_ids

  ## ------ 4. Generate plot via enrichplot ---------------------------------
  tryCatch({
    # Set title if not provided
    plot_title <- if (!is.null(title)) title else "Gene Set Enrichment"

    p_raw <- enrichplot::gseaplot2(
      gsea_obj,
      geneSetID = gene_set_ids,
      title = plot_title,
      subplots = c(1, 2, 3),
      pvalue_table = FALSE,
      rel_heights = c(1.5, 0.5, 0.5),
      color = palette  # Pass unnamed color vector
    )

    ## ------ 5. Styling helper function ------------------------------------
    stylise <- function(p, show_x, show_y, show_legend) {
      p +
        ggplot2::labs(color = NULL) +
        ggplot2::theme_classic(base_size = base_size) +
        ggplot2::theme(
          legend.position = if (show_legend) legend_pos else "none",
          legend.justification = c(1, 1),
          legend.background = ggplot2::element_rect(fill = "white", colour = "grey90"),
          legend.margin = ggplot2::margin(5, 5, 5, 5),
          legend.key.size = ggplot2::unit(0.8, "lines"),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
          axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
          axis.title.x = if (show_x) ggplot2::element_text() else ggplot2::element_blank(),
          axis.title.y = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
          axis.text.x = if (show_x) ggplot2::element_text() else ggplot2::element_blank(),
          axis.text.y = if (show_y) ggplot2::element_text() else ggplot2::element_blank(),
          plot.margin = ggplot2::margin(5, 10, 5, 5)
        )
    }

    # Apply styling to each panel
    p1 <- stylise(p_raw[[1]], show_x = FALSE, show_y = TRUE, show_legend = TRUE)
    p2 <- stylise(p_raw[[2]], show_x = FALSE, show_y = FALSE, show_legend = FALSE)
    p3 <- stylise(p_raw[[3]], show_x = TRUE, show_y = TRUE, show_legend = FALSE)

    ## ------ 6. Legend styling ----------------------------------------------
    # Note: enrichplot::gseaplot2() already uses Description for legend labels.
    # We just ensure colors are properly displayed.
    p1 <- p1 +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(colour = palette)
        )
      )

    ## ------ 7. Combine panels with patchwork ------------------------------
    patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(2, 0.5, 0.5))

  }, error = function(e) {
    stop("Error generating running sum plot: ", e$message)
  })
}
