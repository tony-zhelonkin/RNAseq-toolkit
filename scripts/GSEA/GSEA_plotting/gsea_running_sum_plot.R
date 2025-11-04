#' Safe GSEA running-sum plot (handles 1-N gene-sets)
#'
#' @param gsea_obj     A `gseaResult` object (clusterProfiler)
#' @param gene_set_ids Integer vector (row indices) or character vector
#'                     (gs_name/ID). Default picks the top 5 by |NES|.
#' @param palette      Optional named colour vector.
#' @param labels       Optional named vector of legend labels.
#' @param legend_pos   Legend position inside the first panel.
#' @param base_size    Base font size.
#' @return             A patchwork object (three stacked panels).
#' @export
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids = NULL,
                                  palette      = NULL,
                                  labels       = NULL,
                                  legend_pos   = c(.98, .98),
                                  base_size    = 14,
                                  max_name_length = 40) {

  ## ------ 0. sanity checks (unchanged) ---------------------------------
  stopifnot(methods::is(gsea_obj, "gseaResult"))
  res_df <- gsea_obj@result
  if (is.null(res_df) || nrow(res_df) == 0)
        stop("`gsea_obj` has no result rows – nothing to plot.")
  if (!requireNamespace("enrichplot", quietly = TRUE))
        stop("Package 'enrichplot' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE))
        stop("Package 'patchwork' is required.")

  ## ------ 1. gene-set selection (unchanged, bullet-proof) --------------
  if (is.null(gene_set_ids)) {
      ord <- order(abs(res_df$NES), decreasing = TRUE)
      gene_set_ids <- ord[seq_len(min(5L, length(ord)))]
  }
  if (is.numeric(gene_set_ids)) {
      gene_set_ids <- gene_set_ids[gene_set_ids >= 1 &
                                   gene_set_ids <= nrow(res_df)]
      gene_set_ids <- res_df$ID[gene_set_ids]
  }
  gene_set_ids <- unique(gene_set_ids[!is.na(gene_set_ids)])
  if (!length(gene_set_ids))
      stop("No valid gene-set IDs supplied – nothing to plot.")

  ## ------ 2. palette & legend labels -----------------------------------
  n_sets <- length(gene_set_ids)
  if (is.null(palette)) {
      base_pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                    "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")
      palette  <- base_pal[seq_len(n_sets)]
  } else if (length(palette) < n_sets) {
      palette <- rep(palette, length.out = n_sets)
  }
  names(palette) <- gene_set_ids

  if (is.null(labels)) {
      labels <- res_df$Description[match(gene_set_ids, res_df$ID)]
  }
  ## truncate long labels like in your SynGO plot
  labels <- vapply(labels, function(x) {
      if (nchar(x) > max_name_length)
          paste0(substr(x, 1, max_name_length - 3), "…")
      else x
  }, FUN.VALUE = character(1))
  names(labels) <- gene_set_ids

  ## ------ 3. build raw panels ------------------------------------------
  p_raw <- enrichplot::gseaplot2(gsea_obj,
                                 geneSetID    = gene_set_ids,
                                 subplots     = c(1, 2, 3),
                                 pvalue_table = FALSE,
                                 rel_heights  = c(1.5, .5, .5),
                                 color        = palette)

  ## ------ 4. styling helper  -------------------------------------------
  stylise <- function(p, show_x, show_y, show_legend) {
      p +
        labs(color = NULL) +
        theme_classic(base_size = base_size) +
        theme(
            legend.position      = if (show_legend) legend_pos else "none",
            legend.justification = c(1, 1),
            legend.background    = element_rect(fill = "white", colour = "grey90"),
            legend.margin        = margin(5, 5, 5, 5),
            legend.key.size      = unit(0.8, "lines"),
            axis.line            = element_line(colour = "black", linewidth = .4),
            axis.ticks           = element_line(colour = "black", linewidth = .4),
            axis.title.x = if (show_x) element_text() else element_blank(),
            axis.title.y = if (show_y) element_text() else element_blank(),
            axis.text.x  = if (show_x) element_text() else element_blank(),
            axis.text.y  = if (show_y) element_text() else element_blank(),
            plot.margin  = margin(5, 10, 5, 5)
        )
  }

  p1 <- stylise(p_raw[[1]], FALSE, TRUE,  TRUE)
  p2 <- stylise(p_raw[[2]], FALSE, FALSE, FALSE)
  p3 <- stylise(p_raw[[3]], TRUE,  TRUE,  FALSE)

  ## add legend labels via guides() so they always align with colours
  p1 <- p1 +
        guides(color = guide_legend(override.aes = list(colour = palette),
                                    labels = labels))

  ## ------ 5. stack with patchwork --------------------------------------
  patchwork::wrap_plots(p1, p2, p3, ncol = 1, heights = c(2, 1, 1))
}