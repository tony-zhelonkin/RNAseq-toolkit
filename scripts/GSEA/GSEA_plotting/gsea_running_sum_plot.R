#' Running-sum plot
#'
#' @param gsea_obj   gseaResult
#' @param gene_set_ids numeric or character; 1–5 recommended
#' @param palette    named vector of colours (optional)
#' @param labels     named vector of legend labels (optional)
#' @param legend_pos legend position on the first panel
#' @param base_size  base font size
#' @return patchwork object ready for ggsave()
#' @export
gsea_running_sum_plot <- function(gsea_obj,
                                  gene_set_ids  = 1:5,
                                  palette       = NULL,
                                  labels        = NULL,
                                  legend_pos    = c(.65, .8),
                                  base_size     = 14) {

  stopifnot(methods::is(gsea_obj, "gseaResult"))
  if (!requireNamespace("enrichplot", quietly = TRUE))
      stop("Package 'enrichplot' is required.")
  if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Package 'patchwork' is required.")

  ## ------------------------------------------------------------------ ##
  ## 1  colour / label mapping                                          ##
  ## ------------------------------------------------------------------ ##
  path_ids <- if (is.numeric(gene_set_ids)) gsea_obj@result$ID[gene_set_ids]
              else                          gene_set_ids
  n_sets    <- length(path_ids)

  if (is.null(palette))
      palette <- c("#648FFF", "#785EF0", "#DC267F", "#FE6100",
                   "#FFB000")[seq_len(n_sets)]
  if (is.null(labels))
      labels <- gsea_obj@result$Description[match(path_ids,
                                                  gsea_obj@result$ID)]

  names(palette) <- names(labels) <- path_ids

  ## ------------------------------------------------------------------ ##
  ## 2  raw plots from enrichplot                                       ##
  ## ------------------------------------------------------------------ ##
  p_raw <- enrichplot::gseaplot2(gsea_obj,
                                 geneSetID    = gene_set_ids,
                                 subplots     = c(1, 2, 3),
                                 pvalue_table = FALSE,
                                 rel_heights  = c(1.5, .5, .5))

  ## ------------------------------------------------------------------ ##
  ## 3  styling helper                                                  ##
  ## ------------------------------------------------------------------ ##
  stylise <- function(g, show_x, show_y, show_legend) {
  g +
    scale_color_manual(values = palette, labels = labels) +
    theme_minimal(base_size = base_size) +
    theme(
      ## ── legend ───────────────────────────────────────────────
      legend.position      = if (show_legend) c(0.98, 0.98) else "none",  # top-right
      legend.justification = c(1, 1),   # anchor exactly in the corner
      legend.title         = element_blank(),
      ## ── axes ─────────────────────────────────────────────────
      axis.title.x = if (show_x) element_text(size = base_size) else element_blank(),
      axis.title.y = if (show_y) element_text(size = base_size) else element_blank(),
      axis.text.x  = if (show_x) element_text(size = base_size * .9) else element_blank(),
      axis.text.y  = if (show_y) element_text(size = base_size * .9) else element_blank(),
      ## ── grid lines ──────────────────────────────────────────
      panel.grid.major.x = element_line(colour = "grey80"),  # keep only the major X grid
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}


  ## ------------------------------------------------------------------ ##
  ## 4  apply styling + stack with patchwork                            ##
  ## ------------------------------------------------------------------ ##
  p1 <- stylise(p_raw[[1]], show_x = FALSE, show_y = TRUE,  show_legend = TRUE)
  p2 <- stylise(p_raw[[2]], show_x = FALSE, show_y = FALSE, show_legend = FALSE)
  p3 <- stylise(p_raw[[3]], show_x = TRUE,  show_y = TRUE,  show_legend = FALSE)

  combined <- patchwork::wrap_plots(p1, p2, p3,
                                    ncol    = 1,
                                    heights = c(2, 1, 1))

  return(combined)   # ONE grob — safe for ggsave()
}
