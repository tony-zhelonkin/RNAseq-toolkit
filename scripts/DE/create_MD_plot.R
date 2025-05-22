#' Create a Mean-Difference Plot for Differential Expression Results
#'
# -----------------------------------------------------------------------------
#  create_MD_plot.R   –   compact, publication‑ready mean–difference (MD) plot
# -----------------------------------------------------------------------------
#  This helper extends the classic limma::plotMD idea with a few additions that
#  make systematic shifts pop out instantly while keeping per‑gene resolution.
#  Key features
#  ─────────────
#    • Background "grey cloud" for non‑significant (NS) genes – light & subtle.
#    • Coloured, opaque dots for FDR‑significant genes (Up / Down).
#    • LOESS trend line (dashed) to expose global bias across expression range.
#    • Vertical guide at median expression – catches normalisation issues.
#    • Optional quadrant counts to quantify directionality at a glance.
#    • Smart label selection: top‑N by FDR per side, or any user‑defined list.
#    • Colour‑blind‑safe Okabe‑Ito palette by default.
#' @param fit An MArrayLM object from limma containing the fitted model.
#' @param coef Integer or character, indicating which coefficient/contrast to plot.
#' @param de_results Data frame with DE results matching the fit object (requires 'adj.P.Val').
#'        If NULL, will extract results from fit object.
#' @param fc_cutoff Numeric, absolute log2 fold change cutoff for highlighting (default: 1.0).
#' @param fdr_cutoff Numeric, FDR cutoff for significance (default: 0.05).
#' @param top_n Integer, number of top genes by significance to label on each side (default: 5).
#' @param highlight_gene Character vector, specific gene IDs to highlight (default: NULL).
#' @param label_method Character, method to determine which genes to label:
#'        "top" (top_n genes by significance in each direction),
#'        "sig" (genes significant by both FDR & fold change),
#'        "fdr" (genes significant by FDR only),
#'        "log2fc" (genes significant by fold change only),
#'        "none" (no gene labels).
#' @param max.overlaps Integer, maximum number of label overlaps allowed (default: 10).
#' @param title Character, plot title (default: based on coefficient name).
#' @param color_palette Character vector of colors for up, down, and non-significant.
#' @param show_grid Logical, whether to display grid lines (default: FALSE).
#'
#' @return A ggplot2 object representing the MD plot.
#' @export
#' @import ggplot2 dplyr ggrepel
#'
#' @examples
#' # Basic usage with fit object and default parameters
#' plot_obj <- create_MD_plot(fit, coef = "group_treatment_vs_control")

create_MD_plot <- function(
    fit,                            # MArrayLM from limma
    coef,                           # index or name of contrast
    de_results     = NULL,          # data.frame from topTable (optional)
    fc_cutoff      = 1,             # abs(log2FC) threshold for significance
    fdr_cutoff     = 0.05,          # FDR threshold for significance
    top_n          = 5,             # how many top genes to label per side
    highlight_gene = NULL,          # vector of gene IDs to emphasise
    label_method   = "top",         # "top" | "sig" | "fdr" | "log2fc" | "none"
    max.overlaps   = 10,            # ggrepel argument
    title          = NULL,          # custom title (defaults to coef name)
    color_palette  = c(             # Okabe–Ito
      Up   = "#D55E00",            # orange‑red
      Down = "#0072B2",            # blue
      NS   = "#999999"             # grey
    ),
    show_grid      = FALSE,         # background grid toggle
    show_quadrant_counts = TRUE)    # annotate #points per quadrant
{
  # -------------------- helper::theme ---------------------------------------
  custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
    source("scripts/custom_minimal_theme.R", local = TRUE)
    custom_minimal_theme_with_grid
  } else {
    function() ggplot2::theme_minimal()
  }

  # -------------------- helper::shade / text‑colour -------------------------
  shade <- function(hex, factor = .6) {
    rgb <- grDevices::col2rgb(hex)/255 * factor
    grDevices::rgb(pmax(pmin(rgb,1),0)[1],
                   pmax(pmin(rgb,1),0)[2],
                   pmax(pmin(rgb,1),0)[3])
  }
  text_col <- function(hex) {                       # black on light, white on dark
    rgb <- grDevices::col2rgb(hex)/255
    lum <- 0.299*rgb[1] + 0.587*rgb[2] + 0.114*rgb[3]
    ifelse(lum > .55, "black", "white")
  }

  # -------------------- extract coefficient name ---------------------------
  coef_name <- if (is.numeric(coef)) colnames(fit)[coef] else as.character(coef)
  if (is.null(title)) title <- paste("MD plot:", coef_name)

  # -------------------- fetch DE table -------------------------------------
  if (is.null(de_results)) {
    de_results <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  }
  if (!"AveExpr" %in% colnames(de_results) && !is.null(fit$Amean))
      de_results$AveExpr <- fit$Amean

  req <- c("logFC", "AveExpr", "adj.P.Val")
  stopifnot(all(req %in% colnames(de_results)))

  # -------------------- annotate significance ------------------------------
  df <- dplyr::mutate(de_results,
        sig_fc  = abs(logFC) >= fc_cutoff,
        sig_fdr = adj.P.Val  <= fdr_cutoff,
        status  = dplyr::case_when(sig_fdr & logFC > 0 ~ "Up",
                                   sig_fdr & logFC < 0 ~ "Down",
                                   TRUE                ~ "NS"),
        gene_id = rownames(de_results))

  # -------------------- choose labels --------------------------------------
  get_top <- function(direction) {
    dplyr::filter(df, (direction=="up" & logFC>0) | (direction=="down" & logFC<0)) |>
      dplyr::arrange(adj.P.Val) |>
      dplyr::slice_head(n = top_n)
  }
  lab_df <- switch(label_method,
      top    = dplyr::bind_rows(get_top("up"), get_top("down")),
      sig    = dplyr::filter(df, sig_fc & sig_fdr),
      fdr    = dplyr::filter(df, sig_fdr),
      log2fc = dplyr::filter(df, sig_fc),
      none   = df[0,],
      df[0,])
  if (!is.null(highlight_gene))
      lab_df <- dplyr::bind_rows(lab_df,
                 dplyr::filter(df, gene_id %in% highlight_gene)) |>
                 dplyr::distinct()

  # -------------------- split NS / sig for two‑layer plotting --------------
  df_ns  <- dplyr::filter(df, status == "NS")
  df_sig <- dplyr::filter(df, status != "NS")

  # axis limits -------------------------------------------------------------
  x_pad <- diff(range(df$AveExpr))*0.05
  y_pad <- diff(range(df$logFC)) *0.05

  # dark palette for label text
  dark_pal <- vapply(color_palette, shade, character(1))

  # -------------------- build ggplot ---------------------------------------
  g <- ggplot2::ggplot() +

       ## 1  background cloud
       ggplot2::geom_point(data = df_ns,
                           aes(AveExpr, logFC),
                           colour = color_palette["NS"],
                           size = 1.3, alpha = .25, shape = 16) +

       ## 2  significant dots
       ggplot2::geom_point(data = df_sig,
                           aes(AveExpr, logFC, colour = status),
                           size = 2.2, alpha = .9, shape = 16) +

       ## guides & trend ----------------------------------------------------
       ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
       ggplot2::geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       ggplot2::geom_smooth(data = df, aes(AveExpr, logFC),
                            method = "loess", se = FALSE,
                            colour = "red", linetype = 2, linewidth = .6) +
       ggplot2::geom_vline(xintercept = median(df$AveExpr), colour = "grey60", linetype = 3) +

       ## colour scale ------------------------------------------------------
       ggplot2::scale_colour_manual(values = color_palette, name = "Significance",
            breaks = c("Up","Down"),
            labels = c(sprintf("Up  (FDR ≤ %.2g)", fdr_cutoff),
                       sprintf("Down (FDR ≤ %.2g)", fdr_cutoff))) +

       ggplot2::coord_cartesian(xlim = range(df$AveExpr)+c(-x_pad,x_pad),
                                ylim = range(df$logFC)+c(-y_pad,y_pad)) +
       ggplot2::labs(x = "Average Expression",
                     y = "log2 Fold Change",
                     title = title,
                     caption = sprintf("Dashed: |logFC| ≥ %.1f", fc_cutoff)) +
       custom_minimal_theme_with_grid()

  if (!show_grid) {
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank())
  }

    # quadrant counts ---------------------------------------------------------
  if (show_quadrant_counts && nrow(df_sig)) {
    qx <- median(df$AveExpr); qy <- 0
    quad <- data.frame(
      x     = c(min(df$AveExpr), min(df$AveExpr), max(df$AveExpr), max(df$AveExpr)),
      y     = c(max(df$logFC),   min(df$logFC),   max(df$logFC),   min(df$logFC)),
      hjust = c(0, 0, 1, 1),
      vjust = c(1, 0, 1, 0),
      lbl   = c(sum(df_sig$logFC > 0 & df_sig$AveExpr < qx),
                sum(df_sig$logFC < 0 & df_sig$AveExpr < qx),
                sum(df_sig$logFC > 0 & df_sig$AveExpr > qx),
                sum(df_sig$logFC < 0 & df_sig$AveExpr > qx))
    )
    g <- g + ggplot2::geom_text(
      data = quad,
      mapping = ggplot2::aes(x = x, y = y, label = lbl, hjust = hjust, vjust = vjust),
      colour = "grey30", 
      size = 3
    )
  }

  # labels ------------------------------------------------------------------
  if (nrow(lab_df))
    g <- g + ggrepel::geom_text_repel(
            data  = lab_df,
            aes(AveExpr, logFC, label = gene_id),
            colour    = ifelse(lab_df$gene_id %in% highlight_gene, "black",
                               dark_pal[lab_df$status]),
            fontface  = ifelse(lab_df$gene_id %in% highlight_gene, "bold", "plain"),
            size      = 3.3,
            box.padding = .35, point.padding = .25,
            max.overlaps = max.overlaps,
            show.legend  = FALSE)

  return(g)
}
