#' Create a Mean-Difference Plot for Differential Expression Results
#'
#' Creates a customizable mean-difference (MD) plot highlighting significant genes.
#' MD plots show the relationship between gene expression levels (x-axis) and 
#' log-fold changes (y-axis), which can reveal intensity-dependent biases.
#' 
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
    fit,
    coef,
    de_results    = NULL,
    fc_cutoff     = 1.0,
    fdr_cutoff    = 0.05,
    top_n         = 5,
    highlight_gene= NULL,
    label_method  = "top",
    max.overlaps  = 10, 
    title         = NULL,
    color_palette = c(
      "Up"   = "#D55E00",   # orange-red
      "Down" = "#0173B2",   # blue
      "NS"   = "#7F7F7F"    # grey
    ),
    show_grid     = FALSE
) {
  # ───────────────────────────────── helpers ──────────────────────────
  shade <- function(hex, factor = .6) {
    rgb <- grDevices::col2rgb(hex)/255 * factor
    grDevices::rgb(pmax(pmin(rgb,1),0)[1],
                   pmax(pmin(rgb,1),0)[2],
                   pmax(pmin(rgb,1),0)[3])
  }
  
  custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
    source("scripts/custom_minimal_theme.R", local = TRUE)
    custom_minimal_theme_with_grid
  } else {
    function() ggplot2::theme_minimal()
  }

  # ──────────────────────── extract data ────────────────────────────
  # Convert coefficient index to character if it's numeric
  if (is.numeric(coef) && coef <= length(colnames(fit))) {
    coef_name <- colnames(fit)[coef]
  } else {
    coef_name <- as.character(coef)
  }
  
  # Set default title based on coefficient name
  if (is.null(title)) {
    title <- paste("Mean-Difference Plot:", coef_name)
  }
  
  # Get DE results if not provided
  if (is.null(de_results)) {
    de_results <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
  }
  
  # Make sure we have the average expression
  if (!"AveExpr" %in% colnames(de_results) && "Amean" %in% names(fit)) {
    de_results$AveExpr <- fit$Amean
  }
  
  # Check if we have the required columns
  required_cols <- c("logFC", "AveExpr", "adj.P.Val")
  missing_cols <- required_cols[!required_cols %in% colnames(de_results)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in de_results: ", paste(missing_cols, collapse = ", "))
  }

  # ────────────────── 1. annotate significance ───────────────────────
  df <- dplyr::mutate(de_results,
                     sig_fc = abs(logFC) >= fc_cutoff,
                     sig_fdr = adj.P.Val <= fdr_cutoff,
                     status = dplyr::case_when(
                       sig_fdr & logFC > 0 ~ "Up",
                       sig_fdr & logFC < 0 ~ "Down",
                       TRUE ~ "NS"
                     ),
                     sig_both = sig_fc & sig_fdr,
                     gene_id = rownames(de_results),
                     hover_text = paste(
                       "Gene:", gene_id,
                       "<br>logFC:", round(logFC, 2),
                       "<br>Ave Expr:", round(AveExpr, 2),
                       "<br>FDR:", formatC(adj.P.Val, format = "e", digits = 2)
                     ))

  # ────────────────── 2. label selection ─────────────────────────────
  # Helper to select top_n genes by significance in each direction
  get_top <- function(direction) {
    if (direction == "up") {
      df |>
        dplyr::filter(.data$logFC > 0) |>
        dplyr::arrange(.data$adj.P.Val) |>
        dplyr::slice_head(n = top_n)
    } else {
      df |>
        dplyr::filter(.data$logFC < 0) |>
        dplyr::arrange(.data$adj.P.Val) |>
        dplyr::slice_head(n = top_n)
    }
  }

  if (label_method == "top") {
    lab_df <- dplyr::bind_rows(get_top("up"), get_top("down"))
  } else if (label_method == "sig") {
    lab_df <- df[df$sig_both, ]
  } else if (label_method == "fdr") {
    lab_df <- df[df$sig_fdr, ]
  } else if (label_method == "log2fc") {
    lab_df <- df[df$sig_fc, ]
  } else {
    lab_df <- df[0, ]
  }

  if (!is.null(highlight_gene)) {
    lab_df <- dplyr::bind_rows(lab_df,
                              df[df$gene_id %in% highlight_gene, ]) |>
              dplyr::distinct()
  }

  # ────────────────── 3. colours & limits ─────────────────────────────
  dark_pal <- vapply(color_palette, shade, character(1))
  
  # Calculate axis limits
  x_range <- range(df$AveExpr, na.rm = TRUE)
  x_padding <- (x_range[2] - x_range[1]) * 0.05
  
  y_range <- range(df$logFC, na.rm = TRUE)
  y_padding <- (y_range[2] - y_range[1]) * 0.05

  # ────────────────── 4. build ggplot ────────────────────────────────
  g <- ggplot2::ggplot(df, ggplot2::aes(AveExpr, logFC, colour = status)) +
       ggplot2::geom_point(size = 1.5, alpha = .75) +
       ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
       ggplot2::geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
       ggplot2::scale_colour_manual(name = "Significance",
         values = color_palette,
         breaks = c("Up", "Down", "NS"),
         labels = c(
           "Up"   = sprintf("Up (FDR ≤ %.2g & logFC > %.1f)", fdr_cutoff, fc_cutoff),
           "Down" = sprintf("Down (FDR ≤ %.2g & logFC < -%.1f)", fdr_cutoff, fc_cutoff),
           "NS"   = "Not significant"
         )) +
       ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
       ggplot2::coord_cartesian(
         xlim = c(x_range[1] - x_padding, x_range[2] + x_padding),
         ylim = c(y_range[1] - y_padding, y_range[2] + y_padding)
       ) +
       ggplot2::labs(
         x = "Average Expression",
         y = "log2 Fold Change",
         title = title,
         caption = sprintf("Horizontal lines: solid – logFC = 0; dashed – |logFC| = %.1f", fc_cutoff)
       ) +
       custom_minimal_theme_with_grid()

  if (!show_grid) {
    g <- g + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                           panel.grid.minor = ggplot2::element_blank())
  }

  # ────────────────── 5. labels ───────────────────────────────────────
  if (nrow(lab_df) > 0) {
    g <- g + ggrepel::geom_text_repel(
      data            = lab_df,
      ggplot2::aes(label = gene_id),
      colour          = ifelse(lab_df$gene_id %in% highlight_gene, "black",
                              dark_pal[lab_df$status]),
      fontface        = ifelse(lab_df$gene_id %in% highlight_gene, "bold", "plain"),
      size            = 3.5,
      box.padding     = .4,
      point.padding   = .3,
      max.overlaps    = max.overlaps,
      min.segment.length = 0,
      show.legend     = FALSE)
  }

  return(g)
}
