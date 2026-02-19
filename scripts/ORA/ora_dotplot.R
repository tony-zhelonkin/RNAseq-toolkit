# =============================================================================
# ORA Dotplot Visualization
# =============================================================================
#
# PURPOSE:
#   Create publication-quality dotplots for ORA results.
#   Uses colorblind-safe palettes by default.
#
# USAGE:
#   source("scripts/ORA/ora_dotplot.R")
#   p <- ora_dotplot(ora_result, top_n = 20)
#
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

#' Create dotplot for ORA results
#'
#' @param ora_result enrichResult object from run_ora()
#' @param top_n Number of top terms to show (default 20)
#' @param sort_by Sort by "padj" or "gene_ratio" (default "padj")
#' @param padj_cutoff Filter to terms with padj <= cutoff (default 0.05)
#' @param title Plot title (default NULL = auto-generate)
#' @param show_padj_threshold Show dashed line at padj threshold (default TRUE)
#' @param color_by Color points by "neg_log_padj" or "gene_ratio" (default "neg_log_padj")
#' @param size_by Size points by "gene_count" or "gene_ratio" (default "gene_count")
#' @param strip_prefix Remove common prefixes like "GO:" (default TRUE)
#' @param max_name_length Truncate long names (default 50)
#' @param base_size Base font size (default 11)
#' @return ggplot object
ora_dotplot <- function(ora_result,
                        top_n = 20,
                        sort_by = "padj",
                        padj_cutoff = 0.05,
                        title = NULL,
                        show_padj_threshold = TRUE,
                        color_by = "neg_log_padj",
                        size_by = "gene_count",
                        strip_prefix = TRUE,
                        max_name_length = 50,
                        base_size = 11) {

  # Handle NULL or empty results
  if (is.null(ora_result) || nrow(ora_result@result) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant enrichment",
               size = 5, hjust = 0.5) +
      theme_void() +
      ggtitle(title %||% "ORA Results")
    return(p)
  }

  # Extract and prepare data
  plot_df <- ora_result@result %>%
    as.data.frame() %>%
    dplyr::filter(p.adjust <= padj_cutoff) %>%
    dplyr::mutate(
      # Parse GeneRatio
      gene_ratio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        if (length(parts) == 2) parts[1]/parts[2] else NA
      }),
      neg_log_padj = -log10(p.adjust),
      # Clean names
      term_clean = Description
    )

  if (nrow(plot_df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("No terms with padj <=", padj_cutoff),
               size = 5, hjust = 0.5) +
      theme_void() +
      ggtitle(title %||% "ORA Results")
    return(p)
  }

  # Strip prefixes if requested
  if (strip_prefix) {
    plot_df$term_clean <- gsub("^GO:[0-9]+ ", "", plot_df$term_clean)
    plot_df$term_clean <- gsub("^KEGG_", "", plot_df$term_clean)
  }

  # Truncate long names
  plot_df$term_clean <- substr(plot_df$term_clean, 1, max_name_length)

  # Sort and select top_n
  if (sort_by == "gene_ratio") {
    plot_df <- plot_df %>%
      dplyr::arrange(desc(gene_ratio_numeric)) %>%
      dplyr::slice_head(n = top_n)
  } else {
    plot_df <- plot_df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n)
  }

  # Create factor for ordered y-axis
  plot_df$term_clean <- factor(
    plot_df$term_clean,
    levels = rev(plot_df$term_clean)
  )

  # Build aesthetic mappings
  color_var <- if (color_by == "gene_ratio") "gene_ratio_numeric" else "neg_log_padj"
  size_var <- if (size_by == "gene_ratio") "gene_ratio_numeric" else "Count"

  # Create plot
  p <- ggplot(plot_df, aes(
    x = gene_ratio_numeric,
    y = term_clean,
    color = .data[[color_var]],
    size = .data[[size_var]]
  )) +
    geom_point() +
    scale_size_continuous(
      name = if (size_by == "gene_ratio") "Gene Ratio" else "Gene Count",
      range = c(2, 8)
    ) +
    labs(
      x = "Gene Ratio",
      y = NULL,
      title = title %||% "GO Enrichment"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.y = element_text(size = base_size - 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "#EEEEEE"),
      panel.grid.minor = element_blank()
    )

  # Add color scale
  # Try to use viridis for colorblind-safe palette
  if (requireNamespace("viridis", quietly = TRUE)) {
    p <- p + viridis::scale_color_viridis(
      option = "plasma",
      name = if (color_by == "gene_ratio") "Gene Ratio" else "-log10(padj)",
      direction = 1,
      begin = 0.15,
      end = 0.95
    )
  } else {
    # Fallback to gradient
    p <- p + scale_color_gradient(
      low = "#fee0d2",
      high = "#de2d26",
      name = if (color_by == "gene_ratio") "Gene Ratio" else "-log10(padj)"
    )
  }

  return(p)
}


#' Create faceted dotplot for multiple modules
#'
#' @param ora_list Named list of enrichResult objects (one per module)
#' @param top_n_per_module Top terms per module (default 5)
#' @param ... Additional arguments passed to ora_dotplot processing
#' @return ggplot object
ora_dotplot_facet <- function(ora_list,
                              top_n_per_module = 5,
                              padj_cutoff = 0.05,
                              title = NULL,
                              base_size = 10) {

  # Combine all results
  combined_df <- tibble()

  for (module_name in names(ora_list)) {
    ora_result <- ora_list[[module_name]]

    if (is.null(ora_result) || nrow(ora_result@result) == 0) next

    df <- ora_result@result %>%
      as.data.frame() %>%
      dplyr::filter(p.adjust <= padj_cutoff) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = top_n_per_module) %>%
      dplyr::mutate(
        module = module_name,
        gene_ratio_numeric = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          if (length(parts) == 2) parts[1]/parts[2] else NA
        }),
        neg_log_padj = -log10(p.adjust),
        term_short = substr(Description, 1, 40)
      )

    combined_df <- dplyr::bind_rows(combined_df, df)
  }

  if (nrow(combined_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No enrichment found") +
             theme_void())
  }

  # Plot
  p <- ggplot(combined_df, aes(
    x = module,
    y = term_short,
    color = neg_log_padj,
    size = Count
  )) +
    geom_point() +
    scale_size_continuous(name = "Gene Count", range = c(2, 6)) +
    labs(
      x = NULL,
      y = NULL,
      title = title %||% "Module GO Enrichment"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size - 1),
      axis.text.y = element_text(size = base_size - 2),
      legend.position = "right"
    )

  # Color scale
  if (requireNamespace("viridis", quietly = TRUE)) {
    p <- p + viridis::scale_color_viridis(
      option = "plasma",
      name = "-log10(padj)",
      direction = 1,
      begin = 0.15,
      end = 0.95
    )
  } else {
    p <- p + scale_color_gradient(
      low = "#fee0d2",
      high = "#de2d26",
      name = "-log10(padj)"
    )
  }

  return(p)
}


# =============================================================================
# MESSAGE ON LOAD
# =============================================================================

message("[RNAseq-toolkit] Loaded: ora_dotplot.R (ora_dotplot, ora_dotplot_facet)")
