#' Create a Standard Volcano Plot for Differential Expression Results
#'
#' This function creates a customizable volcano plot from differential expression results,
#' highlighting significant genes based on p-value and fold change thresholds.
#'
#' @param de_results A data frame containing differential expression results with at least
#'        the following columns: gene identifiers as rownames, logFC, and P.Value.
#' @param p_cutoff Numeric, p-value cutoff for significance (default: 0.05).
#' @param fc_cutoff Numeric, fold change cutoff for significance (default: 2.0).
#' @param max.overlaps Integer, maximum number of label overlaps allowed (default: 10).
#' @param label_method Character, method to determine which genes to label:
#'        "sig" (genes significant by both p-value and fold change),
#'        "p" (genes significant by p-value),
#'        "log2fc" (genes significant by fold change),
#'        "none" (no gene labels).
#' @param x_breaks Numeric, interval between x-axis breaks (default: 1).
#' @param title Character, plot title (default: "Volcano Plot").
#' @param color_pallette Character vector, colors for different gene categories in order:
#'        non-significant, significant by fold change, significant by p-value,
#'        significant by both (default: c("gray", "forestgreen", "skyblue", "orange")).
#' @param highlight_gene Character vector, specific gene IDs to highlight (default: NULL).
#'
#' @return A ggplot2 object representing the volcano plot.
#' @export
#'
#' @examples
#' # Assuming de_results is a data frame with columns: logFC, P.Value
#' create_volcano_plot(de_results, p_cutoff = 0.05, fc_cutoff = 1.5)
library(ggplot2)
library(dplyr)
library(ggrepel)

create_volcano_plot <- function(de_results, p_cutoff = 0.05, fc_cutoff = 2.0,
                                max.overlaps = 10, label_method = "sig",
                                x_breaks = 1,
                                title = "Volcano Plot",
                                color_pallette = c("gray", "forestgreen", "skyblue", "orange"),
                                highlight_gene = NULL) {
  # Original data processing remains the same
  de_results <- de_results %>%
    mutate(
      significant_fc = abs(logFC) > fc_cutoff,
      significant_p = P.Value < p_cutoff,
      highlight = significant_fc & significant_p,
      hover_text = paste(
        "Gene:", rownames(de_results),
        "<br>LogFC:", round(logFC, 2),
        "<br>P-value:", formatC(P.Value, format = "e", digits = 2)
      ),
      color = case_when(
        highlight ~ "p-value & Log2FC",
        significant_fc ~ "Log2FC",
        significant_p ~ "p-value",
        TRUE ~ "NS"
      )
    )

  # Determine which genes to label based on label_method
  label_data <- switch(label_method,
    "sig" = de_results[de_results$highlight, ],
    "p" = de_results[de_results$significant_p, ],
    "log2fc" = de_results[de_results$significant_fc, ],
    "none" = NULL
  )

  # Create highlight_data for specifically requested genes
  if (!is.null(highlight_gene)) {
    highlight_data <- de_results[rownames(de_results) %in% highlight_gene, ]
  }

  # Rest of the plotting code remains the same until the ggrepel section
  x_min <- floor(min(de_results$logFC) / x_breaks) * x_breaks
  x_max <- ceiling(max(de_results$logFC) / x_breaks) * x_breaks

  base_plot <- ggplot(de_results, aes(x = logFC, y = -log10(P.Value), color = color, text = hover_text)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
    scale_color_manual(
      values = c(
        "p-value & Log2FC" = color_pallette[4],
        "Log2FC" = color_pallette[2],
        "p-value" = color_pallette[3],
        "NS" = color_pallette[1]
      ),
      breaks = c("p-value & Log2FC", "Log2FC", "p-value", "NS")
    ) +
    labs(
      title = title,
      subtitle = paste("p-value cutoff:", p_cutoff, "| FC cutoff:", fc_cutoff),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    scale_x_continuous(breaks = seq(x_min, x_max, by = x_breaks), limits = c(x_min, x_max)) +
    custom_minimal_theme_with_grid() +
    theme(
      legend.title = element_blank()
    )

  # Add regular labels based on label_method
  if (!is.null(label_data) && nrow(label_data) > 0) {
    base_plot <- base_plot +
      ggrepel::geom_label_repel(
        data = label_data,
        aes(label = rownames(label_data)),
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "black",
        color = "black",
        size = 3.5,
        max.overlaps = max.overlaps
      )
  }

  # Add highlighted gene labels with different style
  if (!is.null(highlight_gene) && nrow(highlight_data) > 0) {
    base_plot <- base_plot +
      ggrepel::geom_label_repel(
        data = highlight_data,
        aes(label = rownames(highlight_data)),
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "#000000",
        color = "#000000",
        fontface = "plain",
        size = 4,
        max.overlaps = Inf
      )
  }

  # Save the plot
  filename <- paste0("3_Results/imgs/volcano/", gsub(" ", "_", title), ".pdf")
  ggsave(
    filename = filename,
    plot = base_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  return(base_plot)
}
