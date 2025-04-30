#' Create a Standard Volcano Plot for Differential Expression Results
#'
#' Create a Standard Volcano Plot for Differential Expression Results
#'
#' Creates a customizable volcano plot from DE results, highlighting significant genes.
#' Assumes `de_results` rownames are gene identifiers.
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param de_results Data frame with DE results (rownames = genes, requires 'logFC', 'P.Value').
#' @param p_cutoff Numeric, p-value cutoff for significance (default: 0.05). Note: Uses raw P.Value column.
#' @param fc_cutoff Numeric, absolute log2 fold change cutoff for significance (default: 2.0).
#' @param max.overlaps Integer, maximum number of label overlaps allowed (default: 10).
#' @param label_method Character, method to determine which genes to label:
#'        "sig" (genes significant by both p-value & fold change),
#'        "p" (genes significant by p-value only),
#'        "log2fc" (genes significant by fold change only),
#'        "none" (no gene labels).
#' @param x_breaks Numeric, interval between x-axis breaks (default: 1).
#' @param title Character, plot title (default: "Volcano Plot").
#' @param color_palette Character vector (length 4), colors for categories:
#'        NS, Log2FC only, p-value only, Both (default: c("gray", "forestgreen", "skyblue", "orange")).
#' @param highlight_gene Character vector, specific gene IDs to highlight with distinct labels (default: NULL).
#'
#' @return A ggplot2 object representing the volcano plot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2 dplyr ggrepel
#' @importFrom magrittr %>%
#'
#' @examples
#' # Assuming de_results is a data frame with columns: logFC, P.Value and gene names as rownames
#' plot_obj <- create_volcano_plot(de_results, p_cutoff = 0.05, fc_cutoff = 1.5)
#' # To save: ggsave("my_volcano.png", plot_obj)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(magrittr) # Ensure pipe operator is available

# Source the custom theme function if it exists
custom_theme_path <- file.path("scripts", "custom_minimal_theme.R")
if (file.exists(custom_theme_path)) {
  source(custom_theme_path)
} else {
  warning("Custom theme file not found at: ", custom_theme_path, ". Using default theme_minimal().")
  # Define a placeholder function if the theme is missing
  custom_minimal_theme_with_grid <- function() theme_minimal()
}


create_volcano_plot <- function(de_results, p_cutoff = 0.05, fc_cutoff = 2.0,
                                max.overlaps = 10, label_method = "sig",
                                x_breaks = 1,
                                title = "Volcano Plot",
                                color_palette = c("gray", "forestgreen", "skyblue", "orange"),
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
      name = "Significance", # Add legend title
      values = c(
        "p-value & Log2FC" = color_palette[4],
        "Log2FC" = color_palette[2],
        "p-value" = color_palette[3],
        "NS" = color_palette[1]
      ),
      labels = c( # Define labels for clarity
        "p-value & Log2FC" = sprintf("p < %.2f & |LogFC| > %.1f", p_cutoff, fc_cutoff),
        "Log2FC" = sprintf("|LogFC| > %.1f", fc_cutoff),
        "p-value" = sprintf("p < %.2f", p_cutoff),
        "NS" = "Not Significant"
      ),
      breaks = c("p-value & Log2FC", "Log2FC", "p-value", "NS") # Ensure order
    ) +
    labs(
      title = title,
      subtitle = paste("p-value cutoff:", p_cutoff, "| FC cutoff:", fc_cutoff),
      x = "Log2 Fold Change",
      y = bquote(-Log[10]~'(P-value)') # Format y-axis label
    ) +
    scale_x_continuous(breaks = seq(x_min, x_max, by = x_breaks), limits = c(x_min, x_max)) +
    custom_minimal_theme_with_grid() # Apply custom theme

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
        max.overlaps = max.overlaps,
        min.segment.length = 0 # Draw segments even if short
      )
  }

  # Add highlighted gene labels with different style (e.g., bold, slightly larger)
  if (!is.null(highlight_gene) && !is.null(highlight_data) && nrow(highlight_data) > 0) {
    base_plot <- base_plot +
      ggrepel::geom_label_repel(
        data = highlight_data,
        aes(label = rownames(highlight_data)),
        box.padding = 0.6, # Slightly more padding
        point.padding = 0.6,
        segment.color = "black",
        color = "black", # Label text color
        fill = alpha("white", 0.7), # Semi-transparent background
        fontface = "bold", # Make highlighted labels bold
        size = 4, # Slightly larger font size
        max.overlaps = Inf, # Ensure these labels are always shown
        min.segment.length = 0
      )
  }

  # Return the plot object instead of saving
  return(base_plot)
}
