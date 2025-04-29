#' Multi-Style Pathway Volcano Plot Function
#'
#' Multi-Style Pathway Volcano Plot Function
#'
#' Creates customized volcano plots highlighting genes from a specific GSEA pathway
#' within differential expression results. Offers multiple styling options.
#' Assumes `de_results` rownames are gene identifiers matching GSEA results.
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param pathway_name Character, name of the pathway to highlight (must match GSEA 'Description').
#' @param gsea_results A GSEA results object containing pathway information.
#' @param de_results Data frame with DE results (rownames = genes, requires 'logFC', 'adj.P.Val').
#' @param p_cutoff Numeric, adjusted p-value cutoff for significance (default: 0.05).
#' @param fc_cutoff Numeric, absolute log2 fold change cutoff for significance (default: 2.0).
#' @param label_method Character, method to determine which genes to label:
#'        "default" (pathway genes significant by both p-value & Log2FC),
#'        "fc" (pathway genes significant by Log2FC),
#'        "p" (pathway genes significant by p-value),
#'        "all" (pathway genes significant by either).
#' @param max_overlaps Integer, max label overlaps allowed for `ggrepel` (default: 100).
#' @param style Character, plotting style: "clean", "claude", or "gpt" (default: "clean").
#'
#' @return A ggplot2 object representing the volcano plot.
#' @export
#' @import ggplot2 dplyr stringr ggrepel
#' @importFrom magrittr %>%
#'
#' @examples
#' # Assuming gsea_results is a GSEA result object and de_results contains DE analysis results
#' analyze_pathway_volcano("HALLMARK_APOPTOSIS", gsea_results, de_results,
#'                         p_cutoff = 0.05, fc_cutoff = 1.5, style = "clean")
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)
library(magrittr) # Ensure pipe operator is available

# Source the custom theme function if it exists
custom_theme_path <- file.path("scripts", "custom_minimal_theme.R")
if (file.exists(custom_theme_path)) {
  source(custom_theme_path)
} else {
  warning("Custom theme file not found at: ", custom_theme_path)
  # Define a placeholder function if the theme is missing
  custom_minimal_theme_with_grid <- function() theme_minimal()
}


# Main function to generate the volcano plot with different styles
analyze_pathway_volcano <- function(pathway_name, gsea_results, de_results,
                                    p_cutoff = 0.05, fc_cutoff = 2.0,
                                    label_method = "default", max_overlaps = 100, style = 'clean') {
  # Step 1: Extract the relevant genes for the specified pathway
  message("Extracting genes for the pathway: ", pathway_name)
  pathway_genes <- as.data.frame(gsea_results) %>%
    filter(Description == pathway_name) %>%
    pull(core_enrichment) %>% 
    str_split("/", simplify = TRUE) %>% 
    as.vector()

  # Step 2: Prepare DE results and identify pathway genes
  message("Preparing DE results...")
  de_results <- de_results %>%
    mutate(
      in_pathway = rownames(.) %in% pathway_genes,
      significant_fc = abs(logFC) > fc_cutoff,
      significant_p = adj.P.Val < p_cutoff,
      highlight = significant_fc & significant_p
    ) %>%
    # Add a color column based on the condition
    mutate(
      color = case_when(
        highlight ~ "p-value & Log2FC",   # Both significant fold change and p-value
        significant_fc ~ "Log2FC",        # Significant fold change only
        significant_p ~ "p-value",        # Significant p-value only
        TRUE ~ "NS"                       # Not significant
      )
    )
  
  # Determine labeling data based on method
  label_data <- switch(
    label_method,
    "default" = de_results %>% filter(in_pathway & highlight),  # Only label pathway genes with both p-value & Log2FC significant
    "fc" = de_results %>% filter(in_pathway & abs(logFC) > fc_cutoff),  # Label pathway genes crossing fold change threshold
    "p" = de_results %>% filter(in_pathway & adj.P.Val < p_cutoff),  # Label pathway genes crossing p-value threshold
    "all" = de_results %>% filter(in_pathway & (abs(logFC) > fc_cutoff | adj.P.Val < p_cutoff)),  # Label all significant pathway genes (either p-value or fold change)
    stop("Invalid label_method. Please use 'default', 'fc', 'p', or 'all'.")
  )

  # Define base plot
  base_plot <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    labs(
      title = paste("Volcano Plot -", pathway_name),
      subtitle = paste("p-value cutoff:", p_cutoff, "| FC cutoff:", fc_cutoff),
      x = "Log2 Fold Change",
      y = "-Log10(Adj. P-value)"
    ) +
    scale_x_continuous(breaks = seq(floor(min(de_results$logFC)), ceiling(max(de_results$logFC)), by = 2))
  
  # Choose plot style based on user input
  plot <- switch(
    style,
    'clean' = plot_style_clean(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    'claude' = plot_style_claude(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    'gpt' = plot_style_gpt(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps),
    stop("Invalid style. Please use 'clean', 'claude', or 'gpt'.")
  )
  
  return(plot)
}


#' Plotting Style: Clean
#'
#' Generates a volcano plot using the 'clean' style. Highlights significant
#' pathway genes based on p-value and fold change. Non-pathway genes are greyed out.
#' Minimal legend.
#'
#' @param base_plot ggplot object, the base volcano plot structure.
#' @param de_results Data frame, prepared DE results with pathway/significance info.
#' @param pathway_genes Character vector, genes belonging to the target pathway.
#' @param label_data Data frame, subset of `de_results` for genes to be labeled.
#' @param p_cutoff Numeric, p-value cutoff used.
#' @param fc_cutoff Numeric, fold change cutoff used.
#' @param max_overlaps Integer, max label overlaps for `ggrepel`.
#' @return A ggplot object with the 'clean' style applied.
#' @keywords internal
plot_style_clean <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  # Define the colors for each condition
  colors <- c(
    "NS" = "#999999",               # Not significant (used for non-pathway genes)
    "p-value & Log2FC" = "#E69F00", # Pathway genes: Both p-value and Log2FC significant
    "p-value" = "#56B4E9"           # Pathway genes: Only p-value significant
  )

  # Filter for pathway-specific data
  pathway_data <- de_results %>% filter(in_pathway)
  # Base plot modifications for clean style
  plot <- base_plot +
    # Non-pathway genes in grey
    geom_point(data = de_results %>% filter(!in_pathway),
               aes(color = "NS"), alpha = 0.6, size = 1) +

    # Pathway genes: highlight based on p-value and fold change significance
    geom_point(data = pathway_data %>% filter(highlight),
               aes(color = "p-value & Log2FC"), alpha = 0.7, size = 1.5) +

    # Pathway genes: highlight based only on p-value significance
    geom_point(data = pathway_data %>% filter(!highlight & significant_p),
               aes(color = "p-value"), alpha = 0.7, size = 1.5) +

    # Set color scale manually
    scale_color_manual(values = colors, name = "Pathway Gene Status") +

    # Apply the custom minimal theme
    custom_minimal_theme_with_grid() +
    theme(legend.position = "none") # Remove the legend

  # Add labels if needed
  if (nrow(label_data) > 0) {
    plot <- plot + geom_label_repel(
      data = label_data,
      aes(label = rownames(label_data)),
      box.padding = 0.5,
      point.padding = 0.5,
      force = 5,
      segment.color = 'grey50',
      max.overlaps = max_overlaps
    )
  }
  
  return(plot)
}


#' Plotting Style: Claude
#'
#' Generates a volcano plot using the 'claude' style. Differentiates pathway
#' genes based on significance status (NS, p-value, Log2FC, both). Includes
#' a detailed legend positioned inside the plot.
#'
#' @inheritParams plot_style_clean
#' @return A ggplot object with the 'claude' style applied.
#' @keywords internal
plot_style_claude <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  colors <- c(
    "NS" = "#999999",               # Pathway gene: Not significant
    "p-value" = "#56B4E9",          # Pathway gene: Significant p-value only
    "Log2FC" = "#0072B2",           # Pathway gene: Significant Log2FC only
    "p-value & Log2FC" = "#E69F00"  # Pathway gene: Significant both
  )

  plot <- base_plot +
    # Non-pathway genes (lighter grey)
    geom_point(data = de_results %>% filter(!in_pathway),
               color = "#CCCCCC", alpha = 0.4, size = 0.6) +
    # Pathway genes colored by significance status
    geom_point(data = de_results %>% filter(in_pathway),
               aes(color = color), alpha = 0.7, size = 1.5) +
    # Manual color scale and legend setup
    scale_color_manual(values = colors,
                       name = "Pathway Gene Significance",
                       labels = c("NS" = "NS",
                                  "p-value" = sprintf("p < %.2f", p_cutoff),
                                  "Log2FC" = sprintf("|LogFC| > %.1f", fc_cutoff),
                                  "p-value & Log2FC" = "Both")) +
    # Apply custom theme and adjust legend position/appearance
    custom_minimal_theme_with_grid() +
    theme(
      legend.position = c(0.95, 0.95), # Position legend top-right
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "grey80", size = 0.5), # Box around legend
      legend.key = element_rect(fill = "white", color = NA) # Ensure legend keys have white background
    )

  # Add labels if needed
  if (nrow(label_data) > 0) {
    plot <- plot + geom_label_repel(
      data = label_data,
      aes(label = rownames(label_data)),
      box.padding = 0.5,
      point.padding = 0.5,
      force = 5,
      segment.color = 'grey50',
      max.overlaps = max_overlaps
    )
  }
  
  return(plot)
}


#' Plotting Style: GPT
#'
#' Generates a volcano plot using the 'gpt' style. Similar coloring to 'claude'
#' but uses alpha transparency for points and overrides alpha in the legend
#' for clarity.
#'
#' @inheritParams plot_style_clean
#' @return A ggplot object with the 'gpt' style applied.
#' @keywords internal
plot_style_gpt <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  colors <- c(
    "NS" = "#999999",               # Pathway gene: Not significant
    "Log2FC" = "#0072B2",           # Pathway gene: Significant Log2FC only
    "p-value" = "#56B4E9",          # Pathway gene: Significant p-value only
    "p-value & Log2FC" = "#E69F00"  # Pathway gene: Significant both
  )

  plot <- base_plot +
    # Non-pathway genes (using alpha)
    geom_point(data = de_results %>% filter(!in_pathway),
               aes(color = color), alpha = 0.3, size = 0.6) +
    # Pathway genes (using alpha)
    geom_point(data = de_results %>% filter(in_pathway),
               aes(color = color), alpha = 0.6, size = 1.5) +
    # Manual color scale
    scale_color_manual(values = colors,
                       name = "Pathway Gene Significance",
                       labels = c("NS" = "NS",
                                  "p-value" = sprintf("p < %.2f", p_cutoff),
                                  "Log2FC" = sprintf("|LogFC| > %.1f", fc_cutoff),
                                  "p-value & Log2FC" = "Both")) +
    # Apply custom theme and adjust legend position
    custom_minimal_theme_with_grid() +
    theme(legend.position = c(0.95, 0.95), # Position legend top-right
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = "white", color = "grey80", size = 0.5)) +
    # Override alpha in the legend guides for clarity
    guides(
      color = guide_legend(override.aes = list(alpha = 1, size = 3)), # Make legend points opaque and larger
      alpha = "none" # Hide the alpha guide if it appears
    )

  # Add labels if needed
  if (nrow(label_data) > 0) {
    plot <- plot + geom_label_repel(
      data = label_data,
      aes(label = rownames(label_data)),
      box.padding = 0.5,
      point.padding = 0.5,
      force = 5,
      segment.color = 'grey50',
      max.overlaps = max_overlaps
    )
  }
  
  return(plot)
}
