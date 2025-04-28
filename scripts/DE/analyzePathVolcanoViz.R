##################
### Combined Volcano Plot Function with Multiple Styles
##################

library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# Custom theme with minimal style and grid (used in multiple styles)
custom_minimal_theme_with_grid <- function() {
  theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey80", size = 0.5),
      panel.grid.minor = element_line(color = "grey90", size = 0.25),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
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

# Style 1: Clean style (original first style)
plot_style_clean <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  # Define the colors for each condition
  colors <- c(
    "NS" = "#999999",               # Not significant
    "p-value & Log2FC" = "#E69F00", # Both p-value and Log2FC significant
    "p-value" = "#56B4E9"           # Only p-value significant
  )
  
  # Filter for pathway-specific data with different highlight criteria
  pathway_data <- de_results %>% filter(in_pathway)
  
  # Base plot with grey points for non-pathway genes
  plot <- base_plot +
    # Non-pathway genes in grey
    geom_point(data = de_results %>% filter(!in_pathway), 
               aes(color = "NS", alpha = 0.6), size = 1) +
    
    # Pathway genes: highlight based on p-value and fold change significance
    geom_point(data = pathway_data %>% filter(highlight), 
               aes(color = "p-value & Log2FC", alpha = 0.7), size = 1.5) +
    
    # Pathway genes: highlight based only on p-value significance
    geom_point(data = pathway_data %>% filter(!highlight & significant_p), 
               aes(color = "p-value", alpha = 0.7), size = 1.5) +
    
    # Set color scale manually and remove legend
    scale_color_manual(values = colors) +
    
    # Apply the custom minimal theme (assuming you have this function defined)
    custom_minimal_theme_with_grid() +
    theme(legend.position = "none")  # Remove the legend entirely
  
  # If there are genes to label, apply label repelling
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


# Style 2: Claude style (second version with enhanced volcano-style legend)
plot_style_claude <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  colors <- c(
    "NS" = "#999999", "p-value" = "#56B4E9", "Log2FC" = "#0072B2", "p-value & Log2FC" = "#E69F00"
  )
  
  plot <- base_plot +
    geom_point(data = de_results %>% filter(!in_pathway), 
               color = "#CCCCCC", alpha = 0.4, size = 0.6) +
    geom_point(data = de_results %>% filter(in_pathway), 
               aes(color = color), alpha = 0.7, size = 1.5) +
    scale_color_manual(values = colors, 
                       name = "Significance",
                       labels = c("NS" = "NS", "p-value" = "p-value",
                                  "Log2FC" = "Log2FC", "p-value & Log2FC" = "p-value & Log2FC")) +
    custom_minimal_theme_with_grid() +
    theme(
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = "white", color = "black")
    )
  
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

# Style 3: GPT style (third version with no alpha in legend)
plot_style_gpt <- function(base_plot, de_results, pathway_genes, label_data, p_cutoff, fc_cutoff, max_overlaps) {
  colors <- c("NS" = "#999999", "Log2FC" = "#0072B2", "p-value" = "#56B4E9", "p-value & Log2FC" = "#E69F00")
  
  plot <- base_plot +
    geom_point(data = de_results %>% filter(!in_pathway), 
               aes(color = color, alpha = 0.3), size = 0.6) +
    geom_point(data = de_results %>% filter(in_pathway), 
               aes(color = color, alpha = 0.6), size = 1.5) +
    scale_color_manual(values = colors) +
    custom_minimal_theme_with_grid()+
    theme(legend.position = c(0.9, 0.9)) +
    guides(
      color = guide_legend(title = "Gene Status", override.aes = list(alpha = 1, size = 4)),
      alpha = "none"  # Remove alpha from legend
    )
  
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
