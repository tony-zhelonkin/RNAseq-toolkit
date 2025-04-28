###############
### Custom GSEA plot function
###############

library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

# Define custom dotplot function for GSEA objects with text processing and sorting options
custom_dotplot <- function(gsea_obj, showCategory = 10, font.size = 10, title = "GSEA Dotplot",
                           replace_ = TRUE, capitalize_1 = TRUE, capitalize_all = FALSE, 
                           filterBy = "qvalue",
                           sortBy = "GeneRatio",
                           q_cut = 0.05,
                           min.dotSize = 2) {
  
  # Extract the result data frame from the GSEA object
  gsea_data <- as.data.frame(gsea_obj@result)
  
  # Calculate the gene count from 'core_enrichment' by counting '/' and adding 1 (number of genes)
  gene_count <- gsea_data %>%
    group_by(ID) %>%
    summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  
  # Merge gene counts with the original GSEA result and calculate GeneRatio
  gsea_data <- left_join(gsea_data, gene_count, by = "ID") %>%
    mutate(GeneRatio = count / setSize)
  
  # Modify Description field based on function arguments
  if (replace_) {
    gsea_data$Description <- gsea_data$Description %>% 
      str_replace_all("_", " ")   # Replace "_" with " " if 'replace' is TRUE
  }
  
  if (capitalize_1) {
    gsea_data$Description <- gsea_data$Description %>%
      str_to_sentence()           # Capitalize only the first word if 'capitalize_1' is TRUE
  }
  
  if (capitalize_all) {
    gsea_data$Description <- gsea_data$Description %>%
      str_to_title()              # Capitalize all words if 'capitalize_all' is TRUE
  }
  
  # Filter for significant pathways (qvalue < 0.01)
  gsea_data_filtered <- gsea_data %>%
    filter(qvalue < q_cut) %>%
    mutate(NES_sign = ifelse(NES > 0, "Positive NES", "Negative NES"))
  
  # Filter logic based on 'filterBy' argument
  if (filterBy == "qvalue") {
    # Sort by qvalue (default behavior)
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(qvalue) %>%           # Sort by qvalue (ascending)
      head(showCategory)
  } else if (filterBy == "NES") {
    # Sort by absolute NES value (strongest enrichment)
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(abs(NES))) %>%    # Sort by absolute NES (descending)
      head(showCategory)
  } else if (filterBy == "NES_positive") {
    # Filter and sort by positive NES values only
    gsea_data_filtered <- gsea_data_filtered %>%
      filter(NES > 0) %>%           # Filter positive NES
      arrange(desc(NES)) %>%        # Sort by NES (descending)
      head(showCategory)
  } else if (filterBy == "NES_negative") {
    # Filter and sort by negative NES values only
    gsea_data_filtered <- gsea_data_filtered %>%
      filter(NES < 0) %>%           # Filter negative NES
      arrange(NES) %>%              # Sort by NES (ascending, more negative first)
      head(showCategory)
  }

 # Sort the filtered data based on the sortBy parameter
  if (sortBy == "GeneRatio") {
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(GeneRatio))
  } else if (sortBy == "qvalue") {
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(qvalue)
  } else {
    warning("Invalid sortBy parameter. Defaulting to GeneRatio.")
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(GeneRatio))
  }

  # Take the top showCategory entries after sorting
  gsea_data_filtered <- gsea_data_filtered %>%
    head(showCategory)

  # Create custom dotplot using ggplot2
  ggplot(gsea_data_filtered, aes(x = GeneRatio, y = reorder(Description, !!sym(sortBy)))) +
    geom_point(aes(size = -log10(qvalue), color = NES_sign)) +
    scale_color_manual(values = c("Positive NES" = "orange", "Negative NES" = "skyblue")) +
    
    # Use scale_size_continuous to set visual size limits for the dots
    scale_size_continuous(range = c(min.dotSize, 10),
                          limits = c(min(-log10(gsea_data_filtered$qvalue)), 
                                     max(-log10(gsea_data_filtered$qvalue))),
                          name = "-log10(qvalue)") +
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      color = "NES",
      size = "-log10(qvalue)"
    ) +
    
    custom_minimal_theme_with_grid() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background with black border
      plot.background = element_rect(fill = "white", color = NA),  # White plot background with no border
      axis.text.y = element_text(size = font.size, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = font.size + 2),
      axis.text.x = element_text(size = font.size),
      legend.position = "right"
    )
}
