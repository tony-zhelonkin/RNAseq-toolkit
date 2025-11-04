library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

custom_dotplot_comparison <- function(gsea_obj_x, gsea_obj_y, pathway_ids, 
                                      font.size = 10, title = "Combined GSEA Dotplot",
                                      replace_ = TRUE, capitalize_1 = FALSE, capitalize_all = FALSE,
                                      min.dotSize = 2, sample_x_name = "Sample X", sample_y_name = "Sample Y",
                                      sortBy = "importance_score") {
  
  # Helper function to process GSEA data
  process_gsea_data <- function(gsea_obj, sample_name) {
    gsea_data <- as.data.frame(gsea_obj@result) %>%
      mutate(
        count = str_count(core_enrichment, "/") + 1,
        GeneRatio = count / setSize,
        Description = if (replace_) str_replace_all(Description, "_", " ") else Description,
        Description = case_when(
          capitalize_all ~ str_to_title(Description),
          capitalize_1 ~ str_to_sentence(Description),
          TRUE ~ Description
        ),
        NES_sign = ifelse(NES > 0, "Positive NES", "Negative NES"),
        Sample = sample_name
      ) %>%
      select(ID, Description, NES, qvalue, GeneRatio, NES_sign, Sample)
  }
 
  # Process both GSEA objects and combine data
  data_x <- process_gsea_data(gsea_obj_x, sample_x_name)
  data_y <- process_gsea_data(gsea_obj_y, sample_y_name)
  combined_data <- bind_rows(data_x, data_y) %>%
    filter(ID %in% pathway_ids)
  
  # Calculate importance score for sorting
  pathway_scores <- combined_data %>%
    group_by(ID, Description) %>%
    summarize(
      mean_GeneRatio = mean(GeneRatio),
      mean_qvalue = mean(qvalue),
      importance_score = mean_GeneRatio * -log10(mean_qvalue),
      .groups = "drop"
    )
  
  # Determine the sorting order based on the sortBy parameter
  sort_order <- switch(sortBy,
    "importance_score" = pathway_scores %>%
      arrange(importance_score) %>%
      pull(Description),
    "qvalue" = pathway_scores %>%
      arrange(mean_qvalue) %>%
      pull(Description),
    "x" = data_x %>%
      arrange(GeneRatio) %>%
      pull(Description),
    "y" = data_y %>%
      arrange(GeneRatio) %>%
      pull(Description),
    {
      warning("Invalid sortBy parameter. Defaulting to importance_score.")
      pathway_scores %>%
        arrange(desc(importance_score)) %>%
        pull(Description)
    }
  )

  # Apply the sorting order and handle missing pathways
  combined_data <- combined_data %>%
    mutate(Description = factor(Description, levels = sort_order))
  
  # Check for missing pathways using ID instead of Description
  missing_pathways <- setdiff(pathway_ids, unique(combined_data$ID))
  if (length(missing_pathways) > 0) {
    warning("Some pathway IDs are missing from the dataset: ", paste(missing_pathways, collapse = ", "))
  }
  
  # Create the plot
  ggplot(combined_data, aes(x = GeneRatio, y = Description)) +
    geom_point(aes(size = -log10(qvalue), color = NES_sign)) +
    scale_color_manual(values = c("Positive NES" = "orange", "Negative NES" = "skyblue")) +
    scale_size_continuous(range = c(min.dotSize, 10),
                          limits = c(min(-log10(combined_data$qvalue)), 
                                     max(-log10(combined_data$qvalue))),
                          name = "-log10(qvalue)") +
    facet_wrap(~ Sample, scales = "fixed", ncol = 2) +
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      color = "NES"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background with black border
      plot.background = element_rect(fill = "white", color = NA),  # White plot background with no border
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text.y = element_text(size = font.size, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = font.size + 2),
      axis.text.x = element_text(size = font.size),
      legend.position = "right",
      strip.text = element_text(size = font.size + 1),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25),
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
}
