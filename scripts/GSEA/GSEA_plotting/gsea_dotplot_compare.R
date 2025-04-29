#' Create a Comparative GSEA Dotplot for Two Datasets
#'
#' This function creates a side-by-side dotplot comparison of GSEA results from two samples,
#' allowing for easy visualization of enriched pathways across both datasets.
#'
#' @param gsea_obj_x A GSEA result object from clusterProfiler for the first sample.
#' @param gsea_obj_y A GSEA result object from clusterProfiler for the second sample.
#' @param pathway_ids Character vector of pathway IDs to include in the plot.
#' @param font.size Numeric, base font size for the plot (default: 10).
#' @param title Character, plot title (default: "GSEA Comparison Dotplot").
#' @param replace_ Logical, whether to replace underscores with spaces in descriptions (default: TRUE).
#' @param capitalize_1 Logical, whether to capitalize the first word in descriptions (default: FALSE).
#' @param capitalize_all Logical, whether to capitalize all words in descriptions (default: FALSE).
#' @param min.dotSize Numeric, minimum dot size in the plot (default: 2).
#' @param sample_x_name Character, name for the first sample (default: "Sample X").
#' @param sample_y_name Character, name for the second sample (default: "Sample Y").
#' @param sortBy Character, method to sort pathways: "importance_score", "qvalue", "x", or "y" (default: "importance_score").
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 12).
#' @param height Numeric, height of the saved plot in inches (default: 8).
#' @param dpi Numeric, resolution of the saved plot (default: 300).
#'
#' @return A ggplot2 object representing the comparative GSEA dotplot.
#' @export
#'
#' @examples
#' # Basic usage with pathway IDs
#' gsea_dotplot_compare(gsea_obj_x, gsea_obj_y, pathway_ids)
#'
#' # Customized sample names
#' gsea_dotplot_compare(gsea_obj_x, gsea_obj_y, pathway_ids,
#'                      sample_x_name = "Treatment", sample_y_name = "Control")
#'
#' # Save the plot
#' gsea_dotplot_compare(gsea_obj_x, gsea_obj_y, pathway_ids,
#'                      title = "Treatment vs Control", save_plot = TRUE)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

gsea_dotplot_compare <- function(gsea_obj_x, 
                                 gsea_obj_y, 
                                 pathway_ids, 
                                 font.size = 10, 
                                 title = "GSEA Comparison Dotplot",
                                 replace_ = TRUE, 
                                 capitalize_1 = FALSE, 
                                 capitalize_all = FALSE,
                                 min.dotSize = 2, 
                                 sample_x_name = "Sample X", 
                                 sample_y_name = "Sample Y",
                                 sortBy = "importance_score",
                                 save_plot = FALSE,
                                 output_dir = "plots/",
                                 width = 12,
                                 height = 8,
                                 dpi = 300) {
  
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
      arrange(desc(importance_score)) %>%
      pull(Description),
    "qvalue" = pathway_scores %>%
      arrange(mean_qvalue) %>%
      pull(Description),
    "x" = data_x %>%
      arrange(desc(GeneRatio)) %>%
      pull(Description),
    "y" = data_y %>%
      arrange(desc(GeneRatio)) %>%
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
  p <- ggplot(combined_data, aes(x = GeneRatio, y = Description)) +
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
    custom_minimal_theme_with_grid() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.y = element_text(size = font.size, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = font.size + 2),
      axis.text.x = element_text(size = font.size),
      legend.position = "right",
      strip.text = element_text(size = font.size + 1),
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  # Save the plot if requested
  if (save_plot) {
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create filename from title
    filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_compare.pdf"))
    
    ggsave(
      filename = filename,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
    
    message("Plot saved to: ", filename)
  }
  
  return(p)
}
