#' Create a Comparative GSEA Dotplot for Two Datasets
#'
#' Creates a side-by-side dotplot comparison of GSEA results from two samples/conditions,
#' allowing for easy visualization of commonly enriched pathways. Uses `p.adjust` for dot size.
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param gsea_obj_x A `gseaResult` object from `clusterProfiler` for the first sample/condition.
#' @param gsea_obj_y A `gseaResult` object from `clusterProfiler` for the second sample/condition.
#' @param pathway_ids Character vector of pathway IDs (`ID` column in GSEA results) to include.
#' @param base_font_size Numeric, base font size for plot text elements (default: 10).
#' @param title Character, plot title (default: "GSEA Comparison Dotplot").
#' @param replace_ Logical, if TRUE (default), replace underscores "_" with spaces " " in descriptions.
#' @param capitalize_1 Logical, if TRUE, capitalize the first letter of descriptions (default: FALSE).
#' @param capitalize_all Logical, if TRUE, capitalize the first letter of each word (default: FALSE).
#' @param min.dotSize Numeric, minimum size for the dots (default: 2).
#' @param max.dotSize Numeric, maximum size for the dots (default: 10).
#' @param sample_x_name Character, name for the first sample/condition (facet label) (default: "Sample X").
#' @param sample_y_name Character, name for the second sample/condition (facet label) (default: "Sample Y").
#' @param sortBy Character, method to sort pathways: "importance_score" (default), "p.adjust", "x", or "y".
#'        "p.adjust" sorts by mean p.adjust. "x" sorts by GeneRatio in sample X, "y" by GeneRatio in sample Y.
#' @param pos_color Character, color for positive NES dots (default: "orange").
#' @param neg_color Character, color for negative NES dots (default: "skyblue").
#'
#' @return A ggplot2 object representing the comparative GSEA dotplot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2 tidyr
#' @importFrom dplyr %>% filter select mutate bind_rows group_by summarize arrange desc pull left_join case_when
#' @importFrom stringr str_count str_replace_all str_to_sentence str_to_title
#' @importFrom methods is slot
#' @importFrom stats reorder setNames
#' @importFrom rlang .data
#'
#' @examples
#' # Assuming gsea_res_x and gsea_res_y are valid gseaResult objects
#' # and common_path_ids is a vector of pathway IDs present in both
#' # p <- gsea_dotplot_compare(gsea_res_x, gsea_res_y, common_path_ids,
#' #                           sample_x_name = "Treatment", sample_y_name = "Control")
#' # print(p)
#' # To save: ggsave("my_gsea_compare_dotplot.png", p)

# Function to safely source a script if it exists
source_safe <- function(path) {
  if (file.exists(path)) {
    source(path)
    return(TRUE)
  } else {
    warning("Custom theme script not found: ", path, ". Using default theme_minimal().")
    custom_minimal_theme_with_grid <<- function() theme_minimal() # Define placeholder
    return(FALSE)
  }
}
# Source the custom theme
source_safe("scripts/custom_minimal_theme.R")


gsea_dotplot_compare <- function(gsea_obj_x,
                                 gsea_obj_y,
                                 pathway_ids,
                                 base_font_size = 10, # Renamed
                                 title = "GSEA Comparison Dotplot",
                                 replace_ = TRUE,
                                 capitalize_1 = FALSE,
                                 capitalize_all = FALSE,
                                 min.dotSize = 2,
                                 max.dotSize = 10, # Added max size parameter
                                 sample_x_name = "Sample X",
                                 sample_y_name = "Sample Y",
                                 sortBy = "importance_score",
                                 pos_color = "orange", # Added color param
                                 neg_color = "skyblue") { # Added color param

  # --- Input Validation ---
  validate_gsea_obj <- function(obj, name) {
      if (!methods::is(obj, "gseaResult")) {
          stop(sprintf("Input `%s` must be a gseaResult object.", name))
      }
      if (!methods::.hasSlot(obj, "result") || !is.data.frame(obj@result) || nrow(obj@result) == 0) {
          stop(sprintf("Input `%s` has an invalid or empty result slot.", name))
      }
      required_cols <- c("ID", "Description", "NES", "p.adjust", "core_enrichment", "setSize")
      if (!all(required_cols %in% colnames(obj@result))) {
          stop(sprintf("Input `%s@result` is missing required columns: %s", name,
                       paste(setdiff(required_cols, colnames(obj@result)), collapse=", ")))
      }
  }
  validate_gsea_obj(gsea_obj_x, "gsea_obj_x")
  validate_gsea_obj(gsea_obj_y, "gsea_obj_y")
  if (!is.character(pathway_ids) || length(pathway_ids) == 0) {
      stop("`pathway_ids` must be a non-empty character vector.")
  }
  # ------------------------

  # Helper function to process GSEA data
  process_gsea_data <- function(gsea_obj, sample_name) {
      as.data.frame(gsea_obj@result) %>%
          dplyr::mutate(
              # Calculate gene count from core_enrichment string
              count = stringr::str_count(.data$core_enrichment, "/") + ifelse(nchar(.data$core_enrichment) > 0, 1, 0),
              GeneRatio = .data$count / .data$setSize,
              Description = if (replace_) stringr::str_replace_all(.data$Description, "_", " ") else .data$Description,
              Description = dplyr::case_when(
                  capitalize_all ~ stringr::str_to_title(.data$Description),
                  capitalize_1 ~ stringr::str_to_sentence(.data$Description),
                  TRUE ~ .data$Description
              ),
              NES_sign = ifelse(.data$NES > 0, "Positive NES", "Negative NES"),
              Sample = sample_name,
              # Use p.adjust for significance measure
              negLog10pAdj = -log10(.data$p.adjust)
          ) %>%
          # Select relevant columns, using p.adjust
          dplyr::select(.data$ID, .data$Description, .data$NES, .data$p.adjust, .data$GeneRatio,
                        .data$negLog10pAdj, .data$NES_sign, .data$Sample)
  }

  # Process both GSEA objects
  data_x <- process_gsea_data(gsea_obj_x, sample_x_name)
  data_y <- process_gsea_data(gsea_obj_y, sample_y_name)

  # Combine data and filter for selected pathway IDs
  combined_data <- dplyr::bind_rows(data_x, data_y) %>%
      dplyr::filter(.data$ID %in% pathway_ids)

  # Check if any data remains after filtering
  if (nrow(combined_data) == 0) {
      warning("No data found for the specified pathway IDs. Returning an empty plot.")
      return(ggplot() + labs(title = paste(title, "(No matching pathways)")))
  }

  # Calculate importance score for sorting (using mean p.adjust)
  pathway_scores <- combined_data %>%
      dplyr::group_by(.data$ID, .data$Description) %>%
      dplyr::summarize(
          mean_GeneRatio = mean(.data$GeneRatio, na.rm = TRUE),
          mean_p.adjust = mean(.data$p.adjust, na.rm = TRUE),
          # Calculate score using p.adjust, handle p.adjust=0 case
          importance_score = mean_GeneRatio * (-log10(pmax(.data$mean_p.adjust, .Machine$double.eps))), # Use pmax for stability
          .groups = "drop"
      )

  # Determine the sorting order based on the sortBy parameter
  sort_order <- switch(sortBy,
      "importance_score" = pathway_scores %>% dplyr::arrange(dplyr::desc(.data$importance_score)) %>% dplyr::pull(.data$Description),
      "p.adjust" = pathway_scores %>% dplyr::arrange(.data$mean_p.adjust) %>% dplyr::pull(.data$Description),
      "x" = data_x %>% dplyr::filter(.data$ID %in% pathway_ids) %>% dplyr::arrange(dplyr::desc(.data$GeneRatio)) %>% dplyr::pull(.data$Description),
      "y" = data_y %>% dplyr::filter(.data$ID %in% pathway_ids) %>% dplyr::arrange(dplyr::desc(.data$GeneRatio)) %>% dplyr::pull(.data$Description),
      {
          warning("Invalid sortBy parameter '", sortBy, "'. Defaulting to 'importance_score'.")
          pathway_scores %>% dplyr::arrange(dplyr::desc(.data$importance_score)) %>% dplyr::pull(.data$Description)
      }
  )
  # Ensure unique descriptions in sort order (in case of ties in sorting variable)
  sort_order <- unique(sort_order)

  # Apply the sorting order as factor levels
  combined_data <- combined_data %>%
      dplyr::mutate(Description = factor(.data$Description, levels = rev(sort_order))) # Reverse for ggplot top-to-bottom

  # Check for pathways requested but not found in the data
  missing_pathways <- setdiff(pathway_ids, unique(combined_data$ID))
  if (length(missing_pathways) > 0) {
      warning("Some requested pathway IDs were not found in the provided GSEA results: ",
              paste(missing_pathways, collapse = ", "))
  }

  # Handle cases where p.adjust is 0 or very small for scaling
  combined_data <- combined_data %>%
      dplyr::mutate(negLog10pAdj_capped = pmax(.data$negLog10pAdj, 0)) # Ensure non-negative size

  # Determine size limits, handling potential Inf values
  size_limits <- range(combined_data$negLog10pAdj_capped[is.finite(combined_data$negLog10pAdj_capped)], na.rm = TRUE)
  if (all(!is.finite(size_limits))) size_limits <- c(0, 1) # Default if all are Inf/NA

  # Create the plot
  p <- ggplot(combined_data, aes(x = .data$GeneRatio, y = .data$Description)) +
      geom_point(aes(size = .data$negLog10pAdj_capped, color = .data$NES_sign)) +
      scale_color_manual(name = "NES Sign", values = c("Positive NES" = pos_color, "Negative NES" = neg_color)) +
      scale_size_continuous(name = bquote(-log[10](p.adjust)), # Use p.adjust in label
                            range = c(min.dotSize, max.dotSize),
                            limits = size_limits) +
      facet_wrap(~ Sample, scales = "fixed", ncol = 2) +
      labs(
          title = title,
          x = "Gene Ratio (Genes in Pathway / Set Size)",
          y = NULL # No y-axis label needed
      ) +
      custom_minimal_theme_with_grid() +
      theme(
          axis.text.y = element_text(size = rel(0.9) * base_font_size, hjust = 1), # Relative sizing
          plot.title = element_text(hjust = 0.5, size = rel(1.1) * base_font_size),
          axis.title.x = element_text(size = rel(1) * base_font_size),
          axis.text.x = element_text(size = rel(0.9) * base_font_size),
          legend.title = element_text(size = rel(0.9) * base_font_size),
          legend.text = element_text(size = rel(0.8) * base_font_size),
          legend.position = "right",
          strip.text = element_text(size = rel(1) * base_font_size, face = "bold"),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5) # Match theme
      )

  # Return the plot object
  return(p)
}
