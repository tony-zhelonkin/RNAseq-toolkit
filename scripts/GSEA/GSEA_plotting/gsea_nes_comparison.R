#' Compare GSEA Results Between Two Datasets
#'
#' Creates a scatter plot comparing Enrichment Scores (ES) or Normalized Enrichment
#' Scores (NES) from two GSEA result objects (`gseaResult` from `clusterProfiler`).
#' Categorizes pathways based on score thresholds (defined by `percentile_threshold`)
#' and highlights selected categories. Also extracts lists of common pathway IDs
#' and generates GMT-like dataframes.
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param gsea_obj_x A `gseaResult` object for the first dataset (x-axis).
#' @param gsea_obj_y A `gseaResult` object for the second dataset (y-axis).
#' @param x_label Character, label for the first dataset (default: "Dataset 1").
#' @param y_label Character, label for the second dataset (default: "Dataset 2").
#' @param use_normalized Logical, if TRUE (default), use NES. If FALSE, use raw Enrichment Score.
#' @param percentile_threshold Numeric, percentile used to define score thresholds for
#'        categorization (e.g., 0.95 means top 5% and bottom 5%). Default: 0.95.
#' @param color Character, specifies which categories to highlight distinctly:
#'        "all" (default, all categories), "common" (only common positive/negative),
#'        "distinct" (only dataset-specific positive/negative).
#' @param base_font_size Numeric, base font size for plot text (default: 12).
#' @param max_overlaps Integer, maximum label overlaps allowed by `ggrepel` (default: 20).
#'
#' @return A list containing:
#'   \item{data}{Data frame with combined and categorized results (uses `p.adjust`).}
#'   \item{plot}{The ggplot object for the comparison scatter plot.}
#'   \item{common_pathways}{List of pathway IDs common to both datasets based on categories
#'          (common_pos, common_neg, common_mix, common_all).}
#'   \item{gmt}{List of GMT-like dataframes for common pathway categories.}
#' @export
#' @import ggplot2 ggrepel tidyr
#' @importFrom dplyr %>% filter select mutate bind_rows group_by summarize arrange desc pull left_join case_when slice_max ungroup rename starts_with across sym slice
#' @importFrom stringr str_replace_all
#' @importFrom methods is slot
#' @importFrom stats quantile setNames
#' @importFrom rlang .data :=
#' @importFrom utils head
#'
#' @examples
#' # Assuming gsea_res_x and gsea_res_y are valid gseaResult objects
#' # comparison_results <- gsea_nes_comparison(gsea_res_x, gsea_res_y,
#' #                                           x_label = "Treatment", y_label = "Control")
#' # print(comparison_results$plot)
#' # common_ids <- comparison_results$common_pathways$common_mix

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
gsea_nes_comparison <- function(gsea_obj_x,
                               gsea_obj_y,
                               x_label = "Dataset 1",
                               y_label = "Dataset 2",
                               use_normalized = TRUE,
                               percentile_threshold = 0.95,
                               color = "all", # "all", "common", "distinct"
                               base_font_size = 12, # Added
                               max_overlaps = 20) {

  # --- Input Validation ---
  validate_gsea_obj <- function(obj, name) {
      if (!methods::is(obj, "gseaResult")) stop(sprintf("Input `%s` must be a gseaResult object.", name))
      if (!methods::.hasSlot(obj, "result") || !is.data.frame(obj@result) || nrow(obj@result) == 0) stop(sprintf("Input `%s` has an invalid or empty result slot.", name))
      score_col <- ifelse(use_normalized, "NES", "enrichmentScore")
      required <- c("ID", "Description", score_col, "p.adjust", "core_enrichment")
      if (!all(required %in% colnames(obj@result))) stop(sprintf("Input `%s@result` missing required columns: %s", name, paste(setdiff(required, colnames(obj@result)), collapse=", ")))
  }
  validate_gsea_obj(gsea_obj_x, "gsea_obj_x")
  validate_gsea_obj(gsea_obj_y, "gsea_obj_y")
  if (!color %in% c("all", "common", "distinct")) {
      warning("Invalid `color` argument. Defaulting to 'all'.")
      color <- "all"
  }
  # ------------------------


  # Helper function to get score column name based on use_normalized
  get_score_col_name <- function(label) {
    paste0(ifelse(use_normalized, "NES", "enrichmentScore"), ".", label)
  }
  score_col_x_base <- ifelse(use_normalized, "NES", "enrichmentScore")
  score_col_y_base <- ifelse(use_normalized, "NES", "enrichmentScore")
  score_col_x_label <- get_score_col_name(x_label)
  score_col_y_label <- get_score_col_name(y_label)


  # Prepare data: join results, select columns, rename using labels
  prepare_data <- function(gsea_x, gsea_y) {
      df_x <- as.data.frame(gsea_x@result)
      df_y <- as.data.frame(gsea_y@result)

      # Select and rename columns before joining
      df_x_sel <- df_x %>%
          dplyr::select(ID, Description,
                        !!score_col_x_label := !!sym(score_col_x_base),
                        !!paste0("p.adjust.", x_label) := p.adjust,
                        !!paste0("core_enrichment.", x_label) := core_enrichment)

      df_y_sel <- df_y %>%
          dplyr::select(ID, # Description might differ, join only on ID
                        !!score_col_y_label := !!sym(score_col_y_base),
                        !!paste0("p.adjust.", y_label) := p.adjust,
                        !!paste0("core_enrichment.", y_label) := core_enrichment)

      # Inner join, keeping description from X
      dplyr::inner_join(df_x_sel, df_y_sel, by = "ID")
  }


  # Categorize data based on score thresholds
  categorize_data <- function(data) {
      # Calculate thresholds based on percentiles
      x_thresholds <- stats::quantile(data[[score_col_x_label]], c(1 - percentile_threshold, percentile_threshold), na.rm = TRUE)
      y_thresholds <- stats::quantile(data[[score_col_y_label]], c(1 - percentile_threshold, percentile_threshold), na.rm = TRUE)

      data %>%
          dplyr::mutate(
              category = dplyr::case_when(
                  .data[[score_col_x_label]] > x_thresholds[2] & .data[[score_col_y_label]] > y_thresholds[2] ~ paste(x_label, "&", y_label, "positive"),
                  .data[[score_col_x_label]] < x_thresholds[1] & .data[[score_col_y_label]] < y_thresholds[1] ~ paste(x_label, "&", y_label, "negative"),
                  .data[[score_col_x_label]] > x_thresholds[2] ~ paste(x_label, "positive"),
                  .data[[score_col_x_label]] < x_thresholds[1] ~ paste(x_label, "negative"),
                  .data[[score_col_y_label]] > y_thresholds[2] ~ paste(y_label, "positive"),
                  .data[[score_col_y_label]] < y_thresholds[1] ~ paste(y_label, "negative"),
                  TRUE ~ "Other"
              )
          )
  }


  # Create the scatter plot
  create_plot <- function(data) {
      # Define color palette
      cbf_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999") # Colorblind friendly
      category_levels <- c(paste(x_label, "&", y_label, "positive"),
                           paste(x_label, "&", y_label, "negative"),
                           paste(x_label, "positive"), paste(x_label, "negative"),
                           paste(y_label, "positive"), paste(y_label, "negative"),
                           "Other")
      names(cbf_palette) <- category_levels

      # Determine which categories to highlight based on 'color' parameter
      if (color == "common") {
          highlight_cats <- category_levels[1:2]
      } else if (color == "distinct") {
          highlight_cats <- category_levels[3:6]
      } else { # color == "all"
          highlight_cats <- category_levels[1:6] # Highlight all except "Other"
      }

      # Add a 'highlight' column for easier plotting/coloring
      plot_data <- data %>%
          dplyr::mutate(highlight = ifelse(.data$category %in% highlight_cats, .data$category, "Other")) %>%
          # Ensure 'Other' is plotted first (underneath)
          dplyr::arrange(.data$highlight == "Other")

      # Determine intercept based on score type
      intercept_val <- ifelse(use_normalized, 0, mean(c(data[[score_col_x_label]], data[[score_col_y_label]]), na.rm=TRUE)) # Mean for raw ES? Or 0? Let's use 0 for simplicity.
      intercept_val <- 0 # Override: Use 0 for both NES and raw ES comparison plots

      p <- ggplot(plot_data, aes(x = .data[[score_col_x_label]], y = .data[[score_col_y_label]])) +
          geom_point(aes(color = .data$highlight), alpha = 0.6, size = 2.5) +
          scale_color_manual(name = "Category", values = cbf_palette,
                             breaks = category_levels, # Ensure all levels are in legend potentially
                             labels = stringr::str_replace_all(category_levels, "&", " & ")) + # Add spaces
          geom_vline(xintercept = intercept_val, linetype = "dashed", color = "grey50") +
          geom_hline(yintercept = intercept_val, linetype = "dashed", color = "grey50") +
          labs(
              title = paste(ifelse(use_normalized, "NES", "Enrichment Score"), "Comparison"),
              subtitle = paste(x_label, "vs", y_label, "| Threshold:", percentile_threshold*100, "th percentile"),
              x = paste(ifelse(use_normalized, "NES", "ES"), "-", x_label),
              y = paste(ifelse(use_normalized, "NES", "ES"), "-", y_label)
          ) +
          custom_minimal_theme_with_grid() +
          theme(
              plot.title = element_text(hjust = 0.5, size = rel(1.1) * base_font_size),
              plot.subtitle = element_text(hjust = 0.5, size = rel(0.9) * base_font_size),
              axis.text = element_text(size = rel(0.9) * base_font_size),
              axis.title = element_text(size = rel(1) * base_font_size),
              legend.title = element_text(size = rel(0.9) * base_font_size),
              legend.text = element_text(size = rel(0.8) * base_font_size),
              legend.position = "right"
          )

      # Add labels for top pathways in highlighted categories
      label_data <- plot_data %>%
          dplyr::filter(.data$highlight != "Other") %>%
          dplyr::group_by(.data$highlight) %>%
          # Select top pathways based on combined distance from origin (or intercept)
          dplyr::slice_max(order_by = abs(.data[[score_col_x_label]] - intercept_val) + abs(.data[[score_col_y_label]] - intercept_val), n = 5) %>%
          dplyr::ungroup()

      p <- p + ggrepel::geom_text_repel(
          data = label_data,
          aes(label = .data$Description, color = .data$highlight),
          size = rel(2.5 / .pt), # Adjust size relative to theme
          box.padding = 0.4,
          point.padding = 0.2,
          segment.color = "grey50",
          segment.alpha = 0.8,
          min.segment.length = 0.1,
          force = 1.5,
          max.overlaps = max_overlaps,
          show.legend = FALSE # Avoid duplicating legend
      )

      return(p)
  }


  # Create GMT-like dataframe (simplified)
  create_gmt <- function(data, pathway_ids_to_include) {
      # Ensure data only contains the pathways we need
      filtered_data <- data %>% dplyr::filter(.data$ID %in% pathway_ids_to_include)

      if (nrow(filtered_data) == 0) return(data.frame(ID=character(), Description=character())) # Return empty df if no pathways

      # Prioritize core enrichment from X, fallback to Y
      core_enrich_col_x <- paste0("core_enrichment.", x_label)
      core_enrich_col_y <- paste0("core_enrichment.", y_label)

      gmt_ready <- filtered_data %>%
          dplyr::mutate(
              # Choose the first non-NA core enrichment string
              core_enrichment = dplyr::case_when(
                  !is.na(.data[[core_enrich_col_x]]) & nchar(.data[[core_enrich_col_x]]) > 0 ~ .data[[core_enrich_col_x]],
                  !is.na(.data[[core_enrich_col_y]]) & nchar(.data[[core_enrich_col_y]]) > 0 ~ .data[[core_enrich_col_y]],
                  TRUE ~ "" # Should not happen if input validation passed
              )
          ) %>%
          dplyr::select(.data$ID, .data$Description, .data$core_enrichment) %>%
          # Split core enrichment string into genes
          dplyr::mutate(genes = strsplit(.data$core_enrichment, "/")) %>%
          dplyr::select(-.data$core_enrichment) %>%
          tidyr::unnest(.data$genes) %>%
          # Remove empty strings if any resulted from splitting
          dplyr::filter(nzchar(.data$genes)) %>%
          # Group to prepare for wide format (basic GMT structure)
          dplyr::group_by(.data$ID, .data$Description) %>%
          # Create a simple list column first
          dplyr::summarise(gene_list = list(unique(.data$genes)), .groups = "drop") %>%
          # Order according to input pathway IDs
           dplyr::slice(match(pathway_ids_to_include, .data$ID))

      # Note: This doesn't create the exact tab-separated wide GMT format,
      # but provides the necessary info (ID, Description, list of genes).
      # Further processing would be needed for true GMT file writing.
      return(gmt_ready)
  }


  # --- Main Execution ---
  joined_data <- prepare_data(gsea_obj_x, gsea_obj_y)

  # Check if join resulted in data
  if (nrow(joined_data) == 0) {
      warning("No common pathway IDs found between the two datasets.")
      return(list(data = data.frame(), plot = ggplot() + labs(title="No Common Pathways"), common_pathways = list(), gmt = list()))
  }

  categorized_data <- categorize_data(joined_data)

  plot <- create_plot(categorized_data)

  # Extract common pathway IDs based on categories
  extract_ids <- function(category_name) {
      categorized_data %>%
          dplyr::filter(.data$category == category_name) %>%
          # Order by combined distance from origin for consistency
          dplyr::mutate(dist_sq = .data[[score_col_x_label]]^2 + .data[[score_col_y_label]]^2) %>%
          dplyr::arrange(dplyr::desc(.data$dist_sq)) %>%
          dplyr::pull(.data$ID)
  }

  common_pos_ids <- extract_ids(paste(x_label, "&", y_label, "positive"))
  common_neg_ids <- extract_ids(paste(x_label, "&", y_label, "negative"))
  common_all_ids <- c(common_pos_ids, common_neg_ids)

  # For common_mix, take top N from positive and negative common lists
  n_mix <- 10 # Number to take from each category for the 'mix' list
  common_mix_ids <- c(utils::head(common_pos_ids, n_mix), utils::head(common_neg_ids, n_mix))
  # Reorder common_mix based on original categorized data's distance metric
  common_mix_ids <- categorized_data %>%
      dplyr::filter(.data$ID %in% common_mix_ids) %>%
      dplyr::mutate(dist_sq = .data[[score_col_x_label]]^2 + .data[[score_col_y_label]]^2) %>%
      dplyr::arrange(dplyr::desc(.data$dist_sq)) %>%
      dplyr::pull(.data$ID)


  # Prepare final results dataframe (using p.adjust)
  results_df <- categorized_data %>%
      # Use p.adjust columns
      dplyr::select(.data$ID, .data$Description, .data$category,
                    !!score_col_x_label, !!score_col_y_label,
                    !!paste0("p.adjust.", x_label), !!paste0("p.adjust.", y_label),
                    !!paste0("core_enrichment.", x_label), !!paste0("core_enrichment.", y_label)) %>%
      # Rename score columns for clarity if needed (already done in prepare_data)
      # dplyr::rename(
      #     !!paste0(score_col_x_base, ".", x_label) := !!sym(score_col_x_label),
      #     !!paste0(score_col_y_base, ".", y_label) := !!sym(score_col_y_label)
      # ) %>%
      # Order by combined distance from origin
      dplyr::mutate(dist_sq = .data[[score_col_x_label]]^2 + .data[[score_col_y_label]]^2) %>%
      dplyr::arrange(dplyr::desc(.data$dist_sq)) %>%
      dplyr::select(-.data$dist_sq) # Remove helper column


  # Create GMT-like dataframes
  gmt_pos <- create_gmt(results_df, common_pos_ids)
  gmt_neg <- create_gmt(results_df, common_neg_ids)
  gmt_mix <- create_gmt(results_df, common_mix_ids)
  gmt_all <- create_gmt(results_df, common_all_ids)


  # Return list object
  list(
      data = results_df,
      plot = plot,
      common_pathways = list(
          common_neg = common_neg_ids,
          common_pos = common_pos_ids,
          common_mix = common_mix_ids, # Top N positive & negative common
          common_all = common_all_ids  # All positive & negative common
      ),
      gmt = list(
          gmt_neg = gmt_neg,
          gmt_pos = gmt_pos,
          gmt_mix = gmt_mix,
          gmt_all = gmt_all
      )
  )
}
