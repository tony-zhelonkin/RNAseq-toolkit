library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)  # For separate_rows function

# Main execution
gsea_comparison_plot <- function(x, y, x_label, y_label, 
                                 use_normalized = TRUE, 
                                 percentile_threshold = 0.95,
                                 color = "all",
                                 max_overlaps = 20) { 
  
  # Helper function to get score column name
  get_score_col <- function(data, label) {
    paste0(ifelse(use_normalized, "NES", "enrichmentScore"), ".", label)
  }
  
  # Prepare data
  prepare_data <- function(x, y) {
    score_column <- ifelse(use_normalized, "NES", "enrichmentScore")
    
    if (!all(c(score_column, "ID", "Description", "qvalue", "core_enrichment") %in% names(x)) || 
        !all(c(score_column, "ID", "Description", "qvalue", "core_enrichment") %in% names(y))) {
      stop(paste("Both datasets must contain", score_column, "ID, Description, qvalue, and core_enrichment columns"))
    }
    
    joined_data <- inner_join(x, y, by = "ID", suffix = c(paste0(".", x_label), paste0(".", y_label))) %>%
      select(ID, 
             Description = paste0("Description.", x_label), 
             !!get_score_col(x, x_label), 
             !!get_score_col(y, y_label),
             !!paste0("qvalue.", x_label),
             !!paste0("qvalue.", y_label),
             !!paste0("core_enrichment.", x_label),
             !!paste0("core_enrichment.", y_label))
    
    if (!use_normalized) {
      joined_data <- joined_data %>%
        mutate(across(starts_with("enrichmentScore"), 
                      ~(. - min(.)) / (max(.) - min(.)), 
                      .names = "{.col}"))
    }
    
    joined_data
  }
  
  # Categorize data
  categorize_data <- function(data, x_score_col, y_score_col) {
    x_thresholds <- quantile(data[[x_score_col]], c(1-percentile_threshold, percentile_threshold), na.rm = TRUE)
    y_thresholds <- quantile(data[[y_score_col]], c(1-percentile_threshold, percentile_threshold), na.rm = TRUE)
    
    data %>%
      mutate(
        category = case_when(
          .data[[x_score_col]] > x_thresholds[2] & .data[[y_score_col]] > y_thresholds[2] ~ paste(x_label, "&", y_label, "positive"),
          .data[[x_score_col]] < x_thresholds[1] & .data[[y_score_col]] < y_thresholds[1] ~ paste(x_label, "&", y_label, "negative"),
          .data[[x_score_col]] > x_thresholds[2] ~ paste(x_label, "positive"),
          .data[[x_score_col]] < x_thresholds[1] ~ paste(x_label, "negative"),
          .data[[y_score_col]] > y_thresholds[2] ~ paste(y_label, "positive"),
          .data[[y_score_col]] < y_thresholds[1] ~ paste(y_label, "negative"),
          TRUE ~ "Other"
        )
      )
  }

  # Create plot
  create_plot <- function(data, x_score_col, y_score_col) {
    # Colorblind-friendly palette
    cbf_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#999999")
    names(cbf_palette) <- c(paste(x_label, "&", y_label, "positive"),
                            paste(x_label, "&", y_label, "negative"),
                            paste(x_label, "positive"),
                            paste(x_label, "negative"),
                            paste(y_label, "positive"),
                            paste(y_label, "negative"),
                            "Other")
    
    # Prepare data for plotting based on color argument
    if (color == "common") {
      highlight_categories <- c(paste(x_label, "&", y_label, "positive"), 
                                paste(x_label, "&", y_label, "negative"))
    } else if (color == "distinct") {
      highlight_categories <- c(paste(x_label, "positive"), paste(x_label, "negative"),
                                paste(y_label, "positive"), paste(y_label, "negative"))
    } else {
      highlight_categories <- names(cbf_palette)
    }
    
    data <- data %>%
      mutate(highlight = ifelse(category %in% highlight_categories, category, "Other"))
    
    intercept <- ifelse(use_normalized, 0, 0.5)
    
    p <- ggplot(data, aes(x = .data[[x_score_col]], y = .data[[y_score_col]])) +
      # Plot "Other" points first
      geom_point(data = subset(data, highlight == "Other"), 
                 aes(color = highlight), alpha = 0.3) +
      # Plot highlighted points on top
      geom_point(data = subset(data, highlight != "Other"), 
                 aes(color = highlight), alpha = 0.7) +
      scale_color_manual(values = cbf_palette, 
                         breaks = intersect(names(cbf_palette), unique(data$highlight))) +
      geom_vline(xintercept = intercept, linetype = "dashed", color = "black", alpha = 0.5) +
      geom_hline(yintercept = intercept, linetype = "dashed", color = "black", alpha = 0.5) +
      labs(
        title = paste(ifelse(use_normalized, "NES", "Enrichment Score"), 
                      "Outlier Analysis between", x_label, "and", y_label, "threshold: NES >",  percentile_threshold, "percentile"),
        x = paste(ifelse(use_normalized, "NES", "Enrichment Score"), "(", x_label, ")"),
        y = paste(ifelse(use_normalized, "NES", "Enrichment Score"), "(", y_label, ")"),
        color = "Category"
      ) +
      custom_minimal_theme_with_grid() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),  # White background with black border
        plot.background = element_rect(fill = "white", color = NA),  # White plot background with no border
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.position = "right"
      )
    
    top_pathways <- data %>%
      filter(highlight != "Other") %>%
      group_by(highlight) %>%
      slice_max(order_by = abs(.data[[x_score_col]]) + abs(.data[[y_score_col]]), n = 5) %>%
      ungroup()
    
    p + geom_text_repel(
      data = top_pathways,
      aes(label = Description, color = highlight),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.1,
      segment.color = "grey50",
      segment.alpha = 1,
      min.segment.length = 0,
      force = 1,
      max.overlaps = max_overlaps,
      show.legend = FALSE
    )
  }

  # Create GMT-formatted dataframe
  create_gmt <- function(data, pathway_ids) {
    data %>%
      filter(ID %in% pathway_ids) %>%
      mutate(
        core_enrichment = ifelse(
          !is.na(!!sym(paste0("core_enrichment.", x_label))),
          !!sym(paste0("core_enrichment.", x_label)),
          !!sym(paste0("core_enrichment.", y_label))
        )
      ) %>%
      select(ID, Description, core_enrichment) %>%
      separate_rows(core_enrichment, sep = "/") %>%
      group_by(ID, Description) %>%
      summarise(genes = list(core_enrichment), .groups = "drop") %>%
      unnest(genes) %>%
      group_by(ID, Description) %>%
      mutate(col_num = row_number() + 2) %>%
      pivot_wider(names_from = col_num, values_from = genes, names_prefix = "Gene_") %>%
      ungroup() %>%
      # Ensure the order matches the input pathway_ids
      slice(match(pathway_ids, ID))
  }

  # Main execution
  joined_data <- prepare_data(x, y)
  x_score_col <- get_score_col(x, x_label)
  y_score_col <- get_score_col(y, y_label)
  
  categorized_data <- categorize_data(joined_data, x_score_col, y_score_col)
  
  plot <- create_plot(categorized_data, x_score_col, y_score_col)
  
  # Extract common pathways using the categorized data and order by absolute NES
  common_pos <- categorized_data %>%
    filter(category == paste(x_label, "&", y_label, "positive")) %>%
    mutate(abs_nes = abs(.data[[x_score_col]]) + abs(.data[[y_score_col]])) %>%
    arrange(desc(abs_nes)) %>%
    pull(ID)
  
  common_neg <- categorized_data %>%
    filter(category == paste(x_label, "&", y_label, "negative")) %>%
    mutate(abs_nes = abs(.data[[x_score_col]]) + abs(.data[[y_score_col]])) %>%
    arrange(desc(abs_nes)) %>%
    pull(ID)
  
  common_all <- c(common_pos, common_neg)
  
  common_mix <- bind_rows(
    categorized_data %>% 
      filter(category == paste(x_label, "&", y_label, "positive")) %>% 
      slice_max(order_by = abs(.data[[x_score_col]]) + abs(.data[[y_score_col]]), n = 10),
    categorized_data %>% 
      filter(category == paste(x_label, "&", y_label, "negative")) %>% 
      slice_max(order_by = abs(.data[[x_score_col]]) + abs(.data[[y_score_col]]), n = 10)
  ) %>% 
    arrange(desc(abs(.data[[x_score_col]]) + abs(.data[[y_score_col]]))) %>%
    pull(ID)
  
  # Results DF creation
  results_df <- categorized_data %>%
    filter(category != "Other") %>%
    select(ID, Description, category, 
           !!sym(x_score_col), !!sym(y_score_col),
           !!sym(paste0("qvalue.", x_label)), !!sym(paste0("qvalue.", y_label)),
           !!sym(paste0("core_enrichment.", x_label)), !!sym(paste0("core_enrichment.", y_label))) %>%
    rename(
      !!paste0("NES.", x_label) := !!sym(x_score_col),
      !!paste0("NES.", y_label) := !!sym(y_score_col)
    ) %>%
    arrange(desc(abs(!!sym(paste0("NES.", x_label))) + abs(!!sym(paste0("NES.", y_label)))))
  
  # Create GMT-formatted dataframes for each category
  gmt_pos <- create_gmt(results_df, common_pos)
  gmt_neg <- create_gmt(results_df, common_neg)
  gmt_mix <- create_gmt(results_df, common_mix)
  gmt_all <- create_gmt(results_df, common_all)
  
  # List object 
  list(
    data = results_df, 
    plot = plot, 
    common_pathways = list(
      common_neg = common_neg,
      common_pos = common_pos,
      common_mix = common_mix,
      common_all = common_all
    ),
    gmt = list(
      gmt_neg = gmt_neg,
      gmt_pos = gmt_pos,
      gmt_mix = gmt_mix,
      gmt_all = gmt_all
    )
  )
}
