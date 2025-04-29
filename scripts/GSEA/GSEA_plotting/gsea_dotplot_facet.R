#' Create Faceted GSEA Dotplot (Up/Down Regulated)
#'
#' Creates a dotplot of GSEA results, faceted by enrichment direction (positive/negative NES).
#' Displays the top `showCategory` pathways for each direction based on NES, filtered by `padj_cutoff`.
#' Dot size represents significance (`-log10(p.adjust)`).
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param showCategory Integer, the number of top categories to display in each facet
#'        (ranked by absolute NES after filtering) (default: 10).
#' @param base_font_size Numeric, base font size for plot text elements (default: 10).
#' @param title Character, plot title (default: "GSEA Faceted Dotplot").
#' @param padj_cutoff Numeric, adjusted p-value cutoff for filtering significant pathways
#'        (default: 0.05). Uses the `p.adjust` column.
#' @param replace_ Logical, if TRUE (default), replace underscores "_" with spaces " " in descriptions.
#' @param capitalize_1 Logical, if TRUE, capitalize the first letter of descriptions (default: FALSE).
#' @param capitalize_all Logical, if TRUE, capitalize the first letter of each word (default: FALSE).
#' @param min.dotSize Numeric, minimum size for the dots (default: 2).
#' @param max.dotSize Numeric, maximum size for the dots (default: 10).
#' @param pos_color Character, color for positive NES dots (default: "orange").
#' @param neg_color Character, color for negative NES dots (default: "skyblue").
#'
#' @return A ggplot2 object representing the faceted GSEA dotplot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter arrange desc mutate group_by slice_max ungroup bind_rows case_when
#' @importFrom stringr str_replace_all str_to_sentence str_to_title
#' @importFrom methods is slot
#' @importFrom stats reorder
#' @importFrom rlang .data
#' @importFrom utils head
#'
#' @examples
#' # Assuming gsea_res is a valid gseaResult object
#' # Basic usage:
#' p <- gsea_dotplot_facet(gsea_res)
#' # print(p)
#'
#' # Show top 5, different colors:
#' p2 <- gsea_dotplot_facet(gsea_res, showCategory = 5, pos_color = "firebrick", neg_color = "steelblue")
#' # print(p2)
#' # To save: ggsave("my_gsea_facet_dotplot.png", p2)

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


gsea_dotplot_facet <- function(gsea_obj,
                               showCategory = 10,
                               base_font_size = 10, # Renamed
                               title = "GSEA Faceted Dotplot",
                               padj_cutoff = 0.05, # Renamed
                               replace_ = TRUE,
                               capitalize_1 = FALSE,
                               capitalize_all = FALSE,
                               min.dotSize = 2, # Added
                               max.dotSize = 10, # Added
                               pos_color = "orange", # Added
                               neg_color = "skyblue") { # Added

    # --- Input Validation ---
    if (!methods::is(gsea_obj, "gseaResult")) {
         stop("Input `gsea_obj` must be a gseaResult object from clusterProfiler.")
    }
    if (!methods::.hasSlot(gsea_obj, "result") || !is.data.frame(gsea_obj@result) || nrow(gsea_obj@result) == 0) {
        stop("Input `gsea_obj` has an invalid or empty result slot.")
    }
    required_cols <- c("Description", "NES", "p.adjust")
    if (!all(required_cols %in% colnames(gsea_obj@result))) {
        stop(sprintf("Input `gsea_obj@result` is missing required columns: %s",
                     paste(setdiff(required_cols, colnames(gsea_obj@result)), collapse=", ")))
    }
    # ------------------------

    # Extract results data and add direction
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        dplyr::filter(.data$p.adjust < padj_cutoff) %>% # Use p.adjust
        dplyr::mutate(Direction = ifelse(.data$NES > 0, "Upregulated", "Downregulated"))

    # Check if any pathways remain after filtering
    if (nrow(gsea_data) == 0) {
        warning(sprintf("No pathways found with p.adjust < %f. Returning empty plot.", padj_cutoff))
        return(ggplot() + labs(title = paste(title, "(No significant pathways)")))
    }

    # Process pathway descriptions
    if (replace_) {
        gsea_data$Description <- stringr::str_replace_all(gsea_data$Description, "_", " ")
    }
    if (capitalize_1) {
        gsea_data$Description <- stringr::str_to_sentence(gsea_data$Description)
    }
    if (capitalize_all) {
        gsea_data$Description <- stringr::str_to_title(gsea_data$Description)
    }

    # Get top N for each direction based on absolute NES
    plot_data <- gsea_data %>%
        dplyr::group_by(.data$Direction) %>%
        dplyr::slice_max(order_by = abs(.data$NES), n = showCategory) %>%
        dplyr::ungroup() %>%
        # Order pathways globally by NES for consistent y-axis
        dplyr::arrange(.data$NES) %>%
        dplyr::mutate(Description = factor(.data$Description, levels = unique(.data$Description)))

     # Handle cases where p.adjust is 0 or very small for scaling
    plot_data <- plot_data %>%
        dplyr::mutate(negLog10pAdj = -log10(.data$p.adjust),
                      negLog10pAdj_capped = pmax(.data$negLog10pAdj, 0)) # Ensure non-negative size

    # Determine size limits, handling potential Inf values
    size_values <- plot_data$negLog10pAdj_capped[is.finite(plot_data$negLog10pAdj_capped)]
    size_limits <- if(length(size_values) > 0) range(size_values, na.rm = TRUE) else c(0, 1)


    # Create faceted plot
    p <- ggplot(plot_data, aes(x = .data$NES, y = .data$Description)) + # Plot NES directly
        geom_point(aes(size = .data$negLog10pAdj_capped, color = .data$Direction)) +
        scale_color_manual(values = c("Upregulated" = pos_color, "Downregulated" = neg_color)) +
        scale_size_continuous(name = bquote(-log[10](p.adjust)),
                              range = c(min.dotSize, max.dotSize),
                              limits = size_limits) +
        facet_grid(Direction ~ ., scales = "free_y", space = "free_y") + # Allow y-axes to be independent
        labs(
            title = title,
            x = "Normalized Enrichment Score (NES)", # Updated x-axis label
            y = NULL,
            color = "Direction" # Legend title for color
        ) +
        custom_minimal_theme_with_grid() +
        theme(
            axis.text.y = element_text(size = rel(0.9) * base_font_size),
            axis.title.x = element_text(size = rel(1) * base_font_size),
            axis.text.x = element_text(size = rel(0.9) * base_font_size),
            plot.title = element_text(hjust = 0.5, size = rel(1.1) * base_font_size),
            legend.title = element_text(size = rel(0.9) * base_font_size),
            legend.text = element_text(size = rel(0.8) * base_font_size),
            strip.text = element_text(size = rel(1) * base_font_size, face = "bold"),
            strip.background = element_rect(fill = "grey90", color = "grey70"), # Adjusted strip background
            panel.spacing = unit(0.5, "lines") # Reduced panel spacing
        )

    # Return the plot object
    return(p)
}
