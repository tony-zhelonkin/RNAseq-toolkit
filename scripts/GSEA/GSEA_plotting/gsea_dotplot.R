#' Create Customizable GSEA Dotplot
#'
#' Generates a dotplot from GSEA results, showing GeneRatio vs. Pathway Description.
#' Dot size represents significance (`-log10(p.adjust)`), and color indicates NES sign.
#' Includes options for filtering, sorting, text formatting, and wrapping.
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param showCategory Integer, number of categories to display after filtering (default: 10).
#' @param base_font_size Numeric, base font size for plot text elements (default: 10).
#' @param title Character, plot title (default: "GSEA Dotplot").
#' @param replace_ Logical, if TRUE (default), replace underscores "_" with spaces " " in descriptions.
#' @param capitalize_1 Logical, if TRUE, capitalize the first letter of descriptions (default: TRUE).
#' @param capitalize_all Logical, if TRUE, capitalize the first letter of each word (default: FALSE).
#' @param filterBy Character, method to rank pathways before selecting top `showCategory`:
#'        "p.adjust" (ranks by p.adjust, ascending), "NES" (ranks by absolute NES, descending),
#'        "NES_positive" (filters for NES > 0, ranks by NES descending),
#'        "NES_negative" (filters for NES < 0, ranks by NES ascending). Default: "p.adjust".
#' @param sortBy Character, method to sort the final `showCategory` pathways for display:
#'        "GeneRatio" (descending) or "p.adjust" (ascending). Default: "GeneRatio".
#' @param padj_cutoff Numeric, adjusted p-value cutoff for initially filtering pathways
#'        (default: 0.05). Uses the `p.adjust` column.
#' @param min.dotSize Numeric, minimum size for the dots (default: 2).
#' @param max.dotSize Numeric, maximum size for the dots (default: 10).
#' @param wrap_text Logical, whether to wrap long pathway descriptions (default: TRUE).
#' @param wrap_width Integer, maximum width for text wrapping (default: 40).
#' @param pos_color Character, color for positive NES dots (default: "orange").
#' @param neg_color Character, color for negative NES dots (default: "skyblue").
#'
#' @return A ggplot2 object representing the GSEA dotplot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2 scales
#' @importFrom dplyr %>% filter select mutate arrange desc left_join group_by summarise case_when sym
#' @importFrom stringr str_count str_replace_all str_to_sentence str_to_title
#' @importFrom methods is slot
#' @importFrom stats reorder setNames
#' @importFrom rlang .data := !!
#' @importFrom utils head
#'
#' @examples
#' # Assuming gsea_res is a valid gseaResult object
#' # Basic usage (filter by p.adjust, show top 10, sort by GeneRatio):
#' p <- gsea_dotplot(gsea_res)
#' # print(p)
#'
#' # Filter by positive NES, show top 15, sort by p.adjust:
#' p2 <- gsea_dotplot(gsea_res, showCategory = 15, filterBy = "NES_positive", sortBy = "p.adjust")
#' # print(p2)
#' # To save: ggsave("my_gsea_dotplot.png", p2)

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


#' Smart text wrapping function (internal)
smart_wrap <- function(text, width = 40) {
    words <- unlist(strsplit(text, " "))
    if (length(words) == 0) return(text) # Handle empty strings
    total_chars <- nchar(text) - (length(words) - 1) # Approx chars without spaces

    if (total_chars > width) {
        # Try to find a split point near the middle
        char_count <- 0
        split_point <- 0
        for (i in 1:length(words)) {
            char_count <- char_count + nchar(words[i])
            if (char_count >= total_chars / 2 && i < length(words)) {
                split_point <- i
                break
            }
        }
        if (split_point == 0) split_point <- max(1, length(words) %/% 2) # Fallback split

        first_half <- paste(words[1:split_point], collapse = " ")
        second_half <- paste(words[(split_point + 1):length(words)], collapse = " ")
        return(paste(first_half, second_half, sep = "\n"))
    }
    return(text)
}


gsea_dotplot <- function(gsea_obj,
                         showCategory = 10,
                         base_font_size = 10, # Renamed
                         title = "GSEA Dotplot",
                         replace_ = TRUE,
                         capitalize_1 = TRUE,
                         capitalize_all = FALSE,
                         filterBy = "p.adjust", # Changed default, renamed qvalue
                         sortBy = "GeneRatio",
                         padj_cutoff = 0.05, # Renamed
                         min.dotSize = 2,
                         max.dotSize = 10, # Added
                         wrap_text = TRUE,
                         wrap_width = 40,
                         pos_color = "orange", # Added
                         neg_color = "skyblue") { # Added

    # --- Input Validation ---
    if (!methods::is(gsea_obj, "gseaResult")) {
         stop("Input `gsea_obj` must be a gseaResult object from clusterProfiler.")
    }
    if (!methods::.hasSlot(gsea_obj, "result") || !is.data.frame(gsea_obj@result) || nrow(gsea_obj@result) == 0) {
        stop("Input `gsea_obj` has an invalid or empty result slot.")
    }
    required_cols <- c("ID", "Description", "NES", "p.adjust", "core_enrichment", "setSize")
    if (!all(required_cols %in% colnames(gsea_obj@result))) {
        stop(sprintf("Input `gsea_obj@result` is missing required columns: %s",
                     paste(setdiff(required_cols, colnames(gsea_obj@result)), collapse=", ")))
    }
    if (!filterBy %in% c("p.adjust", "NES", "NES_positive", "NES_negative")) {
        warning("Invalid `filterBy` argument. Using 'p.adjust'.")
        filterBy <- "p.adjust"
    }
     if (!sortBy %in% c("GeneRatio", "p.adjust")) {
        warning("Invalid `sortBy` argument. Using 'GeneRatio'.")
        sortBy <- "GeneRatio"
    }
    # ------------------------

    # Extract the result data frame
    gsea_data <- as.data.frame(gsea_obj@result)

    # Calculate GeneRatio
    gsea_data <- gsea_data %>%
        dplyr::mutate(
            count = stringr::str_count(.data$core_enrichment, "/") + ifelse(nchar(.data$core_enrichment) > 0, 1, 0),
            GeneRatio = .data$count / .data$setSize,
            negLog10pAdj = -log10(.data$p.adjust) # Calculate for size mapping
        )

    # Modify Description field
    if (replace_) {
        gsea_data$Description <- stringr::str_replace_all(gsea_data$Description, "_", " ")
    }
    if (capitalize_1) {
        gsea_data$Description <- stringr::str_to_sentence(gsea_data$Description)
    }
    if (capitalize_all) {
        gsea_data$Description <- stringr::str_to_title(gsea_data$Description)
    }
    if (wrap_text) {
        # Use Vectorize for efficiency if many descriptions
        wrap_fun <- Vectorize(function(txt) smart_wrap(txt, width = wrap_width))
        gsea_data$Description <- wrap_fun(gsea_data$Description)
    }

    # --- Filtering and Ranking Logic ---
    # 1. Filter by padj_cutoff
    gsea_data_filtered <- gsea_data %>%
        dplyr::filter(.data$p.adjust < padj_cutoff)

    # Check if any pathways remain
    if (nrow(gsea_data_filtered) == 0) {
        warning(sprintf("No pathways found with p.adjust < %f. Returning empty plot.", padj_cutoff))
        return(ggplot() + labs(title = paste(title, "(No significant pathways)")))
    }

    # 2. Apply NES sign filter if needed
    if (filterBy == "NES_positive") {
        gsea_data_filtered <- gsea_data_filtered %>% dplyr::filter(.data$NES > 0)
    } else if (filterBy == "NES_negative") {
        gsea_data_filtered <- gsea_data_filtered %>% dplyr::filter(.data$NES < 0)
    }

    # Check again if pathways remain after NES filter
    if (nrow(gsea_data_filtered) == 0) {
        warning(sprintf("No pathways found matching filter criteria (padj < %f, filterBy='%s'). Returning empty plot.", padj_cutoff, filterBy))
        return(ggplot() + labs(title = paste(title, "(No matching pathways)")))
    }

    # 3. Arrange by the primary filter metric to select top categories
    if (filterBy == "p.adjust") {
        gsea_data_filtered <- gsea_data_filtered %>% dplyr::arrange(.data$p.adjust)
    } else if (filterBy == "NES" || filterBy == "NES_positive") {
        gsea_data_filtered <- gsea_data_filtered %>% dplyr::arrange(dplyr::desc(abs(.data$NES))) # Rank by abs(NES) for NES filter too
    } else if (filterBy == "NES_negative") {
        gsea_data_filtered <- gsea_data_filtered %>% dplyr::arrange(.data$NES) # Rank by NES ascending for negative
    }

    # 4. Take top 'showCategory'
    gsea_data_filtered <- gsea_data_filtered %>% utils::head(showCategory)

    # 5. Sort the final subset by the 'sortBy' metric for plotting order
    if (sortBy == "GeneRatio") {
        # Order factor levels for ggplot based on GeneRatio descending
        plot_data <- gsea_data_filtered %>%
            dplyr::mutate(Description = factor(.data$Description, levels = rev(unique(.data$Description[order(.data$GeneRatio, decreasing = FALSE)]))))
    } else { # sortBy == "p.adjust"
         # Order factor levels for ggplot based on p.adjust ascending
         plot_data <- gsea_data_filtered %>%
            dplyr::mutate(Description = factor(.data$Description, levels = rev(unique(.data$Description[order(.data$p.adjust, decreasing = TRUE)]))))
    }
    # --- End Filtering/Ranking ---


    # Handle cases where p.adjust is 0 or very small for scaling
    plot_data <- plot_data %>%
        dplyr::mutate(negLog10pAdj_capped = pmax(.data$negLog10pAdj, 0)) # Ensure non-negative size

    # Determine size limits, handling potential Inf values
    size_values <- plot_data$negLog10pAdj_capped[is.finite(plot_data$negLog10pAdj_capped)]
    size_limits <- if(length(size_values) > 0) range(size_values, na.rm = TRUE) else c(0, 1)


    # Create custom dotplot using ggplot2
    p <- ggplot(plot_data, aes(x = .data$GeneRatio, y = .data$Description)) +
        geom_point(aes(size = .data$negLog10pAdj_capped, color = .data$NES > 0)) + # Color by NES sign directly
        scale_color_manual(name = "NES Sign", values = c(`FALSE` = neg_color, `TRUE` = pos_color)) +
        scale_size_continuous(name = bquote(-log[10](p.adjust)),
                              range = c(min.dotSize, max.dotSize),
                              limits = size_limits) +
        labs(
            title = title,
            x = "Gene Ratio",
            y = NULL # No y-axis label
        ) +
        custom_minimal_theme_with_grid() +
        theme(
            axis.text.y = element_text(size = rel(0.9) * base_font_size, hjust = 1),
            plot.title = element_text(hjust = 0.5, size = rel(1.1) * base_font_size),
            axis.title.x = element_text(size = rel(1) * base_font_size),
            axis.text.x = element_text(size = rel(0.9) * base_font_size),
            legend.title = element_text(size = rel(0.9) * base_font_size),
            legend.text = element_text(size = rel(0.8) * base_font_size),
            legend.position = "right"
        )

    # Return the plot object
    return(p)
}
