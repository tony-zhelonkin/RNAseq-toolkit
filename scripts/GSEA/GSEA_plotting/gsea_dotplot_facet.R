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
#' @param base_font_size Numeric, base font size for plot text elements (default: 8).
#' @param title Character, plot title (default: "GSEA Faceted Dotplot").
#' @param padj_cutoff Numeric, adjusted p-value cutoff for filtering significant pathways
#'        (default: 0.05). Uses the `p.adjust` column.
#' @param replace_ Logical, if TRUE (default), replace underscores "_" with spaces " " in descriptions.
#' @param capitalize_1 Logical, if TRUE, capitalize the first letter of descriptions (default: FALSE).
#' @param capitalize_all Logical, if TRUE, capitalize the first letter of each word (default: TRUE).
#' @param preserve_case Logical, if TRUE, preserve case of common abbreviations like DNA, RNA, etc. (default: TRUE).
#' @param min.dotSize Numeric, minimum size for the dots (default: 1).
#' @param max.dotSize Numeric, maximum size for the dots (default: 5).
#' @param pos_color Character, color for positive NES dots (default: "#D73027").
#' @param neg_color Character, color for negative NES dots (default: "#4575B4").
#' @param wrap_text Logical, whether to wrap long pathway descriptions (default: TRUE).
#' @param wrap_width Integer, maximum width for text wrapping (default: 50).
#' @param width Numeric, width of the plot in inches (default: 7).
#' @param height Numeric, height of the plot in inches (default: 7).
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
  base_dir <- "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/Kelsey_Followup"
  full_path <- file.path(base_dir, path)
  if (file.exists(full_path)) {
    source(full_path)
    return(TRUE)
  } else {
    warning("Custom theme script not found: ", full_path, ". Using default theme_minimal().")
    custom_minimal_theme_with_grid <<- function() ggplot2::theme_minimal() # Define placeholder
    return(FALSE)
  }
}
# Source the custom theme
source_safe("R_GSEA_visualisations/scripts/custom_minimal_theme.R")

#' Improved text wrapping function (internal)
smart_wrap <- function(text, width = 50) {
    words <- unlist(strsplit(text, " "))
    if (length(words) == 0) return(text) # Handle empty strings
    total_chars <- nchar(text) - (length(words) - 1) # Approx chars without spaces

    # For very long pathway names, use more aggressive wrapping
    if (total_chars > width * 1.5) {
        # Try to split into three parts for very long names
        char_per_line <- total_chars / 3
        lines <- character(0)
        current_line <- ""
        current_chars <- 0
        
        for (word in words) {
            word_chars <- nchar(word)
            if (current_chars + word_chars > char_per_line && current_chars > 0) {
                lines <- c(lines, current_line)
                current_line <- word
                current_chars <- word_chars
            } else {
                if (current_chars > 0) current_line <- paste(current_line, word)
                else current_line <- word
                current_chars <- current_chars + word_chars
            }
        }
        
        if (current_chars > 0) lines <- c(lines, current_line)
        return(paste(lines, collapse = "\n"))
    } else if (total_chars > width) {
        # For moderately long names, split into two parts
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


gsea_dotplot_facet <- function(gsea_obj,
                               showCategory = 10,
                               base_font_size = 8,
                               title = "GSEA Faceted Dotplot",
                               padj_cutoff = 0.05,
                               replace_ = TRUE,
                               capitalize_1 = FALSE,
                               capitalize_all = TRUE,
                               preserve_case = TRUE,
                               min.dotSize = 1,
                               max.dotSize = 5,
                               pos_color = "#D73027", # Colorblind-friendly red
                               neg_color = "#4575B4", # Colorblind-friendly blue
                               wrap_text = TRUE,
                               wrap_width = 50,
                               width = 7,
                               height = 7) {

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
        dplyr::mutate(Direction = ifelse(.data$NES > 0, "Up", "Down"))

    # Check if any pathways remain after filtering
    if (nrow(gsea_data) == 0) {
        warning(sprintf("No pathways found with p.adjust < %f. Returning empty plot.", padj_cutoff))
        return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No significant pathways)")))
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
    
    # Preserve case for common abbreviations
    if (preserve_case) {
        common_abbr <- c("DNA", "RNA", "ATP", "ADP", "GTP", "GDP", "NADH", "NADPH", 
                         "FADH", "mRNA", "tRNA", "rRNA", "miRNA", "siRNA", "lncRNA",
                         "NF-kB", "TNF", "IL", "IFN", "TGF", "EGF", "VEGF", "IGF",
                         "MHC", "TCR", "BCR", "CD", "NK", "DC", "Th", "Treg")
        
        for (abbr in common_abbr) {
            # Case-insensitive replacement to ensure we catch all variations
            pattern <- paste0("\\b", tolower(abbr), "\\b")
            replacement <- abbr
            gsea_data$Description <- stringr::str_replace_all(
                gsea_data$Description, 
                stringr::regex(pattern, ignore_case = TRUE), 
                replacement
            )
        }
    }
    
    # Apply text wrapping if enabled
    if (wrap_text) {
        # Use Vectorize for efficiency if many descriptions
        wrap_fun <- Vectorize(function(txt) smart_wrap(txt, width = wrap_width))
        gsea_data$Description <- wrap_fun(gsea_data$Description)
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


    # Create faceted plot - use GeneRatio instead of NES for x-axis
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = .data$Description)) +
        # Simplified point geometry without complex aesthetics
        ggplot2::geom_point(
            ggplot2::aes(
                size = .data$negLog10pAdj_capped, 
                color = .data$Direction
            )
        ) +
        ggplot2::scale_color_manual(values = c("Up" = pos_color, "Down" = neg_color)) +
        ggplot2::scale_size_continuous(name = bquote(-log[10](q-value)),
                              range = c(min.dotSize, max.dotSize),
                              limits = size_limits) +
        ggplot2::facet_grid(Direction ~ ., scales = "free_y", space = "free_y") + # Allow y-axes to be independent
        ggplot2::labs(
            title = title,
            x = "Gene Ratio", # Changed x-axis label to match other dotplots
            y = NULL,
            color = "Direction" # Legend title for color
        ) +
        custom_minimal_theme_with_grid(base_size = base_font_size) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(), # Remove grid lines
            # Improved text spacing and margins
            axis.text.y = ggplot2::element_text(
                size = ggplot2::rel(0.9) * base_font_size,
                hjust = 1,
                margin = ggplot2::margin(r = 10) # More space between y-axis text and plot
            ),
            axis.title.x = ggplot2::element_text(
                size = ggplot2::rel(1) * base_font_size,
                margin = ggplot2::margin(t = 10) # More space above x-axis title
            ),
            axis.text.x = ggplot2::element_text(
                size = ggplot2::rel(0.9) * base_font_size,
                margin = ggplot2::margin(t = 5) # More space above x-axis text
            ),
            plot.title = ggplot2::element_text(
                hjust = 0.5, 
                size = ggplot2::rel(1.1) * base_font_size,
                margin = ggplot2::margin(b = 10) # More space below title
            ),
            legend.title = ggplot2::element_text(size = ggplot2::rel(0.9) * base_font_size),
            legend.text = ggplot2::element_text(size = ggplot2::rel(0.8) * base_font_size),
            legend.margin = ggplot2::margin(l = 10), # More space to the left of legend
            strip.text = ggplot2::element_text(
                size = ggplot2::rel(1) * base_font_size, 
                face = "plain",  # Regular font instead of bold
                margin = ggplot2::margin(b = 5) # More space below strip text
            ),
            strip.background = ggplot2::element_rect(fill = "grey90", color = "grey70"),
            panel.spacing = ggplot2::unit(1, "lines"), # Increased panel spacing
            # Increased overall plot margins
            plot.margin = ggplot2::margin(15, 15, 15, 15) # Larger margins all around
        )

    # Set plot dimensions as attributes for later use
    attr(p, "width") <- width
    attr(p, "height") <- height
    
    # Return the plot object
    return(p)
}
