#' Create GSEA Normalized Enrichment Score (NES) Barplot
#'
#' Creates a horizontal barplot visualizing the NES for the top significant pathways
#' from a `gseaResult` object. Bars are colored based on the sign of NES (positive/negative).
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#'
#' @param gsea_obj A `gseaResult` object from `clusterProfiler`.
#' @param title Character, plot title (default: "GSEA NES Barplot").
#' @param padj_cutoff Numeric, adjusted p-value cutoff for filtering significant pathways
#'        (default: 0.05). Uses the `p.adjust` column from `gsea_obj@result`.
#' @param top_n Integer, the maximum number of top pathways (ranked by absolute NES)
#'        to display after filtering by `padj_cutoff` (default: 30).
#' @param replace_ Logical, if TRUE (default), replace underscores "_" with spaces " "
#'        in pathway descriptions for display.
#' @param capitalize_1 Logical, if TRUE, capitalize the first letter of each pathway
#'        description (default: FALSE). Uses `stringr::str_to_sentence`.
#' @param capitalize_all Logical, if TRUE, capitalize the first letter of each word
#'        in pathway descriptions (default: TRUE). Uses `stringr::str_to_title`.
#' @param preserve_case Logical, if TRUE, preserve case of common abbreviations like DNA, RNA, etc. (default: TRUE).
#' @param pos_color Character, color for positive NES bars (default: "#D73027").
#' @param neg_color Character, color for negative NES bars (default: "#4575B4").
#' @param base_font_size Numeric, base font size for plot text elements (default: 8).
#' @param width Numeric, width of the plot in inches (default: 7).
#' @param height Numeric, height of the plot in inches (default: 5).
#'
#' @return A ggplot2 object representing the GSEA barplot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter arrange desc head mutate
#' @importFrom stringr str_replace_all str_to_sentence str_to_title
#' @importFrom methods is slot
#' @importFrom stats reorder
#' @importFrom utils head
#'
#' @examples
#' # Assuming gsea_res is a valid gseaResult object
#' # Basic usage:
#' p <- gsea_barplot(gsea_res)
#' # print(p)
#'
#' # Display top 10, different colors:
#' p2 <- gsea_barplot(gsea_res, top_n = 10, padj_cutoff = 0.1,
#'                    pos_color = "red", neg_color = "blue")
#' # print(p2)
#' # To save: ggsave("my_gsea_barplot.png", p2)

# Function to safely source a script if it exists
source_safe <- function(path) {
  base_dir <- "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/Kelsey_Followup"
  full_path <- file.path(base_dir, path)
  if (file.exists(full_path)) {
    source(full_path)
    return(TRUE)
  } else {
    warning("Custom theme script not found: ", full_path, ". Using default theme_minimal().")
    # Define a placeholder if theme is missing
    custom_minimal_theme_with_grid <<- function() ggplot2::theme_minimal()
    return(FALSE)
  }
}
# Source the custom theme - needs to be available when the package/script is loaded or run
source_safe("R_GSEA_visualisations/scripts/custom_minimal_theme.R")


gsea_barplot <- function(gsea_obj,
                         title = "GSEA NES Barplot",
                         padj_cutoff = 0.05,
                         top_n = 30,
                         replace_ = TRUE,
                         capitalize_1 = FALSE,
                         capitalize_all = TRUE,
                         preserve_case = TRUE,
                         pos_color = "#D73027", # Colorblind-friendly red
                         neg_color = "#4575B4", # Colorblind-friendly blue
                         base_font_size = 8,
                         width = 7,
                         height = 5) {

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

    # Extract, filter, and prepare data
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        dplyr::filter(p.adjust < padj_cutoff) # Use p.adjust

    # Check if any pathways remain after filtering
    if (nrow(gsea_data) == 0) {
        warning(sprintf("No pathways found with p.adjust < %f. Returning empty plot.", padj_cutoff))
        # Return a blank plot with title
        return(ggplot2::ggplot() + ggplot2::labs(title = paste(title, "(No significant pathways)")))
    }

    gsea_data <- gsea_data %>%
        dplyr::arrange(dplyr::desc(abs(NES))) %>%
        utils::head(top_n)

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

    # Sort by NES for plotting order
    gsea_data <- gsea_data %>% dplyr::arrange(NES)

    # Create the plot
    p <- ggplot2::ggplot(gsea_data, ggplot2::aes(x = stats::reorder(Description, NES), y = NES)) +
        ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = NES > 0)) +
        ggplot2::scale_fill_manual(values = c(`FALSE` = neg_color, `TRUE` = pos_color)) + # Use parameters
        ggplot2::coord_flip() +
        ggplot2::labs(title = title,
             x = NULL, # No x-axis label needed
             y = "NES") + # Simplified y-axis label
        custom_minimal_theme_with_grid(base_size = base_font_size) + # Apply theme with custom base size
        ggplot2::theme(
            legend.position = "none", # No legend needed for fill
            panel.grid = ggplot2::element_blank(), # Remove grid lines
            axis.text.y = ggplot2::element_text(size = base_font_size), 
            axis.title.x = ggplot2::element_text(size = base_font_size),
            axis.text.x = ggplot2::element_text(size = base_font_size),
            plot.title = ggplot2::element_text(hjust = 0.5, size = base_font_size * 1.2),
            plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
        )

    # Set plot dimensions as attributes for later use
    attr(p, "width") <- width
    attr(p, "height") <- height
    
    # Return the plot object
    return(p)
}
