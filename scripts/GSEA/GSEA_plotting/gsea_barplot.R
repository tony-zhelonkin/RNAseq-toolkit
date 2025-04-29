#' Create a GSEA Normalized Enrichment Score Barplot
#'
#' This function creates a barplot visualization of GSEA results, showing the
#' Normalized Enrichment Scores (NES) for top pathways.
#'
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param title Character, plot title (default: "GSEA NES Barplot").
#' @param q_cut Numeric, q-value cutoff for significance (default: 0.05).
#' @param top_n Integer, number of top pathways to display (default: 30).
#' @param replace_ Logical, whether to replace underscores with spaces in descriptions (default: TRUE).
#' @param capitalize_1 Logical, whether to capitalize the first word in descriptions (default: FALSE).
#' @param capitalize_all Logical, whether to capitalize all words in descriptions (default: FALSE).
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 10).
#' @param height Numeric, height of the saved plot in inches (default: 8).
#' @param dpi Numeric, resolution of the saved plot (default: 300).
#'
#' @return A ggplot2 object representing the GSEA barplot.
#' @export
#'
#' @examples
#' # Basic usage
#' gsea_barplot(gsea_obj)
#'
#' # Display more pathways
#' gsea_barplot(gsea_obj, top_n = 50)
#'
#' # Save the plot
#' gsea_barplot(gsea_obj, title = "My GSEA Results", save_plot = TRUE, output_dir = "results/plots/")
library(ggplot2)
library(dplyr)
library(stringr)

gsea_barplot <- function(gsea_obj, 
                         title = "GSEA NES Barplot", 
                         q_cut = 0.05, 
                         top_n = 30,
                         replace_ = TRUE,
                         capitalize_1 = FALSE,
                         capitalize_all = FALSE,
                         save_plot = FALSE,
                         output_dir = "plots/",
                         width = 10, 
                         height = 8,
                         dpi = 300) {
    # Extract significant pathways
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        filter(qvalue < q_cut) %>%
        arrange(desc(abs(NES))) %>%
        head(top_n)
    
    # Process pathway descriptions
    if (replace_) {
        gsea_data$Description <- gsea_data$Description %>% 
            str_replace_all("_", " ")
    }
    
    if (capitalize_1) {
        gsea_data$Description <- gsea_data$Description %>%
            str_to_sentence()
    }
    
    if (capitalize_all) {
        gsea_data$Description <- gsea_data$Description %>%
            str_to_title()
    }
    
    # Sort by NES for the plot
    gsea_data <- gsea_data %>% arrange(NES)
    
    # Create the plot
    p <- ggplot(gsea_data, aes(x = reorder(Description, NES), y = NES)) +
        geom_bar(stat = "identity", aes(fill = NES > 0)) +
        scale_fill_manual(values = c("skyblue", "orange")) +
        coord_flip() +
        labs(title = title,
             x = NULL,
             y = "Normalized Enrichment Score") +
        custom_minimal_theme_with_grid() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(lineheight = 0.8)
        )
    
    # Save the plot if requested
    if (save_plot) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Create filename from title
        filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_barplot.pdf"))
        
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
