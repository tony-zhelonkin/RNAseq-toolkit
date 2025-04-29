#' Create a Faceted GSEA Dotplot with Separate Up/Down Regulation
#'
#' This function creates a faceted dotplot visualization of GSEA results,
#' separating upregulated and downregulated pathways into different facets.
#'
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param showCategory Integer, number of categories to display per facet (default: 10).
#' @param font.size Numeric, base font size for the plot (default: 10).
#' @param title Character, plot title (default: "GSEA Faceted Dotplot").
#' @param q_cut Numeric, q-value cutoff for significance (default: 0.05).
#' @param replace_ Logical, whether to replace underscores with spaces in descriptions (default: TRUE).
#' @param capitalize_1 Logical, whether to capitalize the first word in descriptions (default: FALSE).
#' @param capitalize_all Logical, whether to capitalize all words in descriptions (default: FALSE).
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 12).
#' @param height Numeric, height of the saved plot in inches (default: 8).
#' @param dpi Numeric, resolution of the saved plot (default: 300).
#'
#' @return A ggplot2 object representing the faceted GSEA dotplot.
#' @export
#'
#' @examples
#' # Basic usage
#' gsea_dotplot_facet(gsea_obj)
#'
#' # Display more pathways per facet
#' gsea_dotplot_facet(gsea_obj, showCategory = 15)
#'
#' # Save the plot
#' gsea_dotplot_facet(gsea_obj, title = "My GSEA Results", save_plot = TRUE)
library(ggplot2)
library(dplyr)
library(stringr)

gsea_dotplot_facet <- function(gsea_obj, 
                               showCategory = 10, 
                               font.size = 10, 
                               title = "GSEA Faceted Dotplot",
                               q_cut = 0.05,
                               replace_ = TRUE,
                               capitalize_1 = FALSE,
                               capitalize_all = FALSE,
                               save_plot = FALSE,
                               output_dir = "plots/",
                               width = 12,
                               height = 8,
                               dpi = 300) {
    
    # Extract results data
    gsea_data <- as.data.frame(gsea_obj@result)
    
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
    
    # Get positive and negative results
    pos_data <- gsea_data %>%
        filter(NES > 0, qvalue < q_cut) %>%
        arrange(desc(NES)) %>%
        head(showCategory)
    
    neg_data <- gsea_data %>%
        filter(NES < 0, qvalue < q_cut) %>%
        arrange(NES) %>%
        head(showCategory)
    
    # Combine data
    plot_data <- rbind(pos_data, neg_data) %>%
        mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated"))
    
    # Create faceted plot
    p <- ggplot(plot_data, aes(x = abs(NES), y = reorder(Description, NES))) +
        geom_point(aes(size = -log10(qvalue), color = Direction)) +
        scale_color_manual(values = c("Upregulated" = "orange", "Downregulated" = "skyblue")) +
        scale_size_continuous(name = "-log10(qvalue)") +
        facet_grid(Direction ~ ., scales = "free_y", space = "free") +
        labs(
            title = title,
            x = "Normalized Enrichment Score (absolute value)",
            y = NULL
        ) +
        custom_minimal_theme_with_grid() +
        theme(
            text = element_text(size = font.size),
            strip.text = element_text(size = font.size + 1, face = "bold"),
            strip.background = element_rect(fill = "white", color = "black"),
            panel.spacing = unit(1, "lines")
        )
    
    # Save the plot if requested
    if (save_plot) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Create filename from title
        filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_facet.pdf"))
        
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
