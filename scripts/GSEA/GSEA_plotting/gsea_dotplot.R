#' Create a Customizable GSEA Dotplot
#'
#' This function generates a visually appealing and informative dotplot from GSEA results,
#' with flexible filtering, sorting, and text processing options.
#'
#' @param gsea_obj A GSEA result object from clusterProfiler.
#' @param showCategory Integer, number of categories to display (default: 10).
#' @param font.size Numeric, base font size for the plot (default: 10).
#' @param title Character, plot title (default: "GSEA Dotplot").
#' @param replace_ Logical, whether to replace underscores with spaces in descriptions (default: TRUE).
#' @param capitalize_1 Logical, whether to capitalize the first word in descriptions (default: TRUE).
#' @param capitalize_all Logical, whether to capitalize all words in descriptions (default: FALSE).
#' @param filterBy Character, method to filter results: "qvalue", "NES", "NES_positive", or "NES_negative" (default: "qvalue").
#' @param sortBy Character, method to sort results: "GeneRatio" or "qvalue" (default: "GeneRatio").
#' @param q_cut Numeric, q-value cutoff for significance (default: 0.05).
#' @param min.dotSize Numeric, minimum dot size in the plot (default: 2).
#' @param wrap_text Logical, whether to wrap long pathway descriptions (default: TRUE).
#' @param wrap_width Integer, width for text wrapping (default: 40).
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 10).
#' @param height Numeric, height of the saved plot in inches (default: 7).
#' @param dpi Numeric, resolution of the saved plot (default: 300).
#'
#' @return A ggplot2 object representing the GSEA dotplot.
#' @export
#'
#' @examples
#' # Basic usage
#' gsea_dotplot(gsea_obj)
#'
#' # Customized filtering and sorting
#' gsea_dotplot(gsea_obj, showCategory = 15, filterBy = "NES_positive", sortBy = "GeneRatio")
#'
#' # Save the plot
#' gsea_dotplot(gsea_obj, title = "My GSEA Results", save_plot = TRUE, output_dir = "results/plots/")
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)

#' Smart text wrapping function for pathway descriptions
#'
#' @param text Character string to wrap
#' @param width Integer, maximum width before wrapping (default: 40)
#' @return Character string with newlines inserted for wrapping
#' @keywords internal
smart_wrap <- function(text, width = 40) {
    words <- unlist(strsplit(text, " "))
    total_chars <- nchar(text)
    
    if (total_chars > width) {
        mid_point <- length(words) %/% 2
        first_half <- paste(words[1:mid_point], collapse = " ")
        second_half <- paste(words[(mid_point + 1):length(words)], collapse = " ")
        return(paste(first_half, second_half, sep = "\n"))
    }
    return(text)
}

gsea_dotplot <- function(gsea_obj, 
                         showCategory = 10, 
                         font.size = 10, 
                         title = "GSEA Dotplot",
                         replace_ = TRUE, 
                         capitalize_1 = TRUE, 
                         capitalize_all = FALSE,
                         filterBy = "qvalue",
                         sortBy = "GeneRatio",
                         q_cut = 0.05,
                         min.dotSize = 2,
                         wrap_text = TRUE,
                         wrap_width = 40,
                         save_plot = FALSE,
                         output_dir = "plots/",
                         width = 10,
                         height = 7,
                         dpi = 300) {
  
  # Extract the result data frame from the GSEA object
  gsea_data <- as.data.frame(gsea_obj@result)
  
  # Calculate the gene count from 'core_enrichment' by counting '/' and adding 1 (number of genes)
  gene_count <- gsea_data %>%
    group_by(ID) %>%
    summarise(count = sum(str_count(core_enrichment, "/")) + 1)
  
  # Merge gene counts with the original GSEA result and calculate GeneRatio
  gsea_data <- left_join(gsea_data, gene_count, by = "ID") %>%
    mutate(GeneRatio = count / setSize)
  
  # Modify Description field based on function arguments
  if (replace_) {
    gsea_data$Description <- gsea_data$Description %>% 
      str_replace_all("_", " ")   # Replace "_" with " " if 'replace' is TRUE
  }
  
  if (capitalize_1) {
    gsea_data$Description <- gsea_data$Description %>%
      str_to_sentence()           # Capitalize only the first word if 'capitalize_1' is TRUE
  }
  
  if (capitalize_all) {
    gsea_data$Description <- gsea_data$Description %>%
      str_to_title()              # Capitalize all words if 'capitalize_all' is TRUE
  }
  
  # Apply text wrapping if requested
  if (wrap_text) {
    gsea_data$Description <- sapply(gsea_data$Description, smart_wrap, wrap_width)
  }
  
  # Filter for significant pathways
  gsea_data_filtered <- gsea_data %>%
    filter(qvalue < q_cut) %>%
    mutate(NES_sign = ifelse(NES > 0, "Positive NES", "Negative NES"))
  
  # Filter logic based on 'filterBy' argument
  if (filterBy == "qvalue") {
    # Sort by qvalue (default behavior)
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(qvalue) %>%           # Sort by qvalue (ascending)
      head(showCategory)
  } else if (filterBy == "NES") {
    # Sort by absolute NES value (strongest enrichment)
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(abs(NES))) %>%    # Sort by absolute NES (descending)
      head(showCategory)
  } else if (filterBy == "NES_positive") {
    # Filter and sort by positive NES values only
    gsea_data_filtered <- gsea_data_filtered %>%
      filter(NES > 0) %>%           # Filter positive NES
      arrange(desc(NES)) %>%        # Sort by NES (descending)
      head(showCategory)
  } else if (filterBy == "NES_negative") {
    # Filter and sort by negative NES values only
    gsea_data_filtered <- gsea_data_filtered %>%
      filter(NES < 0) %>%           # Filter negative NES
      arrange(NES) %>%              # Sort by NES (ascending, more negative first)
      head(showCategory)
  }

  # Sort the filtered data based on the sortBy parameter
  if (sortBy == "GeneRatio") {
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(GeneRatio))
  } else if (sortBy == "qvalue") {
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(qvalue)
  } else {
    warning("Invalid sortBy parameter. Defaulting to GeneRatio.")
    gsea_data_filtered <- gsea_data_filtered %>%
      arrange(desc(GeneRatio))
  }

  # Take the top showCategory entries after sorting
  gsea_data_filtered <- gsea_data_filtered %>%
    head(showCategory)

  # Create custom dotplot using ggplot2
  p <- ggplot(gsea_data_filtered, aes(x = GeneRatio, y = reorder(Description, !!sym(sortBy)))) +
    geom_point(aes(size = -log10(qvalue), color = NES_sign)) +
    scale_color_manual(values = c("Positive NES" = "orange", "Negative NES" = "skyblue")) +
    
    # Use scale_size_continuous to set visual size limits for the dots
    scale_size_continuous(range = c(min.dotSize, 10),
                          limits = c(min(-log10(gsea_data_filtered$qvalue)), 
                                     max(-log10(gsea_data_filtered$qvalue))),
                          name = "-log10(qvalue)") +
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      color = "NES",
      size = "-log10(qvalue)"
    ) +
    
    custom_minimal_theme_with_grid() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),  # White background with black border
      plot.background = element_rect(fill = "white", color = NA),  # White plot background with no border
      axis.text.y = element_text(size = font.size, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = font.size + 2),
      axis.text.x = element_text(size = font.size),
      legend.position = "right"
    )

  # Save the plot if requested
  if (save_plot) {
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create filename from title
    filename <- file.path(output_dir, paste0(gsub(" ", "_", title), ".pdf"))
    
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
