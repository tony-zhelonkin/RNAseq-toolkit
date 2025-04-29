#' Create a PCA Plot from a DGEList Object
#'
#' This function performs Principal Component Analysis (PCA) on normalized expression data
#' from a DGEList object and creates a customizable PCA plot.
#'
#' @param DGE_object A DGEList object containing normalized expression data and sample metadata.
#' @param title Character, plot title (default: "PCA Plot").
#' @param limit_buffer Numeric, buffer factor for axis limits (default: 1).
#' @param x_max Numeric, maximum value for x-axis (default: max(abs(pca_data$PC1))).
#' @param y_max Numeric, maximum value for y-axis (default: max(abs(pca_data$PC2))).
#'
#' @return A ggplot2 object representing the PCA plot.
#' @export
#'
#' @examples
#' # Assuming DGE_object is a DGEList object with normalized counts and sample metadata
#' create_pca_plot(DGE_object, title = "PCA of RNA-seq Samples")
source("scripts/custom_minimal_theme.R")

create_pca_plot <- function(
  DGE_object, title = "PCA Plot",
  limit_buffer = 1,
  x_max = max(abs(pca_data$PC1)),
  y_max =max(abs(pca_data$PC2))) {
  # Calculate logCPM
  logCPM <- cpm(DGE_object, log = TRUE, prior.count = 1)

  # Calculate PCA
  pca <- prcomp(t(logCPM))

  # Prepare PCA data frame
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    condition = DGE_object$samples$group,
    organ = DGE_object$samples$organ,   # Extract organ type
    sample_id = DGE_object$samples$id,  # Extract sample ID
    cell_genotype = as.character(DGE_object$samples$group)  # Main group label
  )

  # Define custom shapes: Lung = Inverse Triangle (24), Spleen = Square (15)
  shape_mapping <- c("Lung" = 6, "Spleen" = 0)

  # Calculate variance percentages
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # Calculate axis limits with buffer
  limit_buffer <- limit_buffer

  # Create and return the plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = organ)) +
    geom_point(size = 5) +
    geom_text(aes(label = cell_genotype), vjust = -1.5, hjust = 1, size = 4, fontface = "bold", show.legend = FALSE) +  # Main group text
    #geom_text(aes(label = rownames(DGE_object$samples), color = condition), vjust = 1.5, hjust = 1, size = 3, fontface = "italic", color = "black", show.legend = FALSE) + # Sample ID below the dot
    labs(color = "Condition", shape = "Organ") +
    scale_shape_manual(values = shape_mapping) +  # Apply shape mapping
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                             "#0072B2", "#D55E00", "#CC79A7", "#000000")) +
    xlab(paste0("PC1: ", round(percentVar[1], 1), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar[2], 1), "% variance")) +
    ggtitle(title) +
    custom_minimal_theme_with_grid() +
    xlim(-x_max * limit_buffer, x_max * limit_buffer) +
    ylim(-y_max * limit_buffer, y_max * limit_buffer) +
    theme(
      legend.position = "right",  # Move legend outside to avoid overlaying
      legend.box = "vertical",    
      legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.box.margin = margin(2, 2, 2, 2),
      plot.title = element_text(hjust = 0.2, vjust = 0.2, margin = margin(b = 10)),
      plot.margin = margin(20, 20, 20, 20)
    )

  ggsave(
    filename = paste0("3_Results/imgs/PCA/", gsub(" ", "_", title), ".pdf"),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 300
  )

  return(pca_plot)
}
