#' Create a 3D PCA Plot from a DGEList Object
#'
#' This function performs Principal Component Analysis (PCA) on normalized expression data
#' from a DGEList object and creates an interactive 3D PCA plot using plotly.
#'
#' @param DGE_object A DGEList object containing normalized expression data and sample metadata.
#' @param title Character, plot title (default: "3D PCA Plot").
#'
#' @return A plotly object representing the 3D PCA plot.
#' @export
#'
#' @examples
#' # Assuming DGE_object is a DGEList object with normalized counts and sample metadata
#' create_3d_pca_plot(DGE_object, title = "3D PCA of RNA-seq Samples")
library(edgeR)
library(ggplot2)
library(plotly)

create_3d_pca_plot <- function(DGE_object, title = "3D PCA Plot") {
  # Calculate logCPM
  logCPM <- cpm(DGE_object, log = TRUE, prior.count = 1)

  # Perform PCA
  pca <- prcomp(t(logCPM))

  # Create PCA dataframe with first three principal components
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    condition = DGE_object$samples$group,
    organ = DGE_object$samples$organ,   
    sample_id = rownames(DGE_object$samples),  
    cell_genotype = as.character(DGE_object$samples$group)  
  )

  # Compute variance percentages for axis labels
  percentVar <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)

  # Create 3D PCA plot using plotly
  pca_plot <- plot_ly(
    pca_data,
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~condition, 
    symbol = ~organ, 
    symbols = c("Lung" = "triangle-down", "Spleen" = "square"),  # Custom shapes
    text = ~paste("Sample:", sample_id, "<br>Group:", cell_genotype),  
    marker = list(size = 8)
  ) %>%
    layout(
      title = title,
      scene = list(
        xaxis = list(title = paste0("PC1: ", percentVar[1], "% variance")),
        yaxis = list(title = paste0("PC2: ", percentVar[2], "% variance")),
        zaxis = list(title = paste0("PC3: ", percentVar[3], "% variance"))
      ),
      legend = list(x = 1, y = 0.9)
    )

  return(pca_plot)
}
