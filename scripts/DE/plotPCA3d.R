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
