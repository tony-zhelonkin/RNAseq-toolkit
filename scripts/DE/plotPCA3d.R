#' Create an Interactive 3D PCA Plot from a DGEList Object
#'
#' Performs PCA on logCPM values from an edgeR DGEList object and generates
#' an interactive 3D PCA plot using plotly.
#' Assumes `DGE_object$samples` contains columns 'group' (for color/text) and
#' 'organ' (for symbol). Uses rownames as sample IDs in hover text.
#' Uses hardcoded plotly symbols for 'Lung' ('triangle-down') and 'Spleen' ('square').
#'
#' @param DGE_object A `DGEList` object (from edgeR) containing counts and sample metadata.
#'        Must have `DGE_object$samples` with columns 'group' and 'organ'. Rownames are used as sample IDs.
#' @param title Character, plot title (default: "3D PCA Plot").
#' @param point_size Numeric, size of the markers in the plot (default: 8).
#'
#' @return A plotly object representing the interactive 3D PCA plot.
#' @export
#' @import edgeR stats plotly
#' @importFrom magrittr %>%
#'
#' @examples
#' # Assuming dge_list is a valid DGEList object prepared with edgeR
#' # plot_3d <- create_3d_pca_plot(dge_list, title = "3D PCA of Experiment")
#' # plot_3d # To view the interactive plot
library(edgeR)
library(stats)
library(plotly)
library(magrittr) # For %>% pipe

create_3d_pca_plot <- function(DGE_object, title = "3D PCA Plot", point_size = 8) {

  # Input validation (basic)
  if (!inherits(DGE_object, "DGEList")) {
    stop("Input must be a DGEList object from the edgeR package.")
  }
   if (!all(c("group", "organ") %in% colnames(DGE_object$samples))) {
     warning("Expected columns 'group' and 'organ' not found in DGE_object$samples. Plotting may fail or look incorrect.")
     # Add placeholder columns if missing
     if (!"group" %in% colnames(DGE_object$samples)) DGE_object$samples$group <- "Unknown"
     if (!"organ" %in% colnames(DGE_object$samples)) DGE_object$samples$organ <- "Unknown"
  }

  # Calculate logCPM
  logCPM <- edgeR::cpm(DGE_object, log = TRUE, prior.count = 1)

  # Perform PCA
  # Check for zero-variance genes/samples before PCA
  logCPM_filtered <- logCPM[apply(logCPM, 1, var) > 0, apply(logCPM, 2, var) > 0]
   if(nrow(logCPM_filtered) < 3 || ncol(logCPM_filtered) < 3) {
      stop("Insufficient data variation for 3D PCA after filtering zero-variance rows/columns.")
  }
  pca <- stats::prcomp(t(logCPM_filtered))

  # Ensure at least 3 PCs were computed
  if (ncol(pca$x) < 3) {
      stop("PCA resulted in fewer than 3 principal components. Cannot create 3D plot.")
  }

  # Prepare PCA data frame, matching metadata to PCA'd samples
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    # Ensure metadata matches the samples included in PCA
    DGE_object$samples[colnames(logCPM_filtered), , drop = FALSE]
  )
  # Add sample_id from rownames and ensure group is character for labels
  pca_data$sample_id <- rownames(pca_data)
  pca_data$label_text <- as.character(pca_data$group)


  # Define symbol mapping (plotly symbols)
  symbol_mapping <- c("Lung" = "triangle-down", "Spleen" = "square")
  # Add a default symbol for any other organ types found
  other_organs <- setdiff(unique(pca_data$organ), names(symbol_mapping))
  if(length(other_organs) > 0) {
      symbol_mapping <- c(symbol_mapping, setNames(rep("circle", length(other_organs)), other_organs))
      warning("Found unexpected organ types: ", paste(other_organs, collapse=", "), ". Assigning default symbol ('circle').")
  }

  # Compute variance percentages for axis labels
  percentVar <- round(100 * (pca$sdev^2) / sum(pca$sdev^2), 1)

  # Create 3D PCA plot using plotly
  pca_plot <- plot_ly(
    data = pca_data, # Explicitly pass data
    x = ~PC1, y = ~PC2, z = ~PC3,
    color = ~group, # Use 'group' directly
    symbol = ~organ,
    symbols = symbol_mapping, # Use the defined mapping
    type = 'scatter3d', # Specify type
    mode = 'markers', # Specify mode
    text = ~paste("Sample:", sample_id, "<br>Group:", label_text), # Use defined columns
    hoverinfo = 'text', # Show only text on hover
    marker = list(size = point_size) # Use parameter for size
  ) %>%
    plotly::layout(
      title = list(text = title, x = 0.5, xanchor = 'center'), # Center title
      scene = list(
        xaxis = list(title = sprintf("PC1: %.1f%% variance", percentVar[1])),
        yaxis = list(title = sprintf("PC2: %.1f%% variance", percentVar[2])),
        zaxis = list(title = sprintf("PC3: %.1f%% variance", percentVar[3]))
      ),
      legend = list(
          orientation = "v", # Vertical legend
          x = 1.05, y = 0.9, # Position outside plot area
          title = list(text = "Group") # Add legend title for color
          # Note: Plotly doesn't easily combine symbol legends with color legends automatically
      )
    )

  # Add a separate legend trace for symbols if needed (more complex)
  # For simplicity, relying on hover text and potentially manual annotation if needed.

  return(pca_plot)
}
