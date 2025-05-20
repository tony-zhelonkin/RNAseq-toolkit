#' Create a 2D PCA Plot from a DGEList Object
#'
#' Performs PCA on logCPM values from an edgeR DGEList object and generates a ggplot2 PCA plot.
#' Assumes `DGE_object$samples` contains columns 'group' (for color/label), 'organ' (for shape),
#' and 'id' (optional, currently unused).
#' Requires `custom_minimal_theme_with_grid()` from `scripts/custom_minimal_theme.R`.
#' Uses hardcoded shape mapping for 'Lung' and 'Spleen' and a default color palette.
#'
#' @param DGE_object A `DGEList` object (from edgeR) containing counts and sample metadata.
#'        Must have `DGE_object$samples` with columns 'group' and 'organ'.
#' @param title Character, plot title (default: "PCA Plot").
#' @param xlim_abs Numeric, optional absolute limit for the x-axis (e.g., 50). If NULL (default),
#'        limits are determined automatically from data range.
#' @param ylim_abs Numeric, optional absolute limit for the y-axis (e.g., 30). If NULL (default),
#'        limits are determined automatically from data range.
#' @param point_size Numeric, size of the points in the plot (default: 5).
#' @param label_size Numeric, size of the text labels (default: 4).
#'
#' @return A ggplot2 object representing the PCA plot. The plot is NOT saved automatically.
#' @export
#' @import ggplot2 edgeR stats
#' @importFrom magrittr %>%
#'
#' @examples
#' # Assuming dge_list is a valid DGEList object prepared with edgeR
#' # plot_obj <- create_pca_plot(dge_list, title = "PCA of Experiment")
#' # To save: ggsave("my_pca.png", plot_obj)
library(ggplot2)
library(edgeR) # For cpm() and DGEList object type
library(stats) # For prcomp()
library(magrittr) # For %>%

# Source the custom theme function if it exists
# custom_theme_path <- file.path("scripts", "custom_minimal_theme.R")
# if (file.exists(custom_theme_path)) {
#   source(custom_theme_path)
# } else {
#   warning("Custom theme file not found at: ", custom_theme_path, ". Using default theme_minimal().")
#   custom_minimal_theme_with_grid <- function() theme_minimal()
# }

source_if_present <- function(...) {          # ... = path fragments
  full <- here::here(...)                    # build absolute path first
  if (!file.exists(full)) {                  # ① file there?
    warning("helper not found → ", full); return(invisible(FALSE))
  }
  old <- getwd(); on.exit(setwd(old), add = TRUE)  # ② temp-cd
  setwd(dirname(full))                             # makes relative paths resolve
  source(basename(full), echo = FALSE)
  invisible(TRUE)
}


create_pca_plot <- function(
    DGE_object, title = "PCA Plot",
    xlim_abs = NULL, ylim_abs = NULL,
    point_size = 5, label_size = 4) {

  # Input validation (basic)
  if (!inherits(DGE_object, "DGEList")) {
    stop("Input must be a DGEList object from the edgeR package.")
  }
  if (!all(c("group", "organ") %in% colnames(DGE_object$samples))) {
     warning("Expected columns 'group' and 'organ' not found in DGE_object$samples. Plotting may fail or look incorrect.")
     # Add placeholder columns if missing to avoid immediate error, though plot might be meaningless
     if (!"group" %in% colnames(DGE_object$samples)) DGE_object$samples$group <- "Unknown"
     if (!"organ" %in% colnames(DGE_object$samples)) DGE_object$samples$organ <- "Unknown"
  }
  # Add 'id' column if missing and needed (currently not used in plot)
  if (!"id" %in% colnames(DGE_object$samples)) DGE_object$samples$id <- rownames(DGE_object$samples)


  # Calculate logCPM
  logCPM <- edgeR::cpm(DGE_object, log = TRUE, prior.count = 1)

  # Perform PCA
  # Check for zero-variance genes/samples before PCA
  logCPM_filtered <- logCPM[apply(logCPM, 1, var) > 0, apply(logCPM, 2, var) > 0]
  if(nrow(logCPM_filtered) < 2 || ncol(logCPM_filtered) < 2) {
      stop("Insufficient data variation for PCA after filtering zero-variance rows/columns.")
  }
  pca <- stats::prcomp(t(logCPM_filtered))

  # Prepare PCA data frame, ensuring sample order matches PCA results
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    DGE_object$samples[colnames(logCPM_filtered), , drop = FALSE] # Match metadata to PCA'd samples
  )
  # Use 'group' column directly for labels if it exists, ensure it's character
  pca_data$label_text <- as.character(pca_data$group)

  # Define custom shapes (consider making this a parameter in future)
  # Using numbers: 6 =nabla (down triangle), 0 = square
  shape_mapping <- c("Lung" = 6, "Spleen" = 0)
  # Add a default shape for any other organ types found
  other_organs <- setdiff(unique(pca_data$organ), names(shape_mapping))
  if(length(other_organs) > 0) {
      # Assign default shape (e.g., circle = 1) to other organs
      shape_mapping <- c(shape_mapping, setNames(rep(1, length(other_organs)), other_organs))
      warning("Found unexpected organ types: ", paste(other_organs, collapse=", "), ". Assigning default shape (circle).")
  }


  # Calculate variance percentages
  percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

  # Determine axis limits
  # Use provided limits if not NULL, otherwise calculate from data range
  x_limit_val <- if (!is.null(xlim_abs)) xlim_abs else max(abs(pca_data$PC1)) * 1.1 # Add 10% buffer
  y_limit_val <- if (!is.null(ylim_abs)) ylim_abs else max(abs(pca_data$PC2)) * 1.1 # Add 10% buffer

  # Define color palette (consider making this a parameter)
  color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#000000") # Add more if needed

  # Create the plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, shape = organ)) +
    geom_point(size = point_size) +
    # Use label_text column for labels
    geom_text(aes(label = label_text), vjust = -1.5, hjust = 0.5, size = label_size, fontface = "bold", show.legend = FALSE) +
    labs(color = "Group", shape = "Organ") + # Use more generic legend titles
    scale_shape_manual(values = shape_mapping) +
    scale_color_manual(values = color_palette) + # Apply color palette
    xlab(sprintf("PC1: %.1f%% variance", percentVar[1])) +
    ylab(sprintf("PC2: %.1f%% variance", percentVar[2])) +
    ggtitle(title) +
    custom_minimal_theme_with_grid() +
    xlim(-x_limit_val, x_limit_val) + # Apply calculated/provided limits
    ylim(-y_limit_val, y_limit_val) +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.background = element_rect(color = "grey80", fill = "white", linewidth = 0.5), # Match theme
      legend.box.margin = margin(2, 2, 2, 2),
      plot.title = element_text(hjust = 0.5), # Center title
      plot.margin = margin(10, 10, 10, 10) # Adjust margins
    ) +
    coord_fixed() # Maintain aspect ratio 1:1 for PCA

  # Return the plot object
  return(pca_plot)
}
