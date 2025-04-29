#' Create a Heatmap of Pathway Scores
#'
#' This function creates a heatmap visualization of pathway scores (average expression
#' of pathway genes) across samples, with customizable annotations and formatting.
#'
#' @param scores_matrix Matrix of pathway scores, with samples as rows and pathways as columns.
#' @param sample_order Character vector specifying the order of samples in the heatmap.
#' @param annotation_col Data frame of sample annotations for column annotation bars.
#' @param annotation_colors List of colors for annotation bars.
#' @param title Character, title for the heatmap (default: "Pathway Scores Heatmap").
#' @param gaps_col Numeric vector of column indices where to place gaps (default: NULL).
#' @param db_type Character, database type for pathway name cleaning (default: NULL).
#'        Options: "hallmark", "gobp", "gomf", "gocc", "kegg", "reactome", "biocarta", etc.
#' @param scale Character, scaling parameter passed to pheatmap: "row", "column", or "none" (default: "row").
#' @param color_palette Function, color palette to use (default: colorRampPalette(c("navy", "white", "red"))(50)).
#' @param save_plot Logical, whether to save the plot to a file (default: FALSE).
#' @param output_dir Character, directory to save the plot (default: "plots/").
#' @param width Numeric, width of the saved plot in inches (default: 10).
#' @param height Numeric, height of the saved plot in inches (default: 8).
#'
#' @return The pheatmap object (invisibly).
#' @export
#'
#' @examples
#' # Basic usage
#' gsea_scores_heatmap(pathway_scores, sample_order, annotation_col, annotation_colors)
#'
#' # With database-specific formatting
#' gsea_scores_heatmap(pathway_scores, sample_order, annotation_col, annotation_colors,
#'                    title = "KEGG Pathway Scores", db_type = "kegg")
#'
#' # Save the plot
#' gsea_scores_heatmap(pathway_scores, sample_order, annotation_col, annotation_colors,
#'                    save_plot = TRUE, output_dir = "results/gsea/")
library(pheatmap)
library(dplyr)
library(stringr)

gsea_scores_heatmap <- function(scores_matrix,
                               sample_order = NULL,
                               annotation_col,
                               annotation_colors,
                               title = "Pathway Scores Heatmap",
                               gaps_col = NULL,
                               db_type = NULL,
                               scale = "row",
                               color_palette = colorRampPalette(c("navy", "white", "red"))(50),
                               save_plot = FALSE,
                               output_dir = "plots/",
                               width = 10,
                               height = 8) {
  
  # Validate inputs
  if (!is.matrix(scores_matrix) && !is.data.frame(scores_matrix)) {
    stop("scores_matrix must be a matrix or data frame")
  }
  
  if (is.null(rownames(scores_matrix))) {
    stop("scores_matrix must have rownames (sample IDs)")
  }
  
  if (is.null(colnames(scores_matrix))) {
    stop("scores_matrix must have colnames (pathway IDs)")
  }
  
  # Make a copy of the scores matrix to avoid modifying the original
  plot_matrix <- scores_matrix
  
  # Reorder samples if sample_order is provided
  if (!is.null(sample_order)) {
    # Check if all samples in sample_order exist in the matrix
    missing_samples <- sample_order[!sample_order %in% rownames(plot_matrix)]
    if (length(missing_samples) > 0) {
      warning("Some samples in sample_order are not in the matrix: ", 
             paste(missing_samples, collapse=", "))
    }
    
    # Keep only samples that exist in both sample_order and matrix
    valid_order <- sample_order[sample_order %in% rownames(plot_matrix)]
    if (length(valid_order) > 0) {
      plot_matrix <- plot_matrix[valid_order, , drop = FALSE]
    } else {
      warning("No valid samples in sample_order. Using original order.")
    }
  }
  
  # Clean pathway names based on database type
  if (!is.null(db_type)) {
    clean_names <- colnames(plot_matrix)
    
    # Database-specific prefix removal
    prefixes <- list(
      hallmark = "^HALLMARK_",
      gobp = "^GOBP_",
      gomf = "^GOMF_",
      gocc = "^GOCC_",
      kegg = "^KEGG_",
      reactome = "^REACTOME_",
      biocarta = "^BIOCARTA_",
      wiki = "^WP_",
      grtd = "^GTRD_",
      perturb = "^CGP_"
    )
    
    if (db_type %in% names(prefixes)) {
      clean_names <- gsub(prefixes[[db_type]], "", clean_names)
    }
    
    # Replace underscores with spaces for better readability
    clean_names <- gsub("_", " ", clean_names)
    
    colnames(plot_matrix) <- clean_names
  }
  
  # Ensure annotation_col matches the samples in the matrix
  annotation_col_subset <- annotation_col[rownames(plot_matrix), , drop=FALSE]
  
  # Create the heatmap
  p <- pheatmap(
    t(plot_matrix),  # Transpose so pathways are rows, samples are columns
    scale = scale,
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    cluster_cols = FALSE,  # Don't cluster samples to preserve order
    show_rownames = TRUE,
    show_colnames = FALSE,
    annotation_col = annotation_col_subset,
    annotation_colors = annotation_colors,
    color = color_palette,
    main = title,
    gaps_col = gaps_col,
    fontsize = 9,
    fontsize_row = 8
  )
  
  # Save the plot if requested
  if (save_plot) {
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create filename from title
    filename <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9_]+", "_", title), ".pdf"))
    
    # Save the plot
    pdf(filename, width = width, height = height)
    print(p)
    dev.off()
    
    message("Plot saved to: ", filename)
  }
  
  return(invisible(p))
}

#' Visualize Multiple GSEA Database Results as Heatmaps
#'
#' This function creates and saves heatmaps for multiple GSEA database results,
#' with customized formatting for each database.
#'
#' @param gsea_results A list containing GSEA results with a 'scores' component.
#'        The 'scores' should be a list with database names as keys and score matrices as values.
#' @param annotation_col Data frame of sample annotations for column annotation bars.
#' @param annotation_colors List of colors for annotation bars.
#' @param sample_order Character vector specifying the order of samples (default: NULL).
#' @param output_dir Character, directory to save the plots (default: "plots/gsea_heatmaps").
#' @param gaps_col Numeric vector of column indices where to place gaps (default: NULL).
#' @param min_pathways Integer, minimum number of pathways required to create a heatmap (default: 1).
#'
#' @return A list of pheatmap objects, one for each database.
#' @export
#'
#' @examples
#' # Assuming pooled_gsea_results is the output from run_pooled_gsea()
#' gsea_visualize_all_databases(
#'   pooled_gsea_results,
#'   annotation_col = sample_annotations,
#'   annotation_colors = annotation_colors,
#'   sample_order = sample_order,
#'   output_dir = "results/gsea_heatmaps"
#' )
gsea_visualize_all_databases <- function(gsea_results,
                                        annotation_col,
                                        annotation_colors,
                                        sample_order = NULL,
                                        output_dir = "plots/gsea_heatmaps",
                                        gaps_col = NULL,
                                        min_pathways = 1) {
  
  # Validate input
  if (!is.list(gsea_results) || !("scores" %in% names(gsea_results))) {
    stop("gsea_results must be a list with a 'scores' component")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Database-specific plot parameters
  db_params <- list(
    hallmark = list(width = 10, height = 7),
    reactome = list(width = 16, height = 12),
    kegg = list(width = 14, height = 10),
    gobp = list(width = 14, height = 12),
    gomf = list(width = 12, height = 8),
    gocc = list(width = 12, height = 8),
    biocarta = list(width = 14, height = 10),
    wiki = list(width = 14, height = 10),
    grtd = list(width = 14, height = 10),
    perturb = list(width = 14, height = 10)
  )
  
  # List to store heatmap objects
  heatmap_list <- list()
  
  # Process each database
  for (db in names(gsea_results$scores)) {
    # Skip if no scores for this database
    if (is.null(gsea_results$scores[[db]])) {
      message(sprintf("No scores found for database '%s'", db))
      next
    }
    
    # Skip if not enough pathways
    if (ncol(gsea_results$scores[[db]]) < min_pathways) {
      message(sprintf("Database '%s' has only %d pathways (minimum %d required)",
                     db, ncol(gsea_results$scores[[db]]), min_pathways))
      next
    }
    
    message(sprintf("Creating heatmap for %s database with %d pathways",
                   toupper(db), ncol(gsea_results$scores[[db]])))
    
    # Get database-specific parameters or use defaults
    params <- if (db %in% names(db_params)) db_params[[db]] else list(width = 12, height = 8)
    
    # Calculate appropriate height based on number of pathways
    height <- params$height + (ncol(gsea_results$scores[[db]]) - 10) * 0.2
    height <- max(params$height, height)  # Ensure minimum height
    
    # Create title
    title <- paste0("GSEA ", toupper(db), " Pathway Scores")
    
    # Create the heatmap
    p <- gsea_scores_heatmap(
      scores_matrix = gsea_results$scores[[db]],
      sample_order = sample_order,
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      title = title,
      gaps_col = gaps_col,
      db_type = db,
      save_plot = TRUE,
      output_dir = output_dir,
      width = params$width,
      height = height
    )
    
    # Store the heatmap object
    heatmap_list[[db]] <- p
  }
  
  message(sprintf("Completed heatmap creation. Created %d heatmaps.", length(heatmap_list)))
  
  return(heatmap_list)
}
