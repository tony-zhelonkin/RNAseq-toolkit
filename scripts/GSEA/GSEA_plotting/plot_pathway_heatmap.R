# Function to create heatmaps for GSEA results with improved error handling
create_gsea_heatmap <- function(
  gsea_scores_matrix,
  annotation_col,
  ann_colors,
  title = "GSEA Pathway Scores",
  sample_order = NULL,
  db_type = NULL,
  gaps_col = NULL
) {
  # Debug: Print information about input data
  message(sprintf("[DEBUG] create_gsea_heatmap: input matrix has %d samples and %d pathways", 
                 nrow(gsea_scores_matrix), ncol(gsea_scores_matrix)))
  message(sprintf("[DEBUG] First few sample names: %s", 
                 paste(head(rownames(gsea_scores_matrix)), collapse=", ")))
  
  # 1. Make sure all samples from annotation_col are included in the matrix
  all_samples <- rownames(annotation_col)
  missing_samples <- all_samples[!all_samples %in% rownames(gsea_scores_matrix)]
  
  if (length(missing_samples) > 0) {
    message(sprintf("[WARNING] Adding %d missing samples: %s", 
                   length(missing_samples), 
                   paste(missing_samples, collapse=", ")))
    
    empty_rows <- matrix(NA, nrow = length(missing_samples), ncol = ncol(gsea_scores_matrix))
    rownames(empty_rows) <- missing_samples
    colnames(empty_rows) <- colnames(gsea_scores_matrix)
    gsea_scores_matrix <- rbind(gsea_scores_matrix, empty_rows)
  } else {
    message("[INFO] No missing samples detected.")
  }
  
  # 2. If sample_order is provided, reorder the matrix
  if (!is.null(sample_order)) {
    # Check if all samples in sample_order exist in the matrix
    missing_ordered_samples <- sample_order[!sample_order %in% rownames(gsea_scores_matrix)]
    if (length(missing_ordered_samples) > 0) {
      message("[WARNING] Some samples in sample_order are not in the matrix: ", 
             paste(missing_ordered_samples, collapse=", "))
    }
    
    # Keep only samples that exist in both sample_order and matrix
    valid_order <- sample_order[sample_order %in% rownames(gsea_scores_matrix)]
    message(sprintf("[INFO] Reordering matrix using %d samples", length(valid_order)))
    
    # Reorder matrix
    gsea_scores_matrix <- gsea_scores_matrix[valid_order, , drop = FALSE]
  }
  
  # 3. Clean pathway names based on database type
  if (!is.null(db_type)) {
    clean_names <- colnames(gsea_scores_matrix)
    
    if (db_type == "hallmark") {
        clean_names <- gsub("^HALLMARK_", "", clean_names)
    } else if (db_type == "gobp") {
        clean_names <- gsub("^GOBP_", "", clean_names)
    } else if (db_type == "gomf") {
        clean_names <- gsub("^GOMF_", "", clean_names)
    } else if (db_type == "gocc") {
        clean_names <- gsub("^GOCC_", "", clean_names)
    } else if (db_type == "kegg") {
        clean_names <- gsub("^KEGG_MEDICUS_", "", clean_names)
    } else if (db_type == "reactome") {
        clean_names <- gsub("^REACTOME_", "", clean_names)
    } else if (db_type == "biocarta") {
        clean_names <- gsub("^BIOCARTA_", "", clean_names)
    } else if (db_type == "grtd") {
        clean_names <- gsub("^GTRD_", "", clean_names)
    } else if (db_type == "wiki") {
        clean_names <- gsub("^WP_", "", clean_names)
    } else if (db_type == "perturb") {
        clean_names <- gsub("^CGP_", "", clean_names)
    }
    
    # Replace underscores with spaces for better readability
    clean_names <- gsub("_", " ", clean_names)
    
    colnames(gsea_scores_matrix) <- clean_names
  }
  
  # 4. Determine plot width based on database type
  plot_width <- 10  # Default width
  if (!is.null(db_type) && (db_type == "kegg" || db_type == "reactome" || 
                            db_type == "grtd" || db_type == "biocarta" ||
                            db_type == "gobp")) {
      plot_width <- 14  # Wider plot for databases with long pathway names
  }
  
  # Debug: Final matrix dimensions
  message(sprintf("[DEBUG] Final matrix dimensions: %d samples x %d pathways", 
                 nrow(gsea_scores_matrix), ncol(gsea_scores_matrix)))
  
  # Make sure annotation_col matches the samples in the matrix
  annotation_col_subset <- annotation_col[rownames(gsea_scores_matrix), , drop=FALSE]
  
  # Create the heatmap
  p <- pheatmap(
    t(gsea_scores_matrix),  # Transpose so pathways are rows, samples are columns
    scale = "row",          # Scale by row (pathway)
    clustering_method = "complete",
    clustering_distance_rows = "correlation",
    cluster_cols = FALSE,   # Don't cluster samples to keep your order
    show_rownames = TRUE,   # Show pathway names
    show_colnames = FALSE,  # Don't show sample names
    annotation_col = annotation_col_subset,  # Match annotations to matrix
    annotation_colors = ann_colors,
    color = viridis::viridis(100),  # Use viridis color palette
    main = title,
    gaps_col = gaps_col,       # Add visual separator between Th1 and Th17
    fontsize_row = 8        # Adjust font size for pathway names
  )
  
  # Return both the pheatmap object and the plot width for saving
  return(list(
    plot = p,
    width = plot_width
  ))
}

# Function to visualize all GSEA database results with improved error handling
visualize_gsea_results <- function(
  gsea_results,
  annotation_col,
  ann_colors,
  sample_order = NULL,
  output_dir = "results/gsea_heatmaps",
  gaps_col = NULL,
  min_pathways = 1  # Minimum number of pathways required to create a heatmap
) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Updated list of databases to visualize - match exactly with what's in the GSEA results
  databases <- c("hallmark", "wiki", "kegg", "perturb", "grtd", "gobp", "gomf", "gocc", "reactome", "biocarta")
  
  # Create a heatmap for each database
  heatmap_list <- list()
  
  # Set up logging
  log_file <- file.path(output_dir, "heatmap_creation_log.txt")
  log_conn <- file(log_file, open = "w")
  on.exit(close(log_conn))
  
  # Function to log messages
  log_message <- function(msg, type = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full_msg <- sprintf("[%s] [%s] %s", timestamp, type, msg)
    message(full_msg)
    writeLines(full_msg, log_conn)
  }
  
  log_message("Starting GSEA heatmap visualization")
  log_message(sprintf("Output directory: %s", output_dir))
  
  # Check if annotation_col has rownames
  if (is.null(rownames(annotation_col))) {
    log_message("annotation_col does not have rownames. This will cause issues.", "ERROR")
    return(NULL)
  } else {
    log_message(sprintf("annotation_col has %d samples: %s", 
                       nrow(annotation_col),
                       paste(head(rownames(annotation_col), 5), collapse=", ")))
  }
  
  # Process each database
  for (db in databases) {
    tryCatch({
      # Check if database exists in results
      if (is.null(gsea_results$scores[[db]])) {
        log_message(sprintf("Database '%s' not found in gsea_results$scores", db), "WARNING")
        next
      }
      
      # Check if there are any pathways
      if (ncol(gsea_results$scores[[db]]) < min_pathways) {
        log_message(sprintf("Database '%s' has only %d pathways (minimum %d required)", 
                           db, ncol(gsea_results$scores[[db]]), min_pathways), "WARNING")
        next
      }
      
      log_message(sprintf("\nCreating heatmap for %s database with %d pathways", 
                         toupper(db), ncol(gsea_results$scores[[db]])))
      
      # Create title
      title <- paste0("GSEA ", toupper(db), " Pathway Scores")
      
      # Create the heatmap with database type for name cleaning
      heatmap_result <- create_gsea_heatmap(
        gsea_scores_matrix = gsea_results$scores[[db]],
        annotation_col = annotation_col,
        ann_colors = ann_colors,
        title = title,
        sample_order = sample_order,
        db_type = db,
        gaps_col = gaps_col
      )
      
      # Save the heatmap to a PDF file with appropriate width
      pdf_file <- file.path(output_dir, paste0("gsea_", db, "_heatmap.pdf"))
      log_message(sprintf("Saving heatmap to %s", pdf_file))
      
      # Calculate appropriate height based on number of pathways
      height <- 8 + ncol(gsea_results$scores[[db]]) * 0.15
      
      # Create PDF
      pdf(pdf_file, width = heatmap_result$width, height = height)
      print(heatmap_result$plot)
      dev.off()
      
      # Store the heatmap object
      heatmap_list[[db]] <- heatmap_result$plot
      
      log_message(sprintf("Successfully created heatmap for %s database", toupper(db)))
      
    }, error = function(e) {
      log_message(sprintf("Error creating heatmap for %s database: %s", 
                         toupper(db), e$message), "ERROR")
    })
  }
  
  log_message(sprintf("Completed heatmap creation. Created %d heatmaps out of %d databases.", 
                     length(heatmap_list), length(databases)))
  
  return(heatmap_list)
}