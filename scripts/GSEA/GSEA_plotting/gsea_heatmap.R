#' Create a heatmap of GSEA results with expression z-scores
#'
#' @param gsea_obj GSEA result object
#' @param expr_data Expression data matrix (genes in rows, samples in columns)
#' @param sample_annotation Sample annotation data frame 
#' @param ann_colors List of colors for annotations
#' @param padj_cutoff Adjusted p-value cutoff
#' @param n_pathways Number of pathways to show
#' @param max_genes_per_pathway Maximum number of genes to show per pathway
#' @param sample_order Optional vector specifying sample order
#' @param title Main title for the heatmap
#' @param show_pathway_names Whether to show pathway names as row annotations
#'
#' @return List with the matrix and pheatmap object (invisibly)
#' @export
gsea_heatmap <- function(
    gsea_obj,
    expr_data,
    sample_annotation,
    ann_colors = NULL,
    padj_cutoff = 0.05,
    n_pathways = 20,
    max_genes_per_pathway = 25,
    sample_order = NULL,
    title = "GSEA Results Heatmap",
    show_pathway_names = TRUE) {
  
  # Filter for significant pathways
  top_tbl <- as.data.frame(gsea_obj@result) |>
    dplyr::filter(p.adjust < as.numeric(padj_cutoff)) |>
    dplyr::arrange(p.adjust) |>
    utils::head(n_pathways)
  
  # Check if we have any significant pathways
  if (nrow(top_tbl) == 0) {
    message("No significant pathways found with padj_cutoff = ", padj_cutoff)
    return(invisible(NULL))
  }
  
  # Check if expression data is provided
  if (missing(expr_data) || is.null(expr_data)) {
    message("Expression data not provided. Falling back to NES score heatmap.")
    return(gsea_heatmap_nes(gsea_obj, sample_annotation, ann_colors, 
                            padj_cutoff, n_pathways, sample_order, title))
  }
  
  # Clean up pathway descriptions for better labels
  top_tbl$Description <- sapply(top_tbl$Description, function(x) {
    # Strip common prefixes
    x <- gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GO_|^GOBP_|^GOCC_|^GOMF_", "", x)
    # Replace underscores with spaces
    x <- gsub("_", " ", x)
    # Convert to title case
    x <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", tolower(x), perl = TRUE)
    # Smart wrapping
    return(smart_wrap(x, 40))
  })
  
  # Extract leading edge genes for each pathway
  pathway_genes <- list()
  for (i in 1:nrow(top_tbl)) {
    # Extract genes from core_enrichment field
    genes <- unlist(strsplit(top_tbl$core_enrichment[i], "/"))
    
    # Limit number of genes per pathway
    if (length(genes) > max_genes_per_pathway) {
      genes <- genes[1:max_genes_per_pathway]
    }
    
    pathway_genes[[top_tbl$Description[i]]] <- genes
  }
  
  # Filter expression data to include only leading edge genes
  all_genes <- unique(unlist(pathway_genes))
  all_genes <- all_genes[all_genes %in% rownames(expr_data)]
  
  if (length(all_genes) == 0) {
    message("No matching genes found in expression data. Check gene identifiers.")
    return(invisible(NULL))
  }
  
  # Subset expression data
  expr_subset <- expr_data[all_genes, , drop = FALSE]
  
  # Calculate z-scores by row (gene)
  expr_z <- t(scale(t(expr_subset)))
  
  # Prepare row annotation to indicate pathway membership
  row_anno <- data.frame(row.names = all_genes)
  for (pathway in names(pathway_genes)) {
    genes <- pathway_genes[[pathway]]
    row_anno[[pathway]] <- rownames(row_anno) %in% genes
  }
  
  # Convert to numeric for better visualization
  row_anno[] <- lapply(row_anno, as.numeric)
  
  # Create row annotation colors
  pathway_colors <- list()
  for (pathway in names(pathway_genes)) {
    pathway_colors[[pathway]] <- c("0" = "white", "1" = "darkgreen")
  }
  
  # Ensure annotation color levels match the data
  valid_ann_colors <- list()
  if (!is.null(ann_colors)) {
    for (col_name in colnames(sample_annotation)) {
      if (col_name %in% names(ann_colors)) {
        # Get levels in this column of annotation
        if (is.factor(sample_annotation[[col_name]])) {
          levels_in_data <- levels(sample_annotation[[col_name]])
        } else {
          levels_in_data <- unique(as.character(sample_annotation[[col_name]]))
        }
        
        # Remove NA values
        levels_in_data <- levels_in_data[!is.na(levels_in_data)]
        
        # Filter color palette to match existing levels
        color_map <- ann_colors[[col_name]]
        matching_levels <- intersect(names(color_map), levels_in_data)
        
        if (length(matching_levels) > 0) {
          valid_ann_colors[[col_name]] <- color_map[matching_levels]
        }
      }
    }
  }
  
  # Combine all annotation colors
  all_colors <- c(valid_ann_colors, pathway_colors)
  
  # Reorder samples if requested
  if (!is.null(sample_order)) {
    valid_samples <- sample_order[sample_order %in% colnames(expr_z)]
    if (length(valid_samples) > 0) {
      expr_z <- expr_z[, valid_samples, drop = FALSE]
    } else {
      warning("None of the sample_order entries found in the matrix. Using original order.")
    }
  }
  
  # Subset the annotation to only include samples in the matrix
  valid_annotation <- sample_annotation[colnames(expr_z), , drop = FALSE]
  
  # Create the heatmap
  phm <- pheatmap::pheatmap(
    expr_z,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    annotation_col = valid_annotation,
    annotation_row = if (show_pathway_names) row_anno else NULL,
    annotation_colors = all_colors,
    main = title,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    cellwidth = 10,
    cellheight = 8,
    treeheight_row = 30,
    treeheight_col = 0,
    show_rownames = TRUE,
    show_colnames = TRUE
  )
  
  return(invisible(list(matrix = expr_z, heatmap = phm, genes = all_genes)))
}

#' Create a heatmap of GSEA results using NES values (fallback function)
#'
#' @param gsea_obj GSEA result object
#' @param sample_annotation Sample annotation data frame 
#' @param ann_colors List of colors for annotations
#' @param padj_cutoff Adjusted p-value cutoff
#' @param n_pathways Number of pathways to show
#' @param sample_order Optional vector specifying sample order
#' @param title Main title for the heatmap
#'
#' @return List with the matrix and pheatmap object (invisibly)
gsea_heatmap_nes <- function(
    gsea_obj,
    sample_annotation,
    ann_colors = NULL,
    padj_cutoff = 0.05,
    n_pathways = 20,
    sample_order = NULL,
    title = "GSEA Results Heatmap") {
  
  # Filter for significant pathways
  top_tbl <- as.data.frame(gsea_obj@result) |>
    dplyr::filter(p.adjust < as.numeric(padj_cutoff)) |>
    dplyr::arrange(p.adjust) |>
    utils::head(n_pathways)
  
  # Check if we have any significant pathways
  if (nrow(top_tbl) == 0) {
    message("No significant pathways found with padj_cutoff = ", padj_cutoff)
    return(invisible(NULL))
  }
  
  # Extract pathway IDs and NES values
  geneset <- top_tbl$ID
  nes_vec <- top_tbl$NES
  names(nes_vec) <- geneset
  
  # Get pathway descriptions for better labels
  desc_vec <- top_tbl$Description
  names(desc_vec) <- geneset
  
  # Clean up descriptions
  desc_vec <- sapply(desc_vec, function(x) {
    # Strip common prefixes
    x <- gsub("^HALLMARK_|^KEGG_|^REACTOME_|^GO_|^GOBP_|^GOCC_|^GOMF_", "", x)
    # Replace underscores with spaces
    x <- gsub("_", " ", x)
    # Convert to title case
    x <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", tolower(x), perl = TRUE)
    # Smart wrapping
    return(smart_wrap(x, 40))
  })
  
  # Verify sample IDs exist in sample_annotation
  sample_ids <- rownames(sample_annotation)
  if (length(sample_ids) == 0) {
    stop("No row names found in sample_annotation")
  }
  
  # Create the pathway × sample matrix (transposed)
  # Pathways are rows, samples are columns
  mat <- matrix(NA, 
    nrow = length(geneset), 
    ncol = length(sample_ids),
    dimnames = list(desc_vec[geneset], sample_ids)
  )
  
  # Fill the matrix with the same NES value for each pathway across all samples
  for (i in 1:nrow(mat)) {
    mat[i,] <- nes_vec[geneset[i]]
  }
  
  # Reorder samples if requested
  if (!is.null(sample_order)) {
    valid_samples <- sample_order[sample_order %in% colnames(mat)]
    if (length(valid_samples) > 0) {
      mat <- mat[, valid_samples, drop = FALSE]
    } else {
      warning("None of the sample_order entries found in the matrix. Using original order.")
    }
  }
  
  # Define color palette for heatmap
  color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)
  
  # Ensure annotation color levels match the data
  valid_ann_colors <- list()
  if (!is.null(ann_colors)) {
    for (col_name in colnames(sample_annotation)) {
      if (col_name %in% names(ann_colors)) {
        # Get levels in this column of annotation
        if (is.factor(sample_annotation[[col_name]])) {
          levels_in_data <- levels(sample_annotation[[col_name]])
        } else {
          levels_in_data <- unique(as.character(sample_annotation[[col_name]]))
        }
        
        # Remove NA values
        levels_in_data <- levels_in_data[!is.na(levels_in_data)]
        
        # Filter color palette to match existing levels
        color_map <- ann_colors[[col_name]]
        matching_levels <- intersect(names(color_map), levels_in_data)
        
        if (length(matching_levels) > 0) {
          valid_ann_colors[[col_name]] <- color_map[matching_levels]
        }
      }
    }
  }
  
  # Subset the annotation to only include samples in the matrix
  valid_annotation <- sample_annotation[colnames(mat), , drop = FALSE]
  
  # Create pheatmap
  phm <- pheatmap::pheatmap(
    mat,
    color = color_palette,
    breaks = seq(-max(abs(mat)), max(abs(mat)), length.out = 101),
    annotation_col = valid_annotation,
    annotation_colors = valid_ann_colors,
    main = title,
    cluster_rows = TRUE,    # Cluster pathways (rows)
    cluster_cols = FALSE,   # Don't cluster samples (columns)
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    cellwidth = 12,
    cellheight = 12,
    treeheight_row = 30,
    treeheight_col = 0,
    margins = c(8, 8)
  )
  
  return(invisible(list(matrix = mat, heatmap = phm)))
}
