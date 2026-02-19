# -------------------------------------------------------------------------
#  build a NES-matrix from the list returned by run_gsea_analysis()
#  rows   = pathways (top n_pathways per database, union of all dbs)
#  cols   = databases
#  values = NES (set to 0 if pathway not present / not significant)
# -------------------------------------------------------------------------
gsea_to_matrix <- function(gsea_results, n_pathways = 30, padj_cutoff = 0.05) {
  # Ensure padj_cutoff is numeric
  padj_cutoff_num <- as.numeric(padj_cutoff)
  if (is.na(padj_cutoff_num)) {
    warning("padj_cutoff is not numeric, using default value of 0.05")
    padj_cutoff_num <- 0.05
  }

  nes_list <- lapply(gsea_results, function(res) {
      if (is.null(res) || nrow(res@result) == 0) return(NULL)
      res@result |>
        dplyr::filter(p.adjust < padj_cutoff_num) |>
        dplyr::arrange(p.adjust)              |>
        utils::head(n_pathways)               |>
        dplyr::select(ID, NES) |>
        tibble::column_to_rownames("ID")
  })

  # collect all unique pathway IDs that survived any DB
  all_ids <- unique(unlist(lapply(nes_list, rownames)))

  nes_mat <- sapply(nes_list, function(tb) {
                vec <- numeric(length(all_ids)); names(vec) <- all_ids
                if (!is.null(tb)) vec[rownames(tb)] <- tb$NES
                vec
             })
  # give nicer col-names (same as list names)
  colnames(nes_mat) <- names(nes_list)
  return(nes_mat)
}

# -------------------------------------------------------------------------
#  draw and save the heat-map
# -------------------------------------------------------------------------
gsea_heatmap_save <- function(nes_matrix,
                              file,
                              annotation_col = NULL,
                              ann_colors     = NULL,
                              main           = "Pathway NES (top sets)",
                              gaps_col       = NULL,
                              cluster_cols   = FALSE) {

  if (!requireNamespace("pheatmap", quietly = TRUE))
        stop("Please install the 'pheatmap' package")

  # Check if the matrix is empty
  if (nrow(nes_matrix) == 0 || ncol(nes_matrix) == 0) {
    warning("Empty matrix provided to gsea_heatmap_save. No heatmap will be generated.")
    return(NULL)
  }

  # order pathways by max |NES| so the strongest ones sit on top
  ord <- order(apply(abs(nes_matrix), 1, max), decreasing = TRUE)
  nes_matrix <- nes_matrix[ord, , drop = FALSE]

  ## heat-map -------------------------------------------------------------
  tryCatch({
    ph <- pheatmap::pheatmap(
            mat                = nes_matrix,
            color              = colorRampPalette(
                                    c("#2166AC","#F7F7F7","#B35806"))(100),
            scale              = "none",
            clustering_method  = "complete",
            cluster_rows       = TRUE,
            cluster_cols       = cluster_cols,
            annotation_col     = annotation_col,
            annotation_colors  = ann_colors,
            show_rownames      = TRUE,
            show_colnames      = TRUE,
            gaps_col           = gaps_col,
            border_color       = NA,
            main               = main,
            fontsize_row       = 6)

    ## write to pdf ---------------------------------------------------------
    pdf(file, width = 8, height = max(6, nrow(nes_matrix) * 0.15))
    print(ph)
    dev.off()

    message("   ↳ heat-map saved to ", normalizePath(file))
  }, error = function(e) {
    warning("Error generating heatmap: ", e$message)
    return(NULL)
  })
}
