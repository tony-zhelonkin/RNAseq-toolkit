#' Plot all GSEA results for a given analysis
#'
#' @param gsea_list List of GSEA results
#' @param analysis_name Name of the analysis
#' @param out_root Output root directory
#' @param n_pathways Number of pathways to show
#' @param padj_cutoff Adjusted p-value cutoff
#' @param expr_data Expression data matrix for heatmaps (optional)
#' @param sample_annotation Sample annotation data.frame (optional)
#' @param sample_order Order of samples for plots (optional)
#' @param ann_colors Colors for annotations (optional)
#'
#' @return NULL (invisibly)
#' @export
plot_all_gsea_results <- function(
    gsea_list,
    analysis_name,
    out_root,
    n_pathways = 20,
    padj_cutoff = 0.05,
    expr_data = NULL,
    sample_annotation = NULL,
    sample_order = NULL,
    ann_colors = NULL
) {
  if (is.null(gsea_list) || length(gsea_list) == 0) {
    message("No GSEA results found for ", analysis_name)
    return(invisible(NULL))
  }
  
  # Create output directory
  out_dir <- file.path(out_root, analysis_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process each database result
  for (db_name in names(gsea_list)) {
    gsea_obj <- gsea_list[[db_name]]
    
    # Skip if no results or empty
    if (is.null(gsea_obj) || length(gsea_obj@result) == 0) {
      message("Skipping ", db_name, " (no results)")
      next
    }
    
    # Create database-specific directory
    db_dir <- file.path(out_dir, db_name)
    dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get plot parameters based on database type
    plot_params <- get_db_plot_params(db_name)
    
    # Extract results
    results <- as.data.frame(gsea_obj@result)
    
    # If no significant results, create an empty plot with message
    if (nrow(results[results$p.adjust < padj_cutoff, ]) == 0) {
      message("No significant pathways for ", db_name, " in ", analysis_name)
      # Create an empty plot
      p_empty <- ggplot2::ggplot() + 
        ggplot2::labs(title = paste0(db_name, ": No significant pathways found"))
      
      save_gsea_plot(
        p_empty,
        filename = file.path(db_dir, paste0(analysis_name, "_", db_name, "_no_results.pdf")),
        width = plot_params$width,
        height = plot_params$height,
        base_font_size = plot_params$font_size
      )
      next
    }
    # 1. Plot upregulated (positive NES) pathways
    p_up <- gsea_dotplot(
      gsea_obj = gsea_obj,
      filterBy = "NES_positive",
      showCategory = n_pathways,
      padj_cutoff = padj_cutoff,
      title = paste0(db_name, ": Upregulated in ", analysis_name)
    )
    
    save_gsea_plot(
      p_up,
      filename = paste0(analysis_name, "_", db_name, "_up_dot.pdf"),
      width = plot_params$width,
      height = plot_params$height,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 2. Plot downregulated (negative NES) pathways
    p_down <- gsea_dotplot(
      gsea_obj = gsea_obj,
      filterBy = "NES_negative",
      showCategory = n_pathways,
      padj_cutoff = padj_cutoff,
      title = paste0(db_name, ": Downregulated in ", analysis_name)
    )
    
    save_gsea_plot(
      p_down,
      filename = paste0(analysis_name, "_", db_name, "_down_dot.pdf"),
      width = plot_params$width,
      height = plot_params$height,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 3. Plot faceted dotplot (added to match run_gsea_analysis.R)
    p_facet <- gsea_dotplot_facet(
      gsea_obj = gsea_obj,
      showCategory = n_pathways,
      padj_cutoff = padj_cutoff,
      title = paste0(db_name, ": ", analysis_name)
    )
    
    save_gsea_plot(
      p_facet,
      filename = paste0(analysis_name, "_", db_name, "_facet.pdf"),
      width = plot_params$width,
      height = plot_params$height * 1.4,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 4. Plot NES barplot (added to match run_gsea_analysis.R)
    p_bar <- gsea_barplot(
      gsea_obj = gsea_obj,
      top_n = n_pathways,
      padj_cutoff = padj_cutoff,
      title = paste0(db_name, ": ", analysis_name, " NES")
    )
    
    save_gsea_plot(
      p_bar,
      filename = paste0(analysis_name, "_", db_name, "_nes_bar.pdf"),
      width = plot_params$width,
      height = plot_params$height,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 5. Plot running sum plot for top 5 pathways (added to match run_gsea_analysis.R)
    top5 <- order(abs(gsea_obj@result$NES), decreasing = TRUE)[1:5]
    p_running <- gsea_running_sum_plot(
      gsea_obj = gsea_obj,
      gene_set_ids = top5,
      base_size = plot_params$font_size
    )
    
    save_gsea_plot(
      p_running,
      filename = paste0(analysis_name, "_", db_name, "_running_sum.pdf"),
      width = plot_params$width,
      height = plot_params$height * 1.2,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 6. Plot combined results (keeping this one as it was in the original)
    p_combined <- gsea_dotplot(
      gsea_obj = gsea_obj,
      filterBy = "p.adjust",
      showCategory = n_pathways,
      padj_cutoff = padj_cutoff,
      title = paste0(db_name, ": All significant pathways in ", analysis_name)
    )
    
    save_gsea_plot(
      p_combined,
      filename = paste0(analysis_name, "_", db_name, "_combined.pdf"),
      width = plot_params$width,
      height = plot_params$height,
      base_font_size = plot_params$font_size,
      dir = db_dir
    )
    
    # 8. Save a text log of the results
    save_gsea_log(
      gsea_obj = gsea_obj,
      filename = paste0(analysis_name, "_", db_name, "_results.txt"),
      padj_cutoff = padj_cutoff,
      dir = db_dir
    )
  }
  
  return(invisible(NULL))
}
