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
  message("[plot_all_gsea_results] Starting for: ", analysis_name)

  if (is.null(gsea_list) || length(gsea_list) == 0) {
    message("[plot_all_gsea_results] No GSEA results found for ", analysis_name)
    return(invisible(NULL))
  }

  message("[plot_all_gsea_results] Processing ", length(gsea_list), " databases: ",
          paste(names(gsea_list), collapse = ", "))

  # Create output directory
  out_dir <- file.path(out_root, analysis_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("[plot_all_gsea_results] Output directory: ", out_dir)

  # Process each database result
  for (db_name in names(gsea_list)) {
    message("[plot_all_gsea_results] Processing database: ", db_name)
    gsea_obj <- gsea_list[[db_name]]

    # Skip if no results or empty
    if (is.null(gsea_obj)) {
      message("[plot_all_gsea_results] Skipping ", db_name, " (NULL object)")
      next
    }

    if (!methods::is(gsea_obj, "gseaResult")) {
      message("[plot_all_gsea_results] Skipping ", db_name,
              " (not a gseaResult object, class: ", class(gsea_obj)[1], ")")
      next
    }

    if (length(gsea_obj@result) == 0 || nrow(gsea_obj@result) == 0) {
      message("[plot_all_gsea_results] Skipping ", db_name,
              " (empty results: ", nrow(gsea_obj@result), " pathways)")
      next
    }

    message("[plot_all_gsea_results] ", db_name, " has ", nrow(gsea_obj@result), " pathways")
    
    # Create database-specific directory
    db_dir <- file.path(out_dir, db_name)
    dir.create(db_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get plot parameters based on database type
    plot_params <- get_db_plot_params(db_name)
    
    # Extract results
    results <- as.data.frame(gsea_obj@result)
    n_sig <- sum(results$p.adjust < padj_cutoff, na.rm = TRUE)

    # If no significant results, create an empty plot with message
    if (n_sig == 0) {
      message("[plot_all_gsea_results] No significant pathways for ", db_name,
              " in ", analysis_name, " (padj < ", padj_cutoff, ")")
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

    message("[plot_all_gsea_results] Generating plots for ", db_name, " (", n_sig, " significant pathways)")

    # 1. Plot upregulated (positive NES) pathways
    message("[plot_all_gsea_results]   Creating upregulated dotplot...")
    tryCatch({
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
      message("[plot_all_gsea_results]   ✓ Saved upregulated dotplot")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating upregulated dotplot: ", e$message)
    })
    
    # 2. Plot downregulated (negative NES) pathways
    message("[plot_all_gsea_results]   Creating downregulated dotplot...")
    tryCatch({
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
      message("[plot_all_gsea_results]   ✓ Saved downregulated dotplot")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating downregulated dotplot: ", e$message)
    })

    # 3. Plot faceted dotplot
    message("[plot_all_gsea_results]   Creating faceted dotplot...")
    tryCatch({
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
      message("[plot_all_gsea_results]   ✓ Saved faceted dotplot")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating faceted dotplot: ", e$message)
    })

    # 4. Plot NES barplot
    message("[plot_all_gsea_results]   Creating NES barplot...")
    tryCatch({
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
      message("[plot_all_gsea_results]   ✓ Saved NES barplot")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating NES barplot: ", e$message)
    })

    # 5. Plot running sum plot for top 5 pathways
    message("[plot_all_gsea_results]   Creating running sum plots...")
    tryCatch({
      top5 <- order(abs(gsea_obj@result$NES), decreasing = TRUE)[1:min(5, nrow(gsea_obj@result))]
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
      message("[plot_all_gsea_results]   ✓ Saved running sum plots")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating running sum plots: ", e$message)
    })

    # 6. Plot combined results
    message("[plot_all_gsea_results]   Creating combined dotplot...")
    tryCatch({
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
      message("[plot_all_gsea_results]   ✓ Saved combined dotplot")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error creating combined dotplot: ", e$message)
    })

    # 8. Save a text log of the results
    message("[plot_all_gsea_results]   Saving results log...")
    tryCatch({
      save_gsea_log(
        gsea_obj = gsea_obj,
        filename = paste0(analysis_name, "_", db_name, "_results.txt"),
        padj_cutoff = padj_cutoff,
        dir = db_dir
      )
      message("[plot_all_gsea_results]   ✓ Saved results log")
    }, error = function(e) {
      message("[plot_all_gsea_results]   ✗ Error saving results log: ", e$message)
    })

    message("[plot_all_gsea_results] Completed ", db_name)
  }

  message("[plot_all_gsea_results] Finished processing all databases for ", analysis_name)
  return(invisible(NULL))
}
