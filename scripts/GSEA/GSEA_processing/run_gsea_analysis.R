#' Run Comprehensive GSEA Analysis with Visualizations
#'
#' This function performs a comprehensive GSEA analysis on differential expression results,
#' generating various visualizations for different MSigDB databases.
#'
#' @param de_table A data frame containing differential expression results with at least
#'        the following columns: gene identifiers as rownames, logFC, and a ranking metric column.
#' @param analysis_name Character, name for the analysis (used in plot titles and filenames).
#' @param rank_metric Character, name of the column to use for ranking genes (default: "t").
#' @param n_pathways Integer, number of pathways to display in plots (default: 15).
#' @param q_cutoff Numeric, q-value cutoff for significance (default: 0.05).
#' @param save_plots Logical, whether to save plots to files (default: TRUE).
#' @param output_dir Character, directory to save plots (default: "results/gsea/").
#' @param databases List, database configurations to use (default: NULL, uses standard set).
#'
#' @return A list of GSEA results by database.
#' @export
#'
#' @examples
#' # Run analysis with default settings
#' gsea_results <- run_gsea_analysis(de_table, "Treatment_vs_Control")
#'
#' # Run analysis with custom databases
#' db_configs <- list(
#'   HALLMARK = list(category = "H", subcategory = NULL),
#'   KEGG = list(category = "C2", subcategory = "CP:KEGG")
#' )
#' gsea_results <- run_gsea_analysis(de_table, "Treatment_vs_Control", 
#'                                  databases = db_configs)

source("scripts/GSEA/GSEA_processing/run_gsea.R")
library(dplyr)

run_gsea_analysis <- function(de_table, 
                             analysis_name,
                             rank_metric = "t",
                             n_pathways = 15,
                             q_cutoff = 0.05,
                             save_plots = TRUE,
                             output_dir = "results/gsea/",
                             databases = NULL) {
  
  # Define default database configurations if not provided
  if (is.null(databases)) {
      databases <- list(
      hallmark = list(name = "Hallmark MH", db_species = "MM", collection = "MH", subcollection = ""),
      wiki = list(name = "Wikipathways M2", db_species = "MM", collection = "M2", subcollection = "CP:WIKIPATHWAYS"),
      kegg = list(name = "KEGG Homo sapiens C2", db_species = "HS", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
      perturb = list(name = "Perturbation (CGP) M2:CGP", db_species = "MM", collection = "M2", subcollection = "CGP"),
      grtd = list(name = "Regulatory sets TF M3:GTRD", db_species = "MM", collection = "M3", subcollection = "GTRD"),
      gobp = list(name = "GO:BP", db_species = "MM", collection = "M5", subcollection = "GO:BP"),
      gomf = list(name = "GO:MF", db_species = "MM", collection = "M5", subcollection = "GO:MF"),
      gocc = list(name = "GO:CC", db_species = "MM", collection = "M5", subcollection = "GO:CC"),
      reactome = list(name = "Reactome", db_species = "MM", collection = "M2", subcollection = "CP:REACTOME"),
      biocarta = list(name = "Biocarta", db_species = "MM", collection = "M2", subcollection = "CP:BIOCARTA")
    )
  }
  
  # Define database-specific plot parameters
  db_plot_params <- list(
        hallmark = list(width = 10, height = 7, font.size = 10),
        kegg = list(width = 14, height = 10, font.size = 9),
        gobp = list(width = 14, height = 12, font.size = 9),
        gomf = list(width = 12, height = 8, font.size = 9),
        gocc = list(width = 12, height = 8, font.size = 9),
        reactome = list(width = 16, height = 12, font.size = 8),
        biocarta = list(width = 14, height = 10, font.size = 9),
        wiki = list(width = 14, height = 10, font.size = 9)
    )
  
  # Create the output directory if it doesn't exist
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # A list to store all GSEA results by database
  result_list <- list()
  
  # Loop through each database config
  for (db_name in names(databases)) {
    # Get database configuration
    db_config <- databases[[db_name]]
    
    # Get plot parameters (use defaults if not specified)
    params <- if (db_name %in% names(db_plot_params)) {
      db_plot_params[[db_name]]
    } else {
      list(width = 12, height = 8, font.size = 9)
    }
    
    # Run GSEA
    message(paste("Running GSEA for", db_name, "database..."))
    gsea_result <- runGSEA(
      DE_results = de_table,
      rank_metric = rank_metric,
      category = db_config$category,
      subcategory = db_config$subcategory,
      padj_method = "fdr",
      nperm = 100000,
      pvalue_cutoff = 0.05
    )
    
    # Store the GSEA result in our output list
    result_list[[db_name]] <- gsea_result
    
    if (save_plots) {
      # Generate Upregulated (NES>0) dotplot
      message(paste("Creating upregulated dotplot for", db_name, "..."))
      gsea_dotplot(
        gsea_result,
        filterBy = "NES_positive",
        sortBy = "GeneRatio",
        font.size = params$font.size,
        showCategory = n_pathways,
        q_cut = q_cutoff,
        replace_ = TRUE,
        capitalize_1 = FALSE,
        capitalize_all = FALSE,
        min.dotSize = 2,
        title = paste(analysis_name, db_name, "Upregulated"),
        save_plot = TRUE,
        output_dir = output_dir,
        width = params$width,
        height = params$height
      )
      
      # Generate Downregulated (NES<0) dotplot
      message(paste("Creating downregulated dotplot for", db_name, "..."))
      gsea_dotplot(
        gsea_result,
        filterBy = "NES_negative",
        sortBy = "GeneRatio",
        font.size = params$font.size,
        showCategory = n_pathways,
        q_cut = q_cutoff,
        replace_ = TRUE,
        capitalize_1 = FALSE,
        capitalize_all = FALSE,
        min.dotSize = 2,
        title = paste(analysis_name, db_name, "Downregulated"),
        save_plot = TRUE,
        output_dir = output_dir,
        width = params$width,
        height = params$height
      )
      
      # Generate faceted dotplot
      message(paste("Creating faceted dotplot for", db_name, "..."))
      gsea_dotplot_facet(
        gsea_result,
        showCategory = n_pathways,
        font.size = params$font.size,
        title = paste(analysis_name, db_name, "Pathways"),
        q_cut = q_cutoff,
        replace_ = TRUE,
        capitalize_1 = FALSE,
        capitalize_all = FALSE,
        save_plot = TRUE,
        output_dir = output_dir,
        width = params$width,
        height = params$height
      )
      
      # Generate NES barplot
      message(paste("Creating NES barplot for", db_name, "..."))
      gsea_barplot(
        gsea_result,
        top_n = n_pathways * 2,
        title = paste(analysis_name, db_name, "NES"),
        q_cut = q_cutoff,
        replace_ = TRUE,
        capitalize_1 = FALSE,
        capitalize_all = FALSE,
        save_plot = TRUE,
        output_dir = output_dir,
        width = params$width,
        height = params$height
      )
      
      # Generate running sum plots for top pathways
      message(paste("Creating running sum plots for top", min(5, nrow(gsea_result)), "pathways in", db_name, "..."))
      top_pathways <- min(5, nrow(gsea_result))
      if (top_pathways > 0) {
        gsea_running_sum_plot(
          gsea_result,
          gene_set_ids = 1:top_pathways,
          title = paste(analysis_name, db_name, "Running Sum"),
          save_plot = TRUE,
          output_dir = output_dir,
          width = params$width,
          height = params$height / 1.5
        )
      }
    }
  }
  
  message("GSEA analysis completed successfully!")
  return(result_list)
}
