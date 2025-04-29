#' Run Comprehensive GSEA Analysis and Visualization Pipeline
#'
#' Performs GSEA using `run_gsea` for multiple specified MSigDB databases/collections
#' on a given differential expression table. Generates a standard set of plots
#' (dotplots, barplot, running sum) for each database using functions from the
#' `GSEA_plotting` directory.
#'
#' Requires `run_gsea.R` and plotting scripts (`gsea_dotplot.R`, `gsea_dotplot_facet.R`,
#' `gsea_barplot.R`, `gsea_running_sum_plot.R`) to be available in their respective
#' subdirectories within `scripts/GSEA/`.
#'
#' @param de_table A data frame with differential expression results. Must have gene
#'        identifiers as rownames, and columns for 'logFC' and the chosen `rank_metric`.
#' @param analysis_name Character, a name for this analysis run (e.g., "Treatment_vs_Control"),
#'        used in plot titles and output filenames.
#' @param rank_metric Character, the column name in `de_table` containing the metric used
#'        to pre-rank genes for GSEA (default: "t"-statistic).
#' @param species Character, the species name for `msigdbr` (e.g., "Homo sapiens",
#'        "Mus musculus"). Default: "Mus musculus".
#' @param n_pathways Integer, the maximum number of pathways to display in plots (default: 15).
#' @param padj_cutoff Numeric, adjusted p-value cutoff for significance in plots (default: 0.05).
#' @param save_plots Logical, if TRUE (default), save generated plots to `output_dir`.
#' @param output_dir Character, the directory path where plots should be saved if `save_plots` is TRUE.
#'        Defaults to "results/gsea/". Will be created if it doesn't exist.
#' @param databases List, a named list where each element defines a database to analyze.
#'        Each element should be a list with keys 'db_species', 'collection', and
#'        'subcollection' (optional, use "" or NULL if none) as used by `msigdbr` via `run_gsea`.
#'        If NULL (default), uses a standard set for the specified `species`: Hallmark, KEGG,
#'        GO (BP, MF, CC), Reactome. The names of the list elements (e.g., "HALLMARK")
#'        are used for labeling and results storage.
#' @param nperm Integer, number of permutations for GSEA (default: 100000). Passed to `run_gsea`.
#' @param pvalue_cutoff Numeric, p-value cutoff used within `run_gsea` (default: 0.05).
#'
#' @return A named list containing the `gseaResult` objects for each database analyzed.
#'         The names correspond to the keys in the `databases` list.
#' @export
#' @import dplyr
#' @importFrom methods is
#' @importFrom utils head
#'
#' @examples
#' \dontrun{
#' # Assuming de_table is a valid DE results data frame
#' # Run analysis with default settings for Mouse
#' gsea_results <- run_gsea_analysis(de_table, "MyExperiment_Default", species = "Mus musculus")
#'
#' # Run analysis for Human with specific collections
#' db_config_custom <- list(
#'   HALLMARK = list(db_species = "HS", collection = "H"), # Hallmark for Human
#'   REACTOME = list(db_species = "HS", collection = "C2", subcollection = "CP:REACTOME")
#' )
#' gsea_results_custom <- run_gsea_analysis(de_table, "MyExperiment_Custom_HS",
#'                                         species = "Homo sapiens", # This species is passed to msigdbr
#'                                         databases = db_config_custom)
#' }
library(dplyr)

# Function to safely source a script if it exists
source_safe <- function(path) {
  if (file.exists(path)) {
    source(path)
    return(TRUE)
  } else {
    warning("Required script not found: ", path)
    return(FALSE)
  }
}

run_gsea_analysis <- function(de_table,
                             analysis_name,
                             rank_metric = "t",
                             species = "Mus musculus", # Added species parameter
                             n_pathways = 15,
                             padj_cutoff = 0.05, # Renamed from q_cutoff
                             save_plots = TRUE,
                             output_dir = "results/gsea/",
                             databases = NULL,
                             nperm = 100000, # Added GSEA parameters
                             pvalue_cutoff = 0.05) { # Added GSEA parameters

  # --- Source Dependencies ---
  scripts_sourced <- all(
    source_safe("scripts/GSEA/GSEA_processing/run_gsea.R"),
    source_safe("scripts/GSEA/GSEA_plotting/gsea_dotplot.R"),
    source_safe("scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R"),
    source_safe("scripts/GSEA/GSEA_plotting/gsea_barplot.R"),
    source_safe("scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")
  )
  if (!scripts_sourced) {
    stop("One or more required scripts could not be sourced. Please check paths.")
  }
  # -------------------------

  # --- Input Validation ---
   if (!is.data.frame(de_table)) stop("`de_table` must be a data frame.")
   if (is.null(rownames(de_table))) stop("`de_table` must have rownames (gene identifiers).")
   if (!rank_metric %in% colnames(de_table)) stop(sprintf("`rank_metric` column '%s' not found in `de_table`.", rank_metric))
   if (!"logFC" %in% colnames(de_table)) warning("'logFC' column not found in `de_table`. Some plots might require it.")
  # ------------------------

  # Define default database configurations based on species if not provided
  if (is.null(databases)) {
      # Determine db_species abbreviation (simple case)
      db_species_abbr <- if (grepl("musculus", species, ignore.case = TRUE)) "MM" else "HS"
      message(sprintf("Using default databases for species '%s' (assuming db_species='%s').", species, db_species_abbr))
      databases <- list(
          HALLMARK = list(db_species = db_species_abbr, collection = if(db_species_abbr == "MM") "MH" else "H"), # MH for mouse, H for human
          KEGG = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:KEGG"),
          GO_BP = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:BP"),
          GO_MF = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:MF"),
          GO_CC = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:CC"),
          REACTOME = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:REACTOME"), 
          WIKI = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
          PERTURB = list(name = "Perturbation (CGP) M2:CGP", db_species = db_species_abbr, collection = "M2", subcollection = "CGP"),
          GRTD = list(name = "Regulatory sets TF M3:GTRD", db_species = db_species_abbr, collection = "M3", subcollection = "GTRD")
          
      )
  }

  # Define database-specific plot parameters (adjust as needed)
  # Use names from the databases list (converted to lowercase for matching)
  db_plot_params <- list(
        hallmark = list(width = 10, height = 7, font.size = 10),
        kegg = list(width = 14, height = 10, font.size = 9),
        go_bp = list(width = 14, height = 12, font.size = 9), # Match default names
        go_mf = list(width = 12, height = 8, font.size = 9),
        go_cc = list(width = 12, height = 8, font.size = 9),
        reactome = list(width = 16, height = 12, font.size = 8)
        # Add others if defaults change
        # wiki = list(width = 14, height = 10, font.size = 9)
    )
  default_plot_params <- list(width = 12, height = 8, font.size = 9)

  # Create the output directory if saving plots
  if (save_plots) {
      if (!dir.exists(output_dir)) {
          message("Creating output directory: ", output_dir)
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      if (!dir.exists(output_dir)) { # Check again if creation failed
          stop("Failed to create output directory: ", output_dir)
      }
  }
  
  # A list to store all GSEA results by database
  result_list <- list()
  
  # Loop through each database config
  for (db_name in names(databases)) {
    # Get database configuration
    db_config <- databases[[db_name]]

    # --- Validate db_config ---
    required_keys <- c("db_species", "collection")
    if (!all(required_keys %in% names(db_config))) {
        warning(sprintf("Database configuration for '%s' is missing required keys: %s. Skipping.",
                        db_name, paste(setdiff(required_keys, names(db_config)), collapse=", ")))
        result_list[[db_name]] <- NULL # Store NULL for skipped DB
        next # Skip to next database
    }
    # Ensure subcollection exists, defaulting to "" if missing or NULL
    if (!"subcollection" %in% names(db_config) || is.null(db_config$subcollection)) {
        db_config$subcollection <- ""
    }
    # --------------------------

    # Get plot parameters (use defaults if not specified)
    # Match db_name (lowercase) to plot params keys
    params <- if (tolower(db_name) %in% names(db_plot_params)) {
      db_plot_params[[tolower(db_name)]]
    } else {
      default_plot_params
    }

    # Run GSEA using run_gsea function (sourced earlier)
    message(paste("Running GSEA for", db_name, "database..."))
    gsea_result <- tryCatch({
        run_gsea( # Use the sourced function name
            DE_results = de_table,
            rank_metric = rank_metric,
            species = species, # Pass overall species name
            db_species = db_config$db_species, # Pass specific db_species from config
            collection = db_config$collection, # Pass collection from config
            subcollection = db_config$subcollection, # Pass subcollection from config
            padj_method = "fdr", # Keep consistent parameter name
            nperm = nperm,
            pvalue_cutoff = pvalue_cutoff # Pass pvalue_cutoff
            # seed is handled within run_gsea
        )
    }, error = function(e) {
        warning(sprintf("GSEA failed for database '%s': %s", db_name, e$message))
        return(NULL) # Return NULL on error
    })

    # Store the GSEA result in our output list
    result_list[[db_name]] <- gsea_result

    # Proceed with plotting only if GSEA was successful and save_plots is TRUE
    if (!is.null(gsea_result) && methods::is(gsea_result, "gseaResult") && save_plots) {
      # Generate Upregulated (NES>0) dotplot
      message(paste("Creating upregulated dotplot for", db_name, "..."))
       tryCatch({
          gsea_dotplot(
              gsea_result,
              filterBy = "NES_positive",
              sortBy = "GeneRatio",
              font.size = params$font.size,
              showCategory = n_pathways,
              q_cut = padj_cutoff, # Use renamed parameter
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
      }, error = function(e) { warning(sprintf("Upregulated dotplot failed for '%s': %s", db_name, e$message)) })

      # Generate Downregulated (NES<0) dotplot
      message(paste("Creating downregulated dotplot for", db_name, "..."))
      tryCatch({
          gsea_dotplot(
              gsea_result,
              filterBy = "NES_negative",
              sortBy = "GeneRatio",
              font.size = params$font.size,
              showCategory = n_pathways,
              q_cut = padj_cutoff, # Use renamed parameter
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
      }, error = function(e) { warning(sprintf("Downregulated dotplot failed for '%s': %s", db_name, e$message)) })

      # Generate faceted dotplot
      message(paste("Creating faceted dotplot for", db_name, "..."))
      tryCatch({
          gsea_dotplot_facet(
              gsea_result,
              showCategory = n_pathways,
              font.size = params$font.size,
              # title = paste(analysis_name, db_name, "Pathways"), # Title seems redundant with filename
              q_cut = padj_cutoff, # Use renamed parameter
              replace_ = TRUE,
              capitalize_1 = FALSE,
              capitalize_all = FALSE,
              save_plot = TRUE,
              output_dir = output_dir,
              width = params$width,
              height = params$height
          )
      }, error = function(e) { warning(sprintf("Faceted dotplot failed for '%s': %s", db_name, e$message)) })

      # Generate NES barplot
      message(paste("Creating NES barplot for", db_name, "..."))
      tryCatch({
          gsea_barplot(
              gsea_result,
              top_n = n_pathways * 2,
              # title = paste(analysis_name, db_name, "NES"), # Title seems redundant
              q_cut = padj_cutoff, # Use renamed parameter
              replace_ = TRUE,
              capitalize_1 = FALSE,
              capitalize_all = FALSE,
              save_plot = TRUE,
              output_dir = output_dir,
              width = params$width,
              height = params$height
          )
      }, error = function(e) { warning(sprintf("NES barplot failed for '%s': %s", db_name, e$message)) })

      # Generate running sum plots for top pathways
      n_running_sum <- min(5, nrow(gsea_result@result)) # Check result slot
      message(paste("Creating running sum plots for top", n_running_sum, "pathways in", db_name, "..."))
      if (n_running_sum > 0) {
        tryCatch({ # Start tryCatch
            gsea_running_sum_plot(
                gsea_result,
                gene_set_ids = 1:n_running_sum, # Use IDs from the result object
                title = paste(analysis_name, db_name, "Running Sum"),
                save_plot = TRUE,
                output_dir = output_dir,
                width = params$width,
                height = params$height / 1.5
            )
        }, error = function(e) { # Add error handling for tryCatch
             warning(sprintf("Running sum plot failed for '%s': %s", db_name, e$message))
        }) # Close tryCatch
      } else {
          message("Skipping running sum plot for ", db_name, " as no significant pathways found or result invalid.")
      } # Close if (n_running_sum > 0)
    } # Close if (!is.null(gsea_result) && ...)
  } # Close for loop
  
  message("GSEA analysis completed successfully!")
  return(result_list)
}
