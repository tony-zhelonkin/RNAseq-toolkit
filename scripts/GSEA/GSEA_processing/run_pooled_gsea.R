#' Run Pooled GSEA Analysis Across Multiple Contrasts and Databases
#'
#' Orchestrates GSEA across multiple contrasts using `limma` results and multiple
#' MSigDB databases via `run_gsea`. It then aggregates significant pathways across
#' contrasts for each database, extracts core enrichment genes for the top pathways,
#' and calculates pathway scores based on normalized expression data from a `DGEList` object.
#'
#' Requires helper scripts from `scripts/GSEA/GSEA_processing/` to be available.
#'
#' @param fit A `limma` fit object (e.g., from `eBayes(lmFit(...))`).
#' @param contrasts A contrast matrix compatible with `fit` (e.g., from `limma::makeContrasts`).
#' @param DGEobject An `edgeR` `DGEList` object containing normalized expression data
#'        (used for calculating pathway scores). Must contain counts accessible via `cpm()`.
#' @param species Character, the species name for `msigdbr` (e.g., "Homo sapiens",
#'        "Mus musculus"). Default: "Mus musculus". Passed to `run_gsea`.
#' @param top_n Integer, number of top pathways (ranked by min p.adjust across contrasts)
#'        to retrieve genes for and calculate scores (default: 30).
#' @param padj_cutoff Numeric, adjusted p-value cutoff used to determine significance
#'        when pooling pathways across contrasts (default: 0.05).
#' @param gsea_pvalue_cutoff Numeric, p-value cutoff passed to the underlying `run_gsea`
#'        function's `pvalue_cutoff` argument (default: 1).
#' @param rank_metric Character, the column name in `topTable` output to use for ranking
#'        genes within each contrast (default: "t"). Passed to `run_gsea`.
#' @param nperm Integer, number of permutations for GSEA (default: 100000). Passed to `run_gsea`.
#' @param databases List, optional. A named list defining databases to analyze, overriding
#'        the defaults. Structure: `list(DB_NAME = list(db_species="MM", collection="C5", subcollection="GO:BP"))`.
#' @param verbose Logical, if TRUE (default), print progress messages to console and log file.
#' @param log_file Character or NULL, path to a file for logging messages. If NULL (default),
#'        messages are only printed to the console (if `verbose` is TRUE).
#'
#' @return A list containing four named elements:
#'   \item{gsea_results}{A nested list: `list(ContrastName = list(DatabaseName = gseaResult))`.}
#'   \item{pools}{A named list: `list(DatabaseName = character_vector_of_significant_pathway_IDs)`.}
#'   \item{genes}{A named list: `list(DatabaseName = list(PathwayID = character_vector_of_genes))`.}
#'   \item{scores}{A named list: `list(DatabaseName = matrix_of_scores(samples x pathways))`.}
#' @export
#' @import limma edgeR
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # Assuming fit, contrasts, and DGEobject are correctly prepared
#' pooled_results <- run_pooled_gsea(fit, contrasts, DGEobject,
#'                                   species = "Mus musculus",
#'                                   top_n = 20,
#'                                   padj_cutoff = 0.05,
#'                                   log_file = "gsea_pooling_log.txt")
#' # Access results:
#' # head(pooled_results$pools$HALLMARK)
#' # head(pooled_results$genes$HALLMARK$SOME_PATHWAY_ID)
#' # head(pooled_results$scores$HALLMARK)
#' }

# Function to safely source a script if it exists (assuming it's defined elsewhere or copy it here)
source_safe <- function(path) {
  full_path <- file.path(getwd(), path) # Ensure path is relative to WD if needed
  if (file.exists(full_path)) {
    tryCatch({
      source(full_path, local = TRUE) # Source into function environment if possible
      return(TRUE)
    }, error = function(e) {
      warning("Error sourcing script '", full_path, "': ", e$message)
      return(FALSE)
    })
  } else {
    warning("Required script not found: ", full_path)
    return(FALSE)
  }
}


run_pooled_gsea <- function(fit,
                           contrasts,
                           DGEobject,
                           species = "Mus musculus", # Added species
                           top_n = 30,
                           padj_cutoff = 0.05, # Renamed and default changed
                           gsea_pvalue_cutoff = 1, # Added for run_gsea
                           rank_metric = "t", # Added rank_metric
                           nperm = 100000, # Added nperm
                           databases = NULL, # Added databases override
                           verbose = TRUE,
                           log_file = NULL) {

    # --- Source Dependencies ---
    # Ensure helper functions are available in the environment
    required_scripts <- c(
        "scripts/GSEA/GSEA_processing/run_gsea.R",
        "scripts/GSEA/GSEA_processing/get_significant_pathways.R",
        "scripts/GSEA/GSEA_processing/get_pathway_genes_all.R",
        "scripts/GSEA/GSEA_processing/calculate_pathway_scores.R"
    )
    scripts_sourced <- all(sapply(required_scripts, source_safe))
    if (!scripts_sourced) {
        stop("One or more required helper scripts could not be sourced. Please check paths.")
    }
    # -------------------------

    # --- Input Validation ---
    if (!inherits(fit, "MArrayLM")) stop("`fit` must be an MArrayLM object from limma.")
    if (!is.matrix(contrasts)) stop("`contrasts` must be a matrix.")
    if (!inherits(DGEobject, "DGEList")) stop("`DGEobject` must be a DGEList object from edgeR.")
    if (!all(colnames(contrasts) %in% colnames(fit$coefficients))) {
        stop("Some contrast names are not found as coefficients in the fit object.")
    }
    # ------------------------

     # Set up logging
     log_conn <- NULL
     log_file_path <- NULL
     if (!is.null(log_file)) {
         tryCatch({
             log_file_path <- normalizePath(log_file, mustWork = FALSE)
             log_conn <- file(log_file_path, open = "wt") # Use text mode
             # Setup cleanup action to close connection
             on.exit({
                 if (!is.null(log_conn) && isOpen(log_conn)) {
                     close(log_conn)
                     message("Closed log file: ", log_file_path)
                 }
             }, add = TRUE)
         }, error = function(e) {
             warning("Could not open log file '", log_file, "'. Logging to console only. Error: ", e$message)
             log_conn <<- NULL # Ensure log_conn is NULL if file opening failed
         })
     }

     # Logging function
     log_message <- function(msg, type = "INFO") {
         if (!verbose && type != "ERROR" && type != "START" && type != "END") return()
         timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
         full_msg <- sprintf("[%s] [%s] %s", timestamp, type, msg)
         message(full_msg) # Always message to console if verbose or ERROR/START/END
         if (!is.null(log_conn) && isOpen(log_conn)) {
             tryCatch({
                 writeLines(full_msg, log_conn)
                 flush(log_conn) # Ensure message is written immediately
             }, error = function(e) {
                 warning("Failed to write to log file: ", e$message)
                 # Consider setting log_conn to NULL if writing fails repeatedly
             })
         }
     }

     # Start logging
     log_message("=== Starting run_pooled_gsea ===", type = "START")
     log_message(sprintf("Parameters: species='%s', top_n=%d, padj_cutoff=%.3g, gsea_pvalue_cutoff=%.3g, rank_metric='%s', nperm=%d",
                        species, top_n, padj_cutoff, gsea_pvalue_cutoff, rank_metric, nperm))
     if (!is.null(log_file_path)) log_message(paste("Logging to file:", log_file_path))

     # Debug info
     log_message(sprintf("Input DGEobject: %d genes, %d samples.", nrow(DGEobject), ncol(DGEobject)), "DEBUG")
     log_message(sprintf("Input fit object: %d genes, %d coefficients.", nrow(fit), ncol(fit$coefficients)), "DEBUG")
     log_message(sprintf("Input contrasts: %d rows, %d columns.", nrow(contrasts), ncol(contrasts)), "DEBUG")


     # Initialize results structure
     gsea_results <- list()
     contrasts_to_analyze <- colnames(contrasts)
     log_message(paste("Contrasts to analyze:", paste(contrasts_to_analyze, collapse=", ")))

     # Define default databases if not provided by user
     if (is.null(databases)) {
         db_species_abbr <- if (grepl("musculus", species, ignore.case = TRUE)) "MM" else "HS"
         log_message(sprintf("Using default databases for species '%s' (assuming db_species='%s').", species, db_species_abbr))
         databases <- list(
             HALLMARK = list(db_species = db_species_abbr, collection = if(db_species_abbr == "MM") "MH" else "H"),
             KEGG = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:KEGG"),
             GO_BP = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:BP"),
             GO_MF = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:MF"),
             GO_CC = list(db_species = db_species_abbr, collection = "C5", subcollection = "GO:CC"),
             REACTOME = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:REACTOME"),
             WIKI = list(db_species = db_species_abbr, collection = "C2", subcollection = "CP:WIKIPATHWAYS")
         )
     }
     log_message(paste("Databases to analyze:", paste(names(databases), collapse=", ")))


     # --- Run GSEA for each contrast and database ---
    for (contrast in contrasts_to_analyze) {
        log_message(sprintf("\n--- Processing contrast: %s ---", contrast))
        
        # Extract DE results
        log_message(sprintf("  Extracting DE results for %s...", contrast))
        de_results <- topTable(fit, coef = contrast, 
                             sort.by = "t", adjust.method = "fdr", n = Inf)
        de_results <- de_results[rownames(de_results) != "", ]
        log_message(sprintf("  Found %d genes with DE results", nrow(de_results)))
        
        # Initialize results for this contrast
        gsea_results[[contrast]] <- list()
        
        # Run GSEA for each database
        for (db_name in names(databases)) {
            db <- databases[[db_name]]
            log_message(sprintf("  Running GSEA for %s...", db$name))
            
            tryCatch({
                result <- runGSEA(
                    de_results,
                    rank_metric = "t",
                    species = "Mus musculus",
                    category = db$category,
                    subcategory = db$subcategory,
                    padj_method = "fdr",
                    nperm = 100000,
                    pvalue_cutoff = pvalue_cutoff
                )
                
                # Check if any terms were enriched
                if (is.character(result)) {
                    log_message(sprintf("  No terms enriched for %s: %s", db$name, result), "WARNING")
                    gsea_results[[contrast]][[db_name]] <- NULL
                } else if (nrow(result@result) == 0) {
                    log_message(sprintf("  No terms enriched for %s", db$name), "WARNING")
                    gsea_results[[contrast]][[db_name]] <- NULL
                } else {
                    log_message(sprintf("  Found %d enriched terms for %s", nrow(result@result), db$name))
                    gsea_results[[contrast]][[db_name]] <- result
                }
            }, error = function(e) {
                log_message(sprintf("Error in %s GSEA for contrast %s: %s", 
                                   db$name, contrast, e$message), "ERROR")
                gsea_results[[contrast]][[db_name]] <- NULL
            })
        }
    }
    
    log_message("\n--- Identifying significant pathways across contrasts ---")
    
    # Initialize pools and genes lists
    pools <- list()
    genes <- list()
    
    # Collect significant pathways and genes for each database
    for (db_name in names(databases)) {
        tryCatch({
            # Get significant pathways
            db_results <- lapply(gsea_results, `[[`, db_name)
            pools[[db_name]] <- get_significant_pathways(db_results, q_cutoff = pvalue_cutoff)
            log_message(sprintf("Found %d significant pathways for %s", length(pools[[db_name]]), db_name))
            
            # Get pathway genes
            genes[[db_name]] <- get_pathway_genes_all(gsea_results, db_name, top = top_n)
            log_message(sprintf("Collected genes for %d pathways for %s", length(genes[[db_name]]), db_name))
        }, error = function(e) {
            log_message(sprintf("Error collecting pathways/genes for %s: %s", db_name, e$message), "ERROR")
            pools[[db_name]] <- character(0)
            genes[[db_name]] <- list()
        })
    }
    
    log_message("\n--- Calculating pathway scores ---")
    
    # Calculate normalized expression scores
    log_message("  Normalizing expression data (CPM)...")
    norm_counts <- cpm(DGEobject, log = TRUE)
    
    # Calculate scores for each database
    scores <- list()
    for (db_name in names(databases)) {
        tryCatch({
            if (length(genes[[db_name]]) > 0) {
                scores[[db_name]] <- calculate_pathway_scores(norm_counts, genes[[db_name]])
                log_message(sprintf("Calculated scores for %s: %d samples x %d pathways", 
                                   db_name, nrow(scores[[db_name]]), ncol(scores[[db_name]])))
            } else {
                log_message(sprintf("No pathway genes found for %s, skipping score calculation", db_name), "WARNING")
                scores[[db_name]] <- matrix(NA, 
                                          nrow = ncol(DGEobject), 
                                          ncol = 0,
                                          dimnames = list(colnames(DGEobject), character(0)))
            }
        }, error = function(e) {
            log_message(sprintf("Error calculating scores for %s: %s", db_name, e$message), "ERROR")
            scores[[db_name]] <- matrix(NA, 
                                      nrow = ncol(DGEobject), 
                                      ncol = 0,
                                      dimnames = list(colnames(DGEobject), character(0)))
        })
    }
    
    log_message("=== run_pooled_gsea() completed successfully ===")
    
    # Return results
    return(list(
        gsea_results = gsea_results,
        pools = pools,
        genes = genes,
        scores = scores
    ))
}
