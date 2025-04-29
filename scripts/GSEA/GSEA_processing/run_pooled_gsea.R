#' Run Pooled GSEA Analysis Across Multiple Contrasts
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) across multiple contrasts,
#' aggregates significant pathways, and calculates pathway scores.
#'
#' @param fit A limma fit object containing differential expression results.
#' @param contrasts A contrast matrix created with makeContrasts().
#' @param DGEobject A DGEList object containing normalized expression data.
#' @param top_n Integer, number of top pathways to include in the results (default: 30).
#' @param pvalue_cutoff Numeric, p-value cutoff for significance (default: 0.01).
#' @param verbose Logical, whether to print progress messages (default: TRUE).
#' @param log_file Character, path to log file (default: NULL, logs to console only).
#'
#' @return A list containing:
#'   \item{gsea_results}{A nested list of GSEA results by contrast and database.}
#'   \item{pools}{A list of significant pathway IDs by database.}
#'   \item{genes}{A list of pathway genes by database.}
#'   \item{scores}{A list of pathway score matrices by database.}
#' @export
#'
#' @examples
#' # Assuming fit is a limma fit object, contrasts is a contrast matrix,
#' # and DGEobject is a DGEList object
#' pooled_results <- run_pooled_gsea(fit, contrasts, DGEobject, top_n = 20)
source("scripts/GSEA/GSEA_processing/run_gsea.R")
source("scripts/GSEA/GSEA_processing/get_pathway_genes.R")
source("scripts/GSEA/GSEA_processing/get_significant_pathways.R")
source("scripts/GSEA/GSEA_processing/get_pathway_genes_all.R")
source("scripts/GSEA/GSEA_processing/calculate_pathway_scores.R")

run_pooled_gsea <- function(fit, 
                           contrasts, 
                           DGEobject, 
                           top_n = 30,
                           pvalue_cutoff = 0.01,
                           verbose = TRUE,
                           log_file = NULL) {
    
    # Set up logging
    log_conn <- NULL
    if (!is.null(log_file)) {
        log_conn <- file(log_file, open = "w")
        on.exit(if (!is.null(log_conn)) close(log_conn))  # Ensure file is closed when function exits
    }
    
    # Function to log messages to both console and file
    log_message <- function(msg, type = "INFO") {
        if (!verbose && type != "ERROR") return()
        
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        full_msg <- sprintf("[%s] [%s] %s", timestamp, type, msg)
        
        message(full_msg)
        if (!is.null(log_conn)) {
            writeLines(full_msg, log_conn)
        }
    }
    
    # Start logging
    log_message("=== Starting run_pooled_gsea ===")
    log_message(sprintf("Top N pathways: %d, p-value cutoff: %.3g", top_n, pvalue_cutoff))
    
    # Debug information about input data
    if (verbose) {
        log_message(sprintf("DGEobject contains %d genes and %d samples", 
                           nrow(DGEobject), ncol(DGEobject)), "DEBUG")
        log_message(sprintf("Sample names in DGEobject: %s", 
                           paste(head(colnames(DGEobject), 5), collapse=", ")), "DEBUG")
    }
    
    # Initialize results
    gsea_results <- list()
    contrasts_to_analyze <- colnames(contrasts)
    log_message(paste("Contrasts to analyze:", paste(contrasts_to_analyze, collapse=", ")))
    
    # Define databases to analyze
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
    
    # Run GSEA for each contrast
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
