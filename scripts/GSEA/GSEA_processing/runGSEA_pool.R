#' Run Pooled GSEA Analysis Across Multiple Contrasts
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) across multiple contrasts,
#' aggregates significant pathways, and calculates pathway scores.
#'
#' @param fit A limma fit object containing differential expression results.
#' @param contrasts A contrast matrix created with makeContrasts().
#' @param DGErankl A DGEList object containing normalized expression data.
#' @param top_n Integer, number of top pathways to include in the results (default: 30).
#' @param verbose Logical, whether to print progress messages (default: TRUE).
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
#' # and DGErankl is a DGEList object
#' pooled_results <- run_pooled_gsea(fit, contrasts, DGErankl, top_n = 20)
source("scripts/GSEA/GSEA_processing/get_pathway_genes.R")
source("scripts/GSEA/GSEA_processing/get_significant_pathways.R")
source("scripts/GSEA/GSEA_processing/get_pathway_genes_all.R")
source("scripts/GSEA/GSEA_processing/calculate_pathway_scores.R")

run_pooled_gsea <- function(fit, contrasts, DGErankl, top_n = 30, verbose = TRUE) {
    # Run GSEA for all pool contrasts
    if(verbose) message("Starting pooled GSEA analysis...")
    gsea_results <- list()
    contrasts_to_analyze <- colnames(contrasts)
    
    if(verbose) message(paste("Found", length(contrasts_to_analyze), "contrasts to analyze:", paste(contrasts_to_analyze, collapse=", ")))
    
    # Run GSEA analysis for each contrast
    for(contrast in contrasts_to_analyze) {
        if(verbose) message(paste("\n--- Processing contrast:", contrast, "---"))
        
        if(verbose) message("  Extracting DE results...")
        de_results <- topTable(fit, coef = contrast, 
                             sort.by = "t", adjust.method = "fdr", n = Inf)
        de_results <- de_results[rownames(de_results) != "", ]
        if(verbose) message(paste("  Found", nrow(de_results), "genes with DE results"))
        
        tryCatch({
            if(verbose) message("  Running GSEA for Hallmark gene sets...")
            gsea_results[[contrast]][["hallmark"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "H", padj_method = "fdr", nperm = 100000,
                pvalue_cutoff = 0.01
            )
            if(verbose) message(paste("    Completed Hallmark GSEA with", nrow(gsea_results[[contrast]][["hallmark"]]), "pathways"))
            
            if(verbose) message("  Running GSEA for KEGG pathways...")
            gsea_results[[contrast]][["kegg"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C2", subcategory = "CP:KEGG", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
            if(verbose) message(paste("    Completed KEGG GSEA with", nrow(gsea_results[[contrast]][["kegg"]]), "pathways"))
            
            if(verbose) message("  Running GSEA for GO Biological Process terms...")
            gsea_results[[contrast]][["gobp"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C5", subcategory = "GO:BP", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
            if(verbose) message(paste("    Completed GO BP GSEA with", nrow(gsea_results[[contrast]][["gobp"]]), "pathways"))
            
            if(verbose) message("  Running GSEA for Reactome pathways...")
            gsea_results[[contrast]][["reactome"]] <- runGSEA(
                de_results, rank_metric = "t", species = "Mus musculus",
                category = "C2", subcategory = "CP:REACTOME", padj_method = "fdr",
                nperm = 100000, pvalue_cutoff = 0.01
            )
            if(verbose) message(paste("    Completed Reactome GSEA with", nrow(gsea_results[[contrast]][["reactome"]]), "pathways"))
            
        }, error = function(e) {
            message(sprintf("Error in contrast %s: %s", contrast, e$message))
        })
    }
    
    if(verbose) message("\n--- Identifying significant pathways across contrasts ---")
    # Get pools for each database
    pools <- list(
        hallmark = get_significant_pathways(lapply(gsea_results, `[[`, "hallmark")),
        kegg = get_significant_pathways(lapply(gsea_results, `[[`, "kegg")),
        gobp = get_significant_pathways(lapply(gsea_results, `[[`, "gobp")),
        reactome = get_significant_pathways(lapply(gsea_results, `[[`, "reactome"))
    )
    
    if(verbose) {
        message(paste("  Identified", length(pools$hallmark), "significant Hallmark pathways"))
        message(paste("  Identified", length(pools$kegg), "significant KEGG pathways"))
        message(paste("  Identified", length(pools$gobp), "significant GO BP terms"))
        message(paste("  Identified", length(pools$reactome), "significant Reactome pathways"))
    }
    
    if(verbose) message("\n--- Extracting pathway genes ---")
    # Get pathway genes
    genes <- list(
        hallmark = get_pathway_genes_all(gsea_results, "hallmark", top = top_n),
        kegg = get_pathway_genes_all(gsea_results, "kegg", top = top_n),
        gobp = get_pathway_genes_all(gsea_results, "gobp", top = top_n),
        reactome = get_pathway_genes_all(gsea_results, "reactome", top = top_n)
    )
    
    if(verbose) {
        message(paste("  Extracted genes for", length(genes$hallmark), "Hallmark pathways"))
        message(paste("  Extracted genes for", length(genes$kegg), "KEGG pathways"))
        message(paste("  Extracted genes for", length(genes$gobp), "GO BP terms"))
        message(paste("  Extracted genes for", length(genes$reactome), "Reactome pathways"))
    }
    
    if(verbose) message("\n--- Calculating pathway scores ---")
    # Calculate normalized expression scores
    if(verbose) message("  Normalizing expression data (CPM)...")
    norm_counts <- cpm(DGErankl, log = TRUE)
    
    scores <- list(
        hallmark = calculate_pathway_scores(norm_counts, genes$hallmark),
        kegg = calculate_pathway_scores(norm_counts, genes$kegg),
        gobp = calculate_pathway_scores(norm_counts, genes$gobp),
        reactome = calculate_pathway_scores(norm_counts, genes$reactome)
    )
    
    if(verbose) message("  Pathway scores calculated successfully")
    
    if(verbose) message("\nGSEA analysis completed successfully!")
    
    return(list(
        gsea_results = gsea_results,
        pools = pools,
        genes = genes,
        scores = scores
    ))
}
