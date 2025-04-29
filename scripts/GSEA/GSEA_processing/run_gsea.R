#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on differential expression results
#' using gene sets from the Molecular Signatures Database (MSigDB).
#'
#' @param DE_results A data frame containing differential expression results with gene identifiers
#'        as rownames and a column for the ranking metric (e.g., t-statistic).
#' @param rank_metric Character, name of the column to use for ranking genes (default: "t").
#' @param species Character, species name for MSigDB gene sets (default: "Mus musculus").
#' @param category Character, MSigDB collection category (default: "H" for Hallmark gene sets).
#'        Options include: "H" (Hallmark), "C1" (positional), "C2" (curated), "C3" (motif),
#'        "C4" (computational), "C5" (GO), "C6" (oncogenic), "C7" (immunologic), "C8" (cell type).
#' @param subcategory Character, MSigDB subcategory (default: NULL).
#'        For example, "CP:KEGG", "CP:REACTOME", "GO:BP", "GO:MF", "GO:CC".
#' @param pvalue_cutoff Numeric, p-value cutoff for significance (default: 1).
#' @param padj_method Character, method for p-value adjustment (default: "fdr").
#' @param nperm Integer, number of permutations for GSEA (default: 100000).
#' @param seed Integer, random seed for reproducibility (default: 123).
#'
#' @return A GSEA result object from the clusterProfiler package.
#' @export
#'
#' @examples
#' # Assuming DE_results is a data frame with differential expression results
#' gsea_hallmark <- runGSEA(DE_results, rank_metric = "t", category = "H")
#' gsea_kegg <- runGSEA(DE_results, rank_metric = "t", category = "C2", subcategory = "CP:KEGG")
runGSEA <- function(DE_results,
                    rank_metric = "t",
                    species = "Mus musculus",
                    category = "H",
                    subcategory = NULL,
                    pvalue_cutoff = 1,
                    padj_method = "fdr",
                    nperm = 100000,
                    seed = 123) {
  
  library(org.Mm.eg.db)
  library(msigdbr)
  library(clusterProfiler)
  
  # Prepare the ranked gene list from DE results
  ranked_genes <- DE_results[[rank_metric]]
  names(ranked_genes) <- rownames(DE_results)
  
  # Sort genes
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)
  
  # Get gene sets from MSigDB
  msigdb_H <- msigdbr(species = species, 
                      category = category,
                      subcategory = subcategory)
  
  # Perform GSEA
  GSEA_result <- GSEA(ranked_genes, 
                TERM2GENE = msigdb_H[, c("gs_name", "gene_symbol")], 
                pvalueCutoff = pvalue_cutoff, 
                verbose = FALSE,
                pAdjustMethod = padj_method,
                eps = 0,
                by = "fgsea",
                seed = seed,
                nPermSimple = nperm)
  
  return(GSEA_result)
}
