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