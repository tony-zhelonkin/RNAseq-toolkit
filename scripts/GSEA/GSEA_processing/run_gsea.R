#' Run Gene Set Enrichment Analysis (GSEA) using msigdbr and clusterProfiler
#'
#' Performs GSEA on a pre-ranked gene list derived from differential expression results,
#' using gene sets fetched from the Molecular Signatures Database (MSigDB) via the
#' `msigdbr` package and the core analysis via `clusterProfiler::GSEA`.
#' Note: Implicitly requires a relevant organism annotation package (e.g., `org.Mm.eg.db`, `org.Hs.eg.db`) to be installed for `clusterProfiler` functionality
#' it might not be directly used in this specific function's code depending on gene ID types.
#'
#' @param DE_results A data frame containing differential expression results. Must have
#'        gene identifiers as rownames and a column specified by `rank_metric`.
#' @param rank_metric Character, the column name in `DE_results` to use for ranking genes
#'        (default: "t"). Genes are sorted decreasingly by this metric.
#' @param species Character, the species name used by `msigdbr` to fetch gene sets
#'        (e.g., "Mus musculus", "Homo sapiens"). Default: "Mus musculus".
#' @param db_species Character, the species abbreviation used by `msigdbr` (e.g., "MM", "HS").
#'        Default: "MM". Should generally match the `species` argument.
#' @param collection Character, the MSigDB collection code (e.g., "H", "MH", "C2", "C5").
#'        Default: "MH" (Hallmark for Mouse). Use `msigdbr::msigdbr_collections()` to see options.
#' @param subcollection Character, the MSigDB subcollection code (e.g., "CP:KEGG", "GO:BP").
#'        Default: "" (empty string, suitable for collections like Hallmark). Use "" or NULL if no subcollection.
#' @param pvalue_cutoff Numeric, p-value cutoff threshold used by `clusterProfiler::GSEA`
#'        (default: 1, meaning no filtering by nominal p-value at the GSEA step).
#' @param padj_method Character, p-value adjustment method used by `clusterProfiler::GSEA`
#'        (default: "fdr").
#' @param nperm Integer, number of permutations for the GSEA algorithm (`nPermSimple` in `fgsea`).
#'        Default: 100000.
#' @param seed Integer, random seed for reproducibility of permutations (default: 123).
#'
#' @return A `gseaResult` object from the `clusterProfiler` package.
#' @export
#' @import msigdbr clusterProfiler
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#' # Assuming DE_results is a data frame with t-statistics and gene symbols as rownames
#' # Run Hallmark GSEA for Mouse (default)
#' gsea_hallmark_mm <- run_gsea(DE_results)
#'
#' # Run KEGG GSEA for Human
#' gsea_kegg_hs <- run_gsea(DE_results_hs, species = "Homo sapiens", db_species = "HS",
#'                          collection = "C2", subcollection = "CP:KEGG")
#' }
# Note: Renaming function to snake_case for consistency
run_gsea <- function(
    DE_results,
    rank_metric   = "t",
    species       = "Mus musculus",
    db_species    = "MM",
    collection    = "MH",         # Default: Hallmark for Mouse
    subcollection = "",           # Default: empty for Hallmark
    pvalue_cutoff = 1,
    padj_method   = "fdr",
    nperm         = 100000,
    seed          = 123
) {
  # --- Dependencies ---
  # Explicitly load org.Mm.eg.db for now, consider making this conditional/parameterized later
  # if supporting multiple species robustly within this function.
  if (grepl("musculus", species, ignore.case = TRUE)) {
      library(org.Mm.eg.db)
  } else if (grepl("sapiens", species, ignore.case = TRUE)) {
      library(org.Hs.eg.db) # Assumes human annotation package is installed
  } else {
      warning("No specific organism DB loaded for species: ", species, ". GSEA might still work if gene IDs are symbols.")
  }
  library(msigdbr)
  library(clusterProfiler)
  # --------------------

  # --- Input Validation ---
  if (!is.data.frame(DE_results)) stop("`DE_results` must be a data frame.")
  if (is.null(rownames(DE_results))) stop("`DE_results` must have rownames (gene identifiers).")
  if (!rank_metric %in% colnames(DE_results)) stop(sprintf("`rank_metric` column '%s' not found in `DE_results`.", rank_metric))
  # --------------------

  # Prepare ranked gene list
  gene_vector <- DE_results[[rank_metric]]
  # Check for NAs in ranking metric
  if (any(is.na(gene_vector))) {
      warning("NA values found in rank_metric column '", rank_metric, "'. Removing corresponding genes.")
      valid_indices <- !is.na(gene_vector)
      gene_vector <- gene_vector[valid_indices]
      gene_names <- rownames(DE_results)[valid_indices]
  } else {
      gene_names <- rownames(DE_results)
  }
  names(gene_vector) <- gene_names
  ranked_genes <- sort(gene_vector, decreasing = TRUE)

  # Remove infinite values if any
   if (any(is.infinite(ranked_genes))) {
       warning("Infinite values found in ranked gene list. Removing them.")
       ranked_genes <- ranked_genes[!is.infinite(ranked_genes)]
   }

  # Retrieve gene sets using msigdbr
  message(sprintf("Fetching MSigDB sets for species='%s', db_species='%s', collection='%s', subcollection='%s'",
                  species, db_species, collection, subcollection))
  msigdb_df <- msigdbr::msigdbr(
      species       = species,
      # db_species is not a direct argument, species handles this
      category      = collection,    # Use 'category' for msigdbr main collection
      subcategory   = subcollection  # Use 'subcategory' for msigdbr subcollection
  )

  if(nrow(msigdb_df) == 0) {
      stop(sprintf("No gene sets found via msigdbr for species='%s', category='%s', subcategory='%s'.\nCheck available collections with msigdbr::msigdbr_collections()",
                   species, collection, subcollection))
  }

  # Prepare TERM2GENE dataframe
  term2gene_df <- msigdb_df[, c("gs_name", "gene_symbol")] # Assuming input uses gene symbols

  message("Running clusterProfiler::GSEA...")
  stats::set.seed(seed)
  GSEA_result <- clusterProfiler::GSEA(
      geneList = ranked_genes,
      TERM2GENE = term2gene_df,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = padj_method,
      eps = 0, # Recommended setting for fgsea
      by = "fgsea", # Use fgsea implementation
      nPermSimple = nperm, # Pass permutation number to fgsea
      verbose = FALSE # Keep GSEA quiet
  )

  # Check if GSEA result is valid
   if (is.null(GSEA_result) || nrow(GSEA_result@result) == 0) {
       warning(sprintf("GSEA returned no significant results for collection='%s', subcollection='%s'.",
                       collection, subcollection))
       # Return the empty/NULL object as is
   } else {
        message("GSEA completed successfully.")
   }

  return(GSEA_result)
}
