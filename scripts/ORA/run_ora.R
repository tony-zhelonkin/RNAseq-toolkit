# =============================================================================
# ORA: Over-Representation Analysis Functions
# =============================================================================
#
# PURPOSE:
#   Wrapper functions for GO/KEGG over-representation analysis using
#   clusterProfiler. Designed to complement GSEA analysis.
#
# KEY DISTINCTION FROM GSEA:
#   - ORA: Binary (gene is in list or not). Best for module validation.
#   - GSEA: Continuous ranking. Best for detecting subtle pathway shifts.
#
# USAGE:
#   source("scripts/ORA/run_ora.R")
#   result <- run_ora(gene_list, species = "Mus musculus", ont = "BP")
#
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(tibble)
})

#' Run GO Over-Representation Analysis
#'
#' @param gene_list Character vector of gene symbols
#' @param species Species name ("Mus musculus" or "Homo sapiens")
#' @param ont GO ontology: "BP", "MF", "CC", or "ALL"
#' @param pvalue_cutoff P-value cutoff (default 0.05)
#' @param qvalue_cutoff Q-value cutoff (default 0.2)
#' @param min_gs_size Minimum gene set size (default 10)
#' @param max_gs_size Maximum gene set size (default 500)
#' @return enrichResult object or NULL if no enrichment found
#'
#' @examples
#' \dontrun{
#' genes <- c("Irf8", "Batf3", "Clec9a", "Xcr1", "Cd8a")
#' result <- run_ora(genes, species = "Mus musculus", ont = "BP")
#' }
run_ora <- function(gene_list,
                    species = "Mus musculus",
                    ont = "BP",
                    pvalue_cutoff = 0.05,
                    qvalue_cutoff = 0.2,
                    min_gs_size = 10,
                    max_gs_size = 500) {

  # Validate inputs
  if (length(gene_list) < 5) {
    warning("Gene list has fewer than 5 genes - ORA may not be meaningful")
  }

  # Load appropriate organism database
  if (grepl("musculus|mouse", species, ignore.case = TRUE)) {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("org.Mm.eg.db required for mouse. Install via BiocManager::install('org.Mm.eg.db')")
    }
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (grepl("sapiens|human", species, ignore.case = TRUE)) {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("org.Hs.eg.db required for human. Install via BiocManager::install('org.Hs.eg.db')")
    }
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    stop("Species must be 'Mus musculus' or 'Homo sapiens'")
  }

  # Remove NA values
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- unique(gene_list)

  # Run enrichGO
  tryCatch({
    result <- enrichGO(
      gene = gene_list,
      OrgDb = org_db,
      keyType = "SYMBOL",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = TRUE
    )

    if (is.null(result) || nrow(result@result) == 0) {
      message("No significant GO enrichment found")
      return(NULL)
    }

    return(result)

  }, error = function(e) {
    warning("ORA failed: ", conditionMessage(e))
    return(NULL)
  })
}


#' Run KEGG Over-Representation Analysis
#'
#' @param gene_list Character vector of gene symbols
#' @param species Species name ("Mus musculus" or "Homo sapiens")
#' @param pvalue_cutoff P-value cutoff (default 0.05)
#' @param qvalue_cutoff Q-value cutoff (default 0.2)
#' @return enrichResult object or NULL if no enrichment found
run_ora_kegg <- function(gene_list,
                         species = "Mus musculus",
                         pvalue_cutoff = 0.05,
                         qvalue_cutoff = 0.2) {

  # Determine organism code and database
  if (grepl("musculus|mouse", species, ignore.case = TRUE)) {
    org_code <- "mmu"
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("org.Mm.eg.db required")
    }
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (grepl("sapiens|human", species, ignore.case = TRUE)) {
    org_code <- "hsa"
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("org.Hs.eg.db required")
    }
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    stop("Species must be 'Mus musculus' or 'Homo sapiens'")
  }

  # Convert symbols to Entrez IDs (KEGG requires Entrez)
  gene_list <- unique(gene_list[!is.na(gene_list)])

  entrez_ids <- tryCatch({
    AnnotationDbi::mapIds(
      org_db,
      keys = gene_list,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
  }, error = function(e) {
    warning("Failed to map gene symbols to Entrez IDs")
    return(NULL)
  })

  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  if (length(entrez_ids) < 5) {
    warning("Too few genes mapped to Entrez IDs")
    return(NULL)
  }

  # Run enrichKEGG
  tryCatch({
    result <- enrichKEGG(
      gene = entrez_ids,
      organism = org_code,
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff
    )

    if (is.null(result) || nrow(result@result) == 0) {
      message("No significant KEGG enrichment found")
      return(NULL)
    }

    return(result)

  }, error = function(e) {
    warning("KEGG ORA failed: ", conditionMessage(e))
    return(NULL)
  })
}


#' Normalize ORA results to standard tibble format
#'
#' Converts enrichResult object to a standardized tibble for easy
#' comparison and aggregation across multiple analyses.
#'
#' @param ora_result enrichResult object from run_ora() or run_ora_kegg()
#' @param database Database name for labeling (e.g., "GO_BP", "KEGG")
#' @param module Module/gene list name for labeling
#' @param padj_cutoff Filter to significant results (default 1 = no filter)
#' @param format_names Clean up pathway names (default TRUE)
#' @param max_name_length Maximum name length before truncation (default 60)
#' @return tibble with standardized columns
normalize_ora_results <- function(ora_result,
                                   database = "GO",
                                   module = "unknown",
                                   padj_cutoff = 1,
                                   format_names = TRUE,
                                   max_name_length = 60) {

  if (is.null(ora_result) || nrow(ora_result@result) == 0) {
    return(tibble())
  }

  df <- ora_result@result %>%
    as_tibble() %>%
    dplyr::filter(p.adjust <= padj_cutoff) %>%
    dplyr::mutate(
      database = database,
      module = module,
      # Parse GeneRatio (e.g., "5/100" -> 0.05)
      gene_ratio_numeric = sapply(GeneRatio, function(x) {
        parts <- as.numeric(strsplit(x, "/")[[1]])
        if (length(parts) == 2) parts[1]/parts[2] else NA
      }),
      neg_log_padj = -log10(p.adjust),
      # Clean names if requested
      pathway_name = if (format_names) {
        substr(gsub("_", " ", Description), 1, max_name_length)
      } else {
        Description
      }
    ) %>%
    dplyr::select(
      pathway_id = ID,
      pathway_name,
      database,
      module,
      gene_ratio = GeneRatio,
      gene_ratio_numeric,
      pvalue,
      padj = p.adjust,
      qvalue,
      neg_log_padj,
      gene_count = Count,
      genes = geneID
    ) %>%
    dplyr::arrange(padj)

  return(df)
}


#' Run ORA across multiple ontologies
#'
#' Convenience function to run GO ORA across BP, MF, and CC ontologies
#' and combine results.
#'
#' @param gene_list Character vector of gene symbols
#' @param species Species name
#' @param ... Additional arguments passed to run_ora()
#' @return Named list of enrichResult objects
run_ora_all_ontologies <- function(gene_list, species = "Mus musculus", ...) {
  ontologies <- c("BP", "MF", "CC")
  results <- list()

  for (ont in ontologies) {
    message("Running GO:", ont)
    results[[paste0("GO_", ont)]] <- run_ora(gene_list, species, ont = ont, ...)
  }

  return(results)
}


# =============================================================================
# MESSAGE ON LOAD
# =============================================================================

message("[RNAseq-toolkit] Loaded: run_ora.R (run_ora, run_ora_kegg, normalize_ora_results)")
