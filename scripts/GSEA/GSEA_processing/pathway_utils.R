#' Pathway Utility Functions for GSEA Analysis
#'
#' Collection of helper functions for working with gene sets and pathways.
#' Designed for use with clusterProfiler-style TERM2GENE/TERM2NAME data frames.
#'
#' @name pathway_utils
#' @keywords internal
NULL


#' Filter Gene Sets by Size
#'
#' Filters TERM2GENE and TERM2NAME data frames to retain only gene sets
#' within specified size bounds. Essential for GSEA analysis where very
#' small or very large gene sets can produce unreliable results.
#'
#' @param T2G Data frame with columns 'gs_name' (or 'term') and 'gene_symbol' (or 'gene')
#' @param T2N Data frame with columns 'gs_name' (or 'term') and 'description' (or 'name')
#' @param min_size Integer, minimum number of genes per set (default: 5)
#' @param max_size Integer, maximum number of genes per set (default: 500)
#' @param verbose Logical, print filtering statistics (default: TRUE)
#'
#' @return List with filtered T2G and T2N data frames
#'
#' @examples
#' filtered <- filter_pathways_by_size(
#'   T2G = transport_db$T2G,
#'   T2N = transport_db$T2N,
#'   min_size = 10,
#'   max_size = 200
#' )
#'
#' @export
filter_pathways_by_size <- function(
    T2G,
    T2N = NULL,
    min_size = 5,
    max_size = 500,
    verbose = TRUE
) {

  # Detect column names (support both conventions)
  gs_col <- if ("gs_name" %in% colnames(T2G)) "gs_name" else "term"
  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"

  if (!gs_col %in% colnames(T2G) || !gene_col %in% colnames(T2G)) {
    stop("T2G must have columns for gene set names and gene symbols")
  }

  # Count genes per pathway
  pathway_sizes <- table(T2G[[gs_col]])

  # Identify pathways within size bounds
  valid_pathways <- names(pathway_sizes)[
    pathway_sizes >= min_size & pathway_sizes <= max_size
  ]

  n_original <- length(unique(T2G[[gs_col]]))
  n_retained <- length(valid_pathways)

  if (verbose) {
    message(sprintf("Filtering gene sets by size [%d, %d]:", min_size, max_size))
    message(sprintf("  Original: %d gene sets", n_original))
    message(sprintf("  Too small (<%d): %d", min_size, sum(pathway_sizes < min_size)))
    message(sprintf("  Too large (>%d): %d", max_size, sum(pathway_sizes > max_size)))
    message(sprintf("  Retained: %d gene sets (%.1f%%)",
                    n_retained, 100 * n_retained / n_original))
  }

  # Filter T2G
  T2G_filtered <- T2G[T2G[[gs_col]] %in% valid_pathways, ]

  # Filter T2N if provided
  T2N_filtered <- NULL
  if (!is.null(T2N)) {
    gs_col_t2n <- if ("gs_name" %in% colnames(T2N)) "gs_name" else "term"
    T2N_filtered <- T2N[T2N[[gs_col_t2n]] %in% valid_pathways, ]
  }

  return(list(
    T2G = T2G_filtered,
    T2N = T2N_filtered,
    n_original = n_original,
    n_retained = n_retained
  ))
}


#' Get Pathway Size Distribution
#'
#' Returns summary statistics about gene set sizes.
#'
#' @param T2G Data frame with gene set-gene mappings
#'
#' @return Named vector with size statistics
#' @export
get_pathway_size_stats <- function(T2G) {
  gs_col <- if ("gs_name" %in% colnames(T2G)) "gs_name" else "term"

  pathway_sizes <- table(T2G[[gs_col]])

  stats <- c(
    n_pathways = length(pathway_sizes),
    min_size = min(pathway_sizes),
    q1_size = as.numeric(quantile(pathway_sizes, 0.25)),
    median_size = as.numeric(median(pathway_sizes)),
    mean_size = round(mean(pathway_sizes), 1),
    q3_size = as.numeric(quantile(pathway_sizes, 0.75)),
    max_size = max(pathway_sizes),
    n_genes = length(unique(T2G[[if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"]]))
  )

  return(stats)
}


#' Convert Gene Set List to TERM2GENE Format
#'
#' Converts a named list of gene sets to the TERM2GENE data frame format
#' expected by clusterProfiler::GSEA().
#'
#' @param geneset_list Named list where names are pathway IDs and values are character vectors of genes
#' @param gs_col Character, name for the gene set column (default: "gs_name")
#' @param gene_col Character, name for the gene column (default: "gene_symbol")
#'
#' @return Data frame in TERM2GENE format
#'
#' @examples
#' genesets <- list(
#'   "PATHWAY_A" = c("Gene1", "Gene2", "Gene3"),
#'   "PATHWAY_B" = c("Gene4", "Gene5")
#' )
#' T2G <- list_to_term2gene(genesets)
#'
#' @export
list_to_term2gene <- function(
    geneset_list,
    gs_col = "gs_name",
    gene_col = "gene_symbol"
) {
  if (!is.list(geneset_list) || is.null(names(geneset_list))) {
    stop("geneset_list must be a named list")
  }

  # Build data frame
  df <- do.call(rbind, lapply(names(geneset_list), function(pathway) {
    genes <- geneset_list[[pathway]]
    data.frame(
      gs_name = rep(pathway, length(genes)),
      gene_symbol = genes,
      stringsAsFactors = FALSE
    )
  }))

  # Rename columns if requested
  colnames(df) <- c(gs_col, gene_col)

  return(df)
}


#' Convert TERM2GENE to Gene Set List
#'
#' Converts a TERM2GENE data frame to a named list format.
#'
#' @param T2G Data frame with gene set-gene mappings
#'
#' @return Named list of gene sets
#' @export
term2gene_to_list <- function(T2G) {
  gs_col <- if ("gs_name" %in% colnames(T2G)) "gs_name" else "term"
  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"

  split(T2G[[gene_col]], T2G[[gs_col]])
}


#' Create TERM2NAME from TERM2GENE
#'
#' Creates a TERM2NAME data frame from TERM2GENE when descriptions are not available.
#' Uses the gene set name as the description, with optional formatting.
#'
#' @param T2G Data frame with gene set-gene mappings
#' @param format_names Logical, apply name formatting (default: TRUE)
#'
#' @return Data frame with gs_name and description columns
#' @export
create_term2name <- function(T2G, format_names = TRUE) {
  gs_col <- if ("gs_name" %in% colnames(T2G)) "gs_name" else "term"

  unique_sets <- unique(T2G[[gs_col]])

  if (format_names) {
    # Basic formatting: replace underscores, title case
    descriptions <- sapply(unique_sets, function(x) {
      x <- gsub("_", " ", x)
      # Remove common prefixes
      x <- gsub("^(HALLMARK|KEGG|REACTOME|GOBP|GOCC|GOMF|WP|TRANSPORTDB)\\s+", "", x, ignore.case = TRUE)
      # Title case
      paste(sapply(strsplit(tolower(x), " ")[[1]], function(w) {
        if (nchar(w) > 0) paste0(toupper(substr(w, 1, 1)), substr(w, 2, nchar(w))) else w
      }), collapse = " ")
    }, USE.NAMES = FALSE)
  } else {
    descriptions <- unique_sets
  }

  data.frame(
    gs_name = unique_sets,
    description = descriptions,
    stringsAsFactors = FALSE
  )
}


#' Export Gene Sets to GMX Format
#'
#' Exports gene sets to GMX format compatible with GSEA desktop application
#' and other external tools.
#'
#' @param T2G Data frame with gene set-gene mappings
#' @param T2N Data frame with gene set descriptions (optional)
#' @param output_file Character, path to output file
#'
#' @return Invisibly returns the output file path
#' @export
export_to_gmx <- function(T2G, T2N = NULL, output_file) {
  gs_col <- if ("gs_name" %in% colnames(T2G)) "gs_name" else "term"
  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"

  # Get gene lists per pathway
  gene_lists <- split(T2G[[gene_col]], T2G[[gs_col]])
  pathways <- names(gene_lists)

  # Get descriptions
  if (!is.null(T2N)) {
    desc_col <- if ("description" %in% colnames(T2N)) "description" else "name"
    gs_col_t2n <- if ("gs_name" %in% colnames(T2N)) "gs_name" else "term"
    descs <- T2N[[desc_col]][match(pathways, T2N[[gs_col_t2n]])]
    descs[is.na(descs)] <- pathways[is.na(descs)]
  } else {
    descs <- pathways
  }

  # Build GMX matrix
  max_len <- max(sapply(gene_lists, length))
  gmx_mat <- matrix("", nrow = max_len + 2, ncol = length(pathways))

  # Row 1: Descriptions
  gmx_mat[1, ] <- descs
  # Row 2: Pathway IDs
  gmx_mat[2, ] <- pathways
  # Rows 3+: Genes
  for (i in seq_along(gene_lists)) {
    genes <- gene_lists[[i]]
    if (length(genes) > 0) {
      gmx_mat[3:(2 + length(genes)), i] <- genes
    }
  }

  # Write to file
  write.table(
    gmx_mat,
    output_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  message(sprintf("Exported %d gene sets to: %s", length(pathways), output_file))
  invisible(output_file)
}


message("[RNAseq-toolkit] Loaded: pathway_utils.R")
