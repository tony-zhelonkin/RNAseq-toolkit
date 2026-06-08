#' Normalize GSEA Results to Standard DataFrame
#'
#' Extracts and standardizes results from a gseaResult object into a consistent
#' DataFrame format suitable for aggregation and visualization.
#'
#' This is a key component of the "Normalize-Then-Visualize" pattern:
#'   1. Each GSEA analysis (MSigDB, Mitochondria, etc.) produces gseaResult objects
#'   2. This function normalizes them to a common DataFrame schema
#'   3. DataFrames can be rbind'd into a master table
#'   4. All visualization code reads from the master table (single source of truth)
#'
#' @param gsea_obj A gseaResult object from clusterProfiler/fgsea, or a data.frame
#' @param database Character, name of the source database (e.g., "Hallmark", "KEGG", "Mitochondria")
#' @param contrast Character, name of the contrast (e.g., "Treatment_vs_Control")
#' @param padj_cutoff Numeric, filter to significant pathways (default: 1 = no filter)
#' @param format_names Logical, apply smart pathway name formatting (default: TRUE)
#' @param max_name_length Integer, truncate pathway names to this length (default: 80)
#'
#' @return A tibble with standardized columns:
#'   - pathway_id: Original pathway ID
#'   - pathway_name: Cleaned/formatted pathway name
#'   - database: Source database name
#'   - contrast: Contrast name
#'   - NES: Normalized enrichment score
#'   - pvalue: Raw p-value
#'   - padj: Adjusted p-value (FDR)
#'   - set_size: Number of genes in the pathway
#'   - leading_edge_size: Number of genes in leading edge
#'   - gene_ratio: leading_edge_size / set_size
#'   - core_enrichment: Leading edge genes (slash-separated)
#'   - direction: "Up" or "Down" based on NES sign
#'   - neg_log_padj: -log10(padj) for visualization
#'
#' @examples
#' # Normalize Hallmark GSEA results
#' hallmark_df <- normalize_gsea_results(
#'   gsea_obj = gsea_H,
#'   database = "Hallmark",
#'   contrast = "Treatment_vs_Control"
#' )
#'
#' # Combine multiple databases
#' master_df <- bind_rows(
#'   normalize_gsea_results(gsea_H, "Hallmark", "Treatment_vs_Control"),
#'   normalize_gsea_results(gsea_kegg, "KEGG", "Treatment_vs_Control"),
#'   normalize_gsea_results(gsea_mito, "Mitochondria", "Treatment_vs_Control")
#' )
#'
#' @export
normalize_gsea_results <- function(
    gsea_obj,
    database,
    contrast,
    padj_cutoff = 1,
    format_names = TRUE,
    max_name_length = 80,
    atlas_universe = NULL
) {

  # Validate input
  if (is.null(gsea_obj)) {
    warning(sprintf("NULL gsea_obj provided for database '%s', returning empty tibble", database))
    return(empty_gsea_tibble())
  }

  # Extract result slot
  if (methods::is(gsea_obj, "gseaResult")) {
    result_df <- as.data.frame(gsea_obj@result)
  } else if (is.data.frame(gsea_obj)) {
    result_df <- gsea_obj
  } else {
    warning(sprintf("Unexpected object type for database '%s': %s", database, class(gsea_obj)[1]))
    return(empty_gsea_tibble())
  }

  # Check for empty results
  if (nrow(result_df) == 0) {
    message(sprintf("No results in gsea_obj for database '%s'", database))
    return(empty_gsea_tibble())
  }

  # Determine p-value column (different GSEA implementations use different names)
  padj_col <- if ("p.adjust" %in% colnames(result_df)) {
    "p.adjust"
  } else if ("qvalue" %in% colnames(result_df)) {
    "qvalue"
  } else if ("padj" %in% colnames(result_df)) {
    "padj"
  } else {
    stop("Could not find adjusted p-value column (p.adjust, qvalue, or padj)")
  }

  # Filter by padj_cutoff if specified
  if (padj_cutoff < 1) {
    result_df <- result_df[result_df[[padj_col]] < padj_cutoff, ]
    if (nrow(result_df) == 0) {
      message(sprintf("No significant results (padj < %.2f) for database '%s'", padj_cutoff, database))
      return(empty_gsea_tibble())
    }
  }

  # Calculate leading edge size from core_enrichment
  result_df$leading_edge_size <- sapply(result_df$core_enrichment, function(x) {
    if (is.na(x) || x == "") return(0L)
    length(strsplit(as.character(x), "/")[[1]])
  })

  # --- genes_full_set semantics for GSEA rows (MADR-008, atlas-universe fix) ---
  # For unweighted gene-set methods (MSigDB, custom GMT, CoReSh-derived GMT,
  # MitoPathways), `genes_full_set` is the full pathway membership intersected
  # with a *gene universe*. Two universes are possible; which is used depends
  # on the `atlas_universe` argument:
  #
  #   (a) atlas_universe = NULL (legacy default):
  #         universe = names(gsea_obj@geneList) — i.e. the per-contrast ranked
  #         list, which equals the post-filterByExpr DE universe for THIS
  #         contrast. Under per-contrast filterByExpr (1.1_pseudobulk_de.R),
  #         this universe varies with sequencing depth and cluster size
  #         (empirically 26x range across 53 contrasts in cdc1_path), making
  #         genes_full_set silently contrast-dependent and breaking
  #         cross-contrast stability of pathway-explorer's similarity layer.
  #
  #   (b) atlas_universe = <character vector of HGNC symbols> (preferred):
  #         universe = atlas_universe, the union of every contrast's
  #         post-filterByExpr universe (written by 1.1_pseudobulk_de.R as
  #         tables/atlas_gene_universe.txt). genes_full_set then depends only
  #         on pathway membership and the atlas universe — both contrast-
  #         invariant. This is what pathway-explorer's geometry layer needs.
  #
  # `core_enrichment` (above) remains the *leading edge* — the contrast-
  # dependent gene subset that drove the running sum past its extremum. By
  # construction `core_enrichment ⊆ genes_full_set` only when the leading
  # edge sits inside the chosen universe; under (b) some pathway members may
  # have been filtered out of the per-contrast ranked list and would not
  # contribute to core_enrichment, but they DO appear in genes_full_set.
  # Both columns are emitted; they are NOT the same value here (unlike for
  # PROGENy/TF rows where MLM/ULM produce no leading-edge concept and the
  # columns coincide).
  #
  # Source of truth: clusterProfiler stores the materialised TERM2GENE list at
  # gsea_obj@geneSets (named list keyed by pathway ID). Per-contrast universe
  # comes from gsea_obj@geneList; atlas universe comes from the caller.
  # Originating concern: pathway-geometry/synthesis.md P0-1.
  if (methods::is(gsea_obj, "gseaResult")) {
    .gene_universe <- if (!is.null(atlas_universe)) {
      as.character(atlas_universe)
    } else {
      names(gsea_obj@geneList)
    }
    .gene_sets     <- gsea_obj@geneSets
    result_df$genes_full_set <- vapply(
      as.character(result_df$ID),
      function(pid) {
        members <- .gene_sets[[pid]]
        if (is.null(members)) return(NA_character_)
        paste(intersect(members, .gene_universe), collapse = "/")
      },
      character(1)
    )
  } else {
    # data.frame input path: no @geneSets available; consumer should populate
    # genes_full_set upstream or accept NA (data_loader fallback applies).
    result_df$genes_full_set <- NA_character_
  }

  # Calculate gene ratio
  result_df$gene_ratio <- result_df$leading_edge_size / result_df$setSize
  result_df$gene_ratio[is.na(result_df$gene_ratio) | is.infinite(result_df$gene_ratio)] <- 0

  # Format pathway names if requested
  if (format_names) {
    # Try to use format_pathway_name if available
    if (exists("format_pathway_name", mode = "function")) {
      result_df$pathway_name <- format_pathway_name(
        result_df$Description,
        use_formatting = TRUE,
        strip_prefix = TRUE
      )
    } else {
      # Fallback: basic cleaning
      result_df$pathway_name <- clean_pathway_name_basic(
        result_df$Description,
        max_length = max_name_length
      )
    }
  } else {
    result_df$pathway_name <- as.character(result_df$Description)
  }

  # Truncate if still too long
  too_long <- nchar(result_df$pathway_name) > max_name_length
  result_df$pathway_name[too_long] <- paste0(
    substr(result_df$pathway_name[too_long], 1, max_name_length - 3), "..."
  )

  # Build standardized tibble
  normalized <- tibble::tibble(
    pathway_id = as.character(result_df$ID),
    pathway_name = as.character(result_df$pathway_name),
    database = database,
    contrast = contrast,
    NES = as.numeric(result_df$NES),
    pvalue = as.numeric(result_df$pvalue),
    padj = as.numeric(result_df[[padj_col]]),
    set_size = as.integer(result_df$setSize),
    leading_edge_size = as.integer(result_df$leading_edge_size),
    gene_ratio = as.numeric(result_df$gene_ratio),
    core_enrichment = as.character(result_df$core_enrichment),
    genes_full_set = as.character(result_df$genes_full_set),  # MADR-008
    direction = ifelse(result_df$NES > 0, "Up", "Down")
  )

  # Add computed columns
  normalized$neg_log_padj <- -log10(normalized$padj)
  normalized$neg_log_padj[is.infinite(normalized$neg_log_padj)] <- 16  # Cap at 16

  return(normalized)
}


#' Create Empty GSEA Tibble with Correct Schema
#'
#' Returns an empty tibble with the same column structure as normalize_gsea_results().
#' Useful for initializing or when no results are available.
#'
#' @return Empty tibble with standardized columns
#' @export
empty_gsea_tibble <- function() {
  tibble::tibble(
    pathway_id = character(),
    pathway_name = character(),
    database = character(),
    contrast = character(),
    NES = numeric(),
    pvalue = numeric(),
    padj = numeric(),
    set_size = integer(),
    leading_edge_size = integer(),
    gene_ratio = numeric(),
    core_enrichment = character(),
    genes_full_set = character(),  # MADR-008
    direction = character(),
    neg_log_padj = numeric()
  )
}


#' Basic Pathway Name Cleaning
#'
#' Simple pathway name cleaning without full smart capitalization.
#' Used as fallback when format_pathway_name is not available.
#'
#' @param names Character vector of pathway names
#' @param max_length Integer, truncate to this length
#' @param strip_prefix Logical, remove common database prefixes
#'
#' @return Cleaned character vector
#' @keywords internal
clean_pathway_name_basic <- function(names, max_length = 80, strip_prefix = TRUE) {
  cleaned <- as.character(names)

  if (strip_prefix) {
    # Common prefixes to remove
    prefixes <- c(
      "HALLMARK_", "KEGG_", "REACTOME_", "WP_",
      "GOBP_", "GOCC_", "GOMF_", "GO_",
      "MITOPATHWAYS_", "MITOXPLORER_", "MITOCARTA_",
      "TRANSPORTDB_", "GTRD_", "TFT_"
    )

    for (prefix in prefixes) {
      cleaned <- gsub(paste0("^", prefix), "", cleaned)
    }
  }

  # Replace underscores with spaces
  cleaned <- gsub("_", " ", cleaned)

  # Simple title case
  cleaned <- sapply(cleaned, function(x) {
    words <- strsplit(tolower(x), " ")[[1]]
    words <- sapply(words, function(w) {
      if (nchar(w) > 0) {
        paste0(toupper(substr(w, 1, 1)), substr(w, 2, nchar(w)))
      } else {
        w
      }
    })
    paste(words, collapse = " ")
  }, USE.NAMES = FALSE)

  return(cleaned)
}


#' Get Summary Statistics from Normalized GSEA Table
#'
#' Summarizes GSEA results by database, contrast, and direction.
#'
#' @param gsea_df Normalized GSEA tibble from normalize_gsea_results()
#' @param padj_cutoff FDR cutoff for significance (default: 0.05)
#'
#' @return Summary tibble with counts per database and direction
#' @export
summarize_gsea_results <- function(gsea_df, padj_cutoff = 0.05) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required for summarize_gsea_results()")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package required for summarize_gsea_results()")
  }

  gsea_df %>%
    dplyr::filter(padj < padj_cutoff) %>%
    dplyr::group_by(database, contrast, direction) %>%
    dplyr::summarise(
      n_pathways = dplyr::n(),
      mean_NES = mean(NES, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = direction,
      values_from = c(n_pathways, mean_NES),
      values_fill = list(n_pathways = 0, mean_NES = NA_real_)
    )
}


message("[RNAseq-toolkit] Loaded: normalize_gsea.R")
