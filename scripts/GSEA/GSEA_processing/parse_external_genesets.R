#' Parse External Gene Set Databases
#'
#' Functions for parsing external gene set databases (TransportDB, MitoDB, etc.)
#' into clusterProfiler-compatible TERM2GENE/TERM2NAME format.
#'
#' @name parse_external_genesets
#' @keywords internal
NULL


#' Parse TransportDB CSV to Gene Sets
#'
#' Parses TransportDB2.0 CSV format into TERM2GENE/TERM2NAME data frames
#' for use with clusterProfiler::GSEA().
#'
#' TransportDB contains transporter proteins classified by family, subfamily,
#' and transport type. This function groups genes by family for GSEA.
#'
#' @param file Character, path to TransportDB CSV file
#' @param prefix Character, prefix for gene set names (default: "TRANSPORTDB")
#' @param group_by Character, column to group by: "family" (recommended), "subfamily", or "class" (default: "family")
#' @param org_db Annotation database for gene ID conversion (e.g., org.Mm.eg.db, org.Hs.eg.db)
#' @param id_type Character, input ID type in TransportDB (default: "SYMBOL")
#' @param min_size Integer, minimum genes per set to retain (default: 5)
#' @param max_size Integer, maximum genes per set to retain (default: 500)
#'
#' @return List with:
#'   - T2G: TERM2GENE data frame (gs_name, gene_symbol)
#'   - T2N: TERM2NAME data frame (gs_name, description)
#'   - stats: Gene set statistics
#'   - source: Database source info
#'
#' @examples
#' \dontrun{
#' library(org.Mm.eg.db)
#' transport_db <- parse_transportdb(
#'   "00_data/references/TransportDB2.0.csv",
#'   org_db = org.Mm.eg.db
#' )
#' }
#'
#' @export
parse_transportdb <- function(
    file,
    prefix = "TRANSPORTDB",
    group_by = "family",
    org_db = NULL,
    id_type = "SYMBOL",
    min_size = 5,
    max_size = 500
) {

  if (!file.exists(file)) {
    stop("TransportDB file not found: ", file)
  }

  message("Parsing TransportDB: ", basename(file))

  # Read CSV
  raw_data <- read.csv(file, stringsAsFactors = FALSE)

  # Expected columns vary by TransportDB version
  # Common columns: Gene, Symbol, Family, Subfamily, Class, Description
  # Try to detect columns
  gene_col <- intersect(c("Symbol", "Gene_Symbol", "gene_symbol", "SYMBOL", "Gene"), colnames(raw_data))[1]
  family_col <- intersect(c("Family", "family", "TC_Family"), colnames(raw_data))[1]
  desc_col <- intersect(c("Description", "description", "Family_Description"), colnames(raw_data))[1]

  if (is.na(gene_col)) {
    stop("Could not find gene symbol column in TransportDB file")
  }
  if (is.na(family_col)) {
    stop("Could not find family column in TransportDB file")
  }

  message(sprintf("  Using columns: gene='%s', family='%s'", gene_col, family_col))

  # Get unique families
  families <- unique(raw_data[[family_col]])
  families <- families[!is.na(families) & families != ""]

  message(sprintf("  Found %d transporter families", length(families)))

  # Build gene sets by family
  T2G_list <- list()
  T2N_list <- list()

  for (fam in families) {
    # Get genes for this family
    genes <- unique(raw_data[raw_data[[family_col]] == fam, gene_col])
    genes <- genes[!is.na(genes) & genes != ""]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      # Create standardized gene set name
      gs_name <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9]", "_", fam)))

      # Get description if available
      if (!is.na(desc_col)) {
        desc <- unique(raw_data[raw_data[[family_col]] == fam, desc_col])[1]
        if (is.na(desc)) desc <- fam
      } else {
        desc <- fam
      }

      T2G_list[[gs_name]] <- data.frame(
        gs_name = gs_name,
        gene_symbol = genes,
        stringsAsFactors = FALSE
      )

      T2N_list[[gs_name]] <- data.frame(
        gs_name = gs_name,
        description = desc,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine
  T2G <- do.call(rbind, T2G_list)
  T2N <- do.call(rbind, T2N_list)
  rownames(T2G) <- NULL
  rownames(T2N) <- NULL

  # Statistics
  stats <- list(
    n_genesets = length(unique(T2G$gs_name)),
    n_genes = length(unique(T2G$gene_symbol)),
    median_size = median(table(T2G$gs_name)),
    source_file = basename(file)
  )

  message(sprintf("  Created %d gene sets with %d unique genes",
                  stats$n_genesets, stats$n_genes))

  return(list(
    T2G = T2G,
    T2N = T2N,
    stats = stats,
    source = "TransportDB2.0",
    created = Sys.time()
  ))
}


#' Parse Generic TSV/CSV Gene Sets
#'
#' Parses a simple two-column or multi-column file into gene sets.
#' Flexible parser for custom gene set files.
#'
#' @param file Character, path to file
#' @param gs_col Character or integer, column name/index for gene set names
#' @param gene_col Character or integer, column name/index for gene symbols
#' @param desc_col Character or integer, column name/index for descriptions (optional)
#' @param prefix Character, prefix to add to gene set names (default: NULL)
#' @param sep Character, field separator (default: auto-detect)
#' @param header Logical, does file have header? (default: TRUE)
#'
#' @return List with T2G, T2N data frames
#' @export
parse_geneset_file <- function(
    file,
    gs_col = 1,
    gene_col = 2,
    desc_col = NULL,
    prefix = NULL,
    sep = NULL,
    header = TRUE
) {

  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

 # Auto-detect separator
  if (is.null(sep)) {
    first_line <- readLines(file, n = 1)
    if (grepl("\t", first_line)) {
      sep <- "\t"
    } else if (grepl(",", first_line)) {
      sep <- ","
    } else {
      sep <- "\t"  # Default
    }
  }

  # Read file
  data <- read.delim(file, sep = sep, header = header, stringsAsFactors = FALSE)

  # Resolve column names/indices
  if (is.character(gs_col)) gs_col <- which(colnames(data) == gs_col)
  if (is.character(gene_col)) gene_col <- which(colnames(data) == gene_col)
  if (!is.null(desc_col) && is.character(desc_col)) {
    desc_col <- which(colnames(data) == desc_col)
  }

  # Extract columns
  gs_names <- data[[gs_col]]
  genes <- data[[gene_col]]

  # Add prefix if specified
  if (!is.null(prefix)) {
    gs_names <- paste0(prefix, "_", gs_names)
  }

  # Build T2G
  T2G <- data.frame(
    gs_name = gs_names,
    gene_symbol = genes,
    stringsAsFactors = FALSE
  )

  # Remove NA/empty
  T2G <- T2G[!is.na(T2G$gs_name) & !is.na(T2G$gene_symbol), ]
  T2G <- T2G[T2G$gs_name != "" & T2G$gene_symbol != "", ]

  # Build T2N
  if (!is.null(desc_col) && length(desc_col) == 1 && desc_col <= ncol(data)) {
    descriptions <- data[[desc_col]]
    T2N <- unique(data.frame(
      gs_name = gs_names,
      description = descriptions,
      stringsAsFactors = FALSE
    ))
  } else {
    # Use gene set name as description
    T2N <- unique(data.frame(
      gs_name = gs_names,
      description = gsub("_", " ", gs_names),
      stringsAsFactors = FALSE
    ))
  }

  T2N <- T2N[!is.na(T2N$gs_name) & T2N$gs_name != "", ]

  message(sprintf("Parsed %d gene sets with %d genes from: %s",
                  length(unique(T2G$gs_name)),
                  length(unique(T2G$gene_symbol)),
                  basename(file)))

  return(list(
    T2G = T2G,
    T2N = T2N,
    source = basename(file),
    created = Sys.time()
  ))
}


#' Parse GMT File to Gene Sets
#'
#' Parses a GMT (Gene Matrix Transposed) file format commonly used by MSigDB
#' and GSEA tools.
#'
#' GMT format: gene_set_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...
#'
#' @param file Character, path to GMT file
#' @param prefix Character, prefix to add to gene set names (optional)
#'
#' @return List with T2G, T2N data frames
#' @export
parse_gmt <- function(file, prefix = NULL) {

  if (!file.exists(file)) {
    stop("GMT file not found: ", file)
  }

  lines <- readLines(file)

  T2G_list <- list()
  T2N_list <- list()

  for (line in lines) {
    if (nchar(line) == 0) next

    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next

    gs_name <- fields[1]
    description <- fields[2]
    genes <- fields[3:length(fields)]
    genes <- genes[genes != "" & !is.na(genes)]

    if (!is.null(prefix)) {
      gs_name <- paste0(prefix, "_", gs_name)
    }

    if (length(genes) > 0) {
      T2G_list[[gs_name]] <- data.frame(
        gs_name = gs_name,
        gene_symbol = genes,
        stringsAsFactors = FALSE
      )

      T2N_list[[gs_name]] <- data.frame(
        gs_name = gs_name,
        description = if (description == "" || description == "na") gs_name else description,
        stringsAsFactors = FALSE
      )
    }
  }

  T2G <- do.call(rbind, T2G_list)
  T2N <- do.call(rbind, T2N_list)
  rownames(T2G) <- NULL
  rownames(T2N) <- NULL

  message(sprintf("Parsed %d gene sets from GMT: %s",
                  length(unique(T2G$gs_name)), basename(file)))

  return(list(
    T2G = T2G,
    T2N = T2N,
    source = basename(file),
    created = Sys.time()
  ))
}


#' Convert Gene IDs in Gene Sets
#'
#' Converts gene identifiers in T2G from one type to another using
#' an annotation database.
#'
#' @param T2G Data frame with gene set-gene mappings
#' @param org_db Annotation database (e.g., org.Mm.eg.db)
#' @param from_type Input ID type (default: "ENSEMBL")
#' @param to_type Output ID type (default: "SYMBOL")
#' @param drop_unmapped Logical, remove genes that couldn't be mapped (default: TRUE)
#'
#' @return T2G with converted gene IDs
#' @export
convert_geneset_ids <- function(
    T2G,
    org_db,
    from_type = "ENSEMBL",
    to_type = "SYMBOL",
    drop_unmapped = TRUE
) {

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("AnnotationDbi package required for ID conversion")
  }

  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"
  original_ids <- unique(T2G[[gene_col]])

  # Clean IDs (remove version suffix for Ensembl)
  if (from_type == "ENSEMBL") {
    clean_ids <- gsub("\\..*$", "", original_ids)
  } else {
    clean_ids <- original_ids
  }

  # Map IDs
  mapped <- AnnotationDbi::mapIds(
    org_db,
    keys = clean_ids,
    column = to_type,
    keytype = from_type,
    multiVals = "first"
  )

  # Create mapping table
  id_map <- data.frame(
    original = original_ids,
    mapped = mapped[match(if (from_type == "ENSEMBL") gsub("\\..*$", "", original_ids) else original_ids, names(mapped))],
    stringsAsFactors = FALSE
  )

  # Report mapping stats
  n_mapped <- sum(!is.na(id_map$mapped))
  message(sprintf("ID conversion (%s -> %s): %d / %d mapped (%.1f%%)",
                  from_type, to_type, n_mapped, nrow(id_map),
                  100 * n_mapped / nrow(id_map)))

  # Apply mapping
  T2G$original_id <- T2G[[gene_col]]
  T2G[[gene_col]] <- id_map$mapped[match(T2G[[gene_col]], id_map$original)]

  # Drop unmapped if requested
  if (drop_unmapped) {
    n_before <- nrow(T2G)
    T2G <- T2G[!is.na(T2G[[gene_col]]), ]
    message(sprintf("  Dropped %d unmapped entries", n_before - nrow(T2G)))
  }

  return(T2G)
}


#' Parse GMX File to Gene Sets
#'
#' Parses a GMX (Gene Matrix Transposed) file format. GMX is the transpose of GMT:
#' - Row 1: descriptions (one per column)
#' - Row 2: gene set names (one per column)
#' - Row 3+: genes (one gene per row per column, variable length)
#'
#' @param file Character, path to GMX file
#' @param prefix Character, prefix to add to gene set names (optional)
#' @param min_size Integer, minimum genes per set to retain (default: 5)
#' @param max_size Integer, maximum genes per set to retain (default: 500)
#'
#' @return List with T2G, T2N data frames
#' @export
parse_gmx <- function(file, prefix = NULL, min_size = 5, max_size = 500) {

  if (!file.exists(file)) {
    stop("GMX file not found: ", file)
  }

  message("Parsing GMX file: ", basename(file))

  # Read all lines
  lines <- readLines(file, warn = FALSE)

  if (length(lines) < 3) {
    stop("GMX file must have at least 3 rows (descriptions, names, genes)")
  }

  # Parse header rows
  descriptions <- strsplit(lines[1], "\t")[[1]]
  gs_names <- strsplit(lines[2], "\t")[[1]]

  n_sets <- length(gs_names)
  message(sprintf("  Found %d gene sets in header", n_sets))

  # Parse gene rows - handle variable-length columns
  gene_rows <- lapply(lines[3:length(lines)], function(x) {
    strsplit(x, "\t")[[1]]
  })

  # Build gene sets by column
  T2G_list <- list()
  T2N_list <- list()

  for (i in seq_len(n_sets)) {
    gs_name <- gs_names[i]
    description <- if (i <= length(descriptions)) descriptions[i] else gs_name

    # Skip empty names
    if (is.na(gs_name) || gs_name == "") next

    # Collect genes for this column (index i)
    genes <- sapply(gene_rows, function(row) {
      if (i <= length(row)) row[i] else NA
    })
    genes <- genes[!is.na(genes) & genes != ""]

    # Apply size filters
    if (length(genes) < min_size || length(genes) > max_size) next

    # Apply prefix if specified
    if (!is.null(prefix)) {
      final_name <- paste0(prefix, "_", gs_name)
    } else {
      final_name <- gs_name
    }

    T2G_list[[final_name]] <- data.frame(
      gs_name = final_name,
      gene_symbol = genes,
      stringsAsFactors = FALSE
    )

    T2N_list[[final_name]] <- data.frame(
      gs_name = final_name,
      description = if (is.na(description) || description == "") final_name else description,
      stringsAsFactors = FALSE
    )
  }

  if (length(T2G_list) == 0) {
    warning("No gene sets passed size filters")
    return(list(
      T2G = data.frame(gs_name = character(), gene_symbol = character(), stringsAsFactors = FALSE),
      T2N = data.frame(gs_name = character(), description = character(), stringsAsFactors = FALSE),
      source = basename(file),
      created = Sys.time()
    ))
  }

  T2G <- do.call(rbind, T2G_list)
  T2N <- do.call(rbind, T2N_list)
  rownames(T2G) <- NULL
  rownames(T2N) <- NULL

  message(sprintf("  Parsed %d gene sets with %d unique genes",
                  length(unique(T2G$gs_name)),
                  length(unique(T2G$gene_symbol))))

  return(list(
    T2G = T2G,
    T2N = T2N,
    source = basename(file),
    created = Sys.time()
  ))
}


#' Parse mitoXplorer Gene Function File
#'
#' Parses mitoXplorer3.0 mouse_gene_function.txt format into gene sets
#' grouped by mito_process (mitochondrial function category).
#'
#' @param file Character, path to mitoXplorer mouse_gene_function.txt
#' @param prefix Character, prefix for gene set names (default: "MITOXPLORER")
#' @param gene_col Character, column name for gene symbols (default: "MGI_symbol")
#' @param process_col Character, column name for pathway/process (default: "mito_process")
#' @param min_size Integer, minimum genes per set to retain (default: 5)
#' @param max_size Integer, maximum genes per set to retain (default: 500)
#'
#' @return List with T2G, T2N data frames
#' @export
parse_mitoxplorer <- function(
    file,
    prefix = "MITOXPLORER",
    gene_col = "MGI_symbol",
    process_col = "mito_process",
    min_size = 5,
    max_size = 500
) {

  if (!file.exists(file)) {
    stop("mitoXplorer file not found: ", file)
  }

  message("Parsing mitoXplorer: ", basename(file))

  # Read tab-separated file
  data <- read.delim(file, stringsAsFactors = FALSE, header = TRUE)

  # Validate columns exist
  if (!gene_col %in% colnames(data)) {
    stop(sprintf("Gene column '%s' not found. Available: %s",
                 gene_col, paste(colnames(data), collapse = ", ")))
  }
  if (!process_col %in% colnames(data)) {
    stop(sprintf("Process column '%s' not found. Available: %s",
                 process_col, paste(colnames(data), collapse = ", ")))
  }

  # Get unique processes
  processes <- unique(data[[process_col]])
  processes <- processes[!is.na(processes) & processes != ""]

  message(sprintf("  Found %d mitochondrial processes", length(processes)))

  # Build gene sets by process
  T2G_list <- list()
  T2N_list <- list()

  for (proc in processes) {
    # Get genes for this process
    genes <- unique(data[data[[process_col]] == proc, gene_col])
    genes <- genes[!is.na(genes) & genes != ""]

    # Apply size filters
    if (length(genes) < min_size || length(genes) > max_size) next

    # Create standardized gene set name
    gs_name <- paste0(prefix, "_", toupper(gsub("[^A-Za-z0-9]", "_", proc)))

    T2G_list[[gs_name]] <- data.frame(
      gs_name = gs_name,
      gene_symbol = genes,
      stringsAsFactors = FALSE
    )

    T2N_list[[gs_name]] <- data.frame(
      gs_name = gs_name,
      description = proc,
      stringsAsFactors = FALSE
    )
  }

  if (length(T2G_list) == 0) {
    warning("No gene sets passed size filters")
    return(list(
      T2G = data.frame(gs_name = character(), gene_symbol = character(), stringsAsFactors = FALSE),
      T2N = data.frame(gs_name = character(), description = character(), stringsAsFactors = FALSE),
      source = basename(file),
      created = Sys.time()
    ))
  }

  T2G <- do.call(rbind, T2G_list)
  T2N <- do.call(rbind, T2N_list)
  rownames(T2G) <- NULL
  rownames(T2N) <- NULL

  # Stats
  message(sprintf("  Created %d gene sets with %d unique genes",
                  length(unique(T2G$gs_name)),
                  length(unique(T2G$gene_symbol))))

  return(list(
    T2G = T2G,
    T2N = T2N,
    stats = list(
      n_genesets = length(unique(T2G$gs_name)),
      n_genes = length(unique(T2G$gene_symbol)),
      source_file = basename(file)
    ),
    source = "mitoXplorer3.0",
    created = Sys.time()
  ))
}


#' Convert Human Gene Symbols to Mouse Orthologs
#'
#' Uses homologene database to convert human gene symbols to mouse orthologs.
#' Falls back to direct symbol matching for genes not in homologene.
#'
#' @param T2G Data frame with gene set-gene mappings (gs_name, gene_symbol columns)
#' @param drop_unmapped Logical, remove genes that couldn't be mapped (default: TRUE)
#' @param verbose Logical, print mapping statistics (default: TRUE)
#'
#' @return T2G with mouse gene symbols
#' @export
convert_human_to_mouse <- function(T2G, drop_unmapped = TRUE, verbose = TRUE) {

  if (!requireNamespace("homologene", quietly = TRUE)) {
    stop("Package 'homologene' required. Install with: install.packages('homologene')")
  }

  gene_col <- if ("gene_symbol" %in% colnames(T2G)) "gene_symbol" else "gene"
  human_genes <- unique(T2G[[gene_col]])

  if (verbose) message(sprintf("Converting %d human genes to mouse orthologs...", length(human_genes)))

  # Use homologene for conversion
  # human = 9606, mouse = 10090
  mouse_orthologs <- homologene::homologene(human_genes, inTax = 9606, outTax = 10090)

  # Build mapping table
  if (nrow(mouse_orthologs) > 0) {
    # homologene returns: [1] human symbol, [2] mouse symbol, [3] human ID, [4] mouse ID
    # Use column 2 (mouse symbol), NOT column 3 (human ID)
    id_map <- data.frame(
      human = mouse_orthologs[[1]],
      mouse = mouse_orthologs[[2]],  # Fixed: was [[3]] (human ID), now [[2]] (mouse symbol)
      stringsAsFactors = FALSE
    )
    # Remove duplicates (keep first mapping)
    id_map <- id_map[!duplicated(id_map$human), ]
  } else {
    id_map <- data.frame(human = character(), mouse = character(), stringsAsFactors = FALSE)
  }

  # Report mapping stats
  n_mapped <- sum(human_genes %in% id_map$human)
  if (verbose) {
    message(sprintf("  Mapped: %d / %d genes (%.1f%%)",
                    n_mapped, length(human_genes),
                    100 * n_mapped / length(human_genes)))
  }

  # Apply mapping
  T2G$original_human <- T2G[[gene_col]]
  T2G[[gene_col]] <- id_map$mouse[match(T2G[[gene_col]], id_map$human)]

  # Drop unmapped if requested
  if (drop_unmapped) {
    n_before <- nrow(T2G)
    T2G <- T2G[!is.na(T2G[[gene_col]]), ]
    if (verbose) message(sprintf("  Dropped %d unmapped entries", n_before - nrow(T2G)))
  }

  return(T2G)
}


message("[RNAseq-toolkit] Loaded: parse_external_genesets.R")
