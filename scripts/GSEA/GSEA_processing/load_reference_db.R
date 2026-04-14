#' Load Bundled Reference Databases
#'
#' Functions for loading pre-packaged reference gene set databases
#' that ship with the RNAseq-toolkit. Returns T2G/T2N lists compatible
#' with clusterProfiler::GSEA().
#'
#' @name load_reference_db
#' @keywords internal
NULL


# ===========================================================================
# Internal helpers
# ===========================================================================

#' Resolve toolkit root directory
#'
#' Detects the toolkit root from this file's location, or uses an
#' explicit path or environment variable.
#'
#' @param toolkit_dir Character, explicit path (optional)
#' @return Character, absolute path to toolkit root
#' @keywords internal
.resolve_toolkit_dir <- function(toolkit_dir = NULL) {

  if (!is.null(toolkit_dir)) {
    toolkit_dir <- normalizePath(toolkit_dir, mustWork = TRUE)
    return(toolkit_dir)
  }

  # Try environment variable
  env_dir <- Sys.getenv("RNASEQ_TOOLKIT_DIR", unset = "")
  if (nchar(env_dir) > 0 && dir.exists(env_dir)) {
    return(normalizePath(env_dir))
  }

  # Auto-detect from this file's location
  # This file is at: scripts/GSEA/GSEA_processing/load_reference_db.R
  # Toolkit root is 4 levels up
  this_file <- NULL

  # Try sys.frame approach (works when source()'d)
  for (i in seq_len(sys.nframe())) {
    f <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(f) && grepl("load_reference_db\\.R$", f)) {
      this_file <- normalizePath(f)
      break
    }
  }

  # Try srcref approach
  if (is.null(this_file)) {
    srcref <- tryCatch(
      attr(body(.resolve_toolkit_dir), "srcfile")$filename,
      error = function(e) NULL
    )
    if (!is.null(srcref) && file.exists(srcref)) {
      this_file <- normalizePath(srcref)
    }
  }

  if (!is.null(this_file)) {
    toolkit_root <- normalizePath(file.path(dirname(this_file), "..", "..", ".."))
    if (file.exists(file.path(toolkit_root, "data", "references", "METADATA.yaml"))) {
      return(toolkit_root)
    }
  }

  stop(
    "Cannot auto-detect RNAseq-toolkit directory.\n",
    "Please provide toolkit_dir explicitly, or set RNASEQ_TOOLKIT_DIR environment variable."
  )
}


#' Map database name to file paths
#' @keywords internal
.db_registry <- function() {
  list(
    mitopathways = list(
      dir       = "mitocarta3.0",
      rds_file  = "mito_mitopathways.rds",
      raw_file  = "raw/MitoPathways3.0.gmx",
      parser    = "parse_gmx",
      prefix    = "MITOPATHWAYS",
      convert   = TRUE,
      label     = "MitoPathways 3.0"
    ),
    mitoxplorer = list(
      dir       = "mitoxplorer3.0",
      rds_file  = "mito_mitoxplorer.rds",
      raw_file  = "raw/mouse_gene_function.txt",
      parser    = "parse_mitoxplorer",
      prefix    = "MITOXPLORER",
      convert   = FALSE,
      label     = "mitoXplorer 3.0"
    ),
    mito_unified = list(
      dir       = "mitochondria_unified",
      rds_file  = "unified_mito_pathways.rds",
      raw_file  = NULL,
      parser    = NULL,
      prefix    = NULL,
      convert   = FALSE,
      label     = "Unified Mitochondrial Pathways",
      is_composite = TRUE,
      components = c("mitopathways", "mitoxplorer")
    ),
    transportdb = list(
      dir       = "transportdb",
      rds_file  = "transportdb_genesets.rds",
      raw_file  = "raw/TransportDB2.0.csv",
      parser    = "parse_transportdb",
      prefix    = "TRANSPORTDB",
      convert   = FALSE,
      label     = "TransportDB 2.0"
    )
  )
}


# ===========================================================================
# Public API
# ===========================================================================

#' Load a Bundled Reference Database
#'
#' Loads a pre-packaged gene set database from the RNAseq-toolkit's
#' bundled references. Returns a T2G/T2N list compatible with
#' clusterProfiler::GSEA().
#'
#' @param database Character. One of: "mitopathways", "mitoxplorer",
#'   "mito_unified", "transportdb"
#' @param species Character. Species directory name (default: "Mus_musculus")
#' @param toolkit_dir Character. Path to RNAseq-toolkit root.
#'   Auto-detected if NULL.
#' @param rebuild Logical. If TRUE, re-parse from raw files even if
#'   processed RDS exists. Requires parser dependencies (default: FALSE).
#'
#' @return List with:
#'   \item{T2G}{data.frame with columns gs_name, gene_symbol}
#'   \item{T2N}{data.frame with columns gs_name, description}
#'   \item{source}{Character, database provenance string}
#'   \item{created}{POSIXct, build timestamp}
#'
#' @examples
#' \dontrun{
#' db <- load_reference_db("mito_unified")
#' gsea_result <- GSEA(ranked_genes, TERM2GENE = db$T2G, TERM2NAME = db$T2N)
#' }
#'
#' @export
load_reference_db <- function(
    database,
    species = "Mus_musculus",
    toolkit_dir = NULL,
    rebuild = FALSE
) {

  # Check for GATOM (not bundled)
  if (tolower(database) == "gatom") {
    stop(
      "GATOM network files are not bundled (too large, ~24MB).\n",
      "Use download_gatom_references() to fetch them:\n",
      "  download_gatom_references(dest_dir = '00_data/references/gatom')"
    )
  }

  # Validate database name
  registry <- .db_registry()
  if (!database %in% names(registry)) {
    available <- paste(names(registry), collapse = ", ")
    stop(
      sprintf("Unknown database: '%s'\n", database),
      sprintf("Available databases: %s\n", available),
      "For GATOM, use download_gatom_references() instead."
    )
  }

  toolkit_dir <- .resolve_toolkit_dir(toolkit_dir)
  db_info <- registry[[database]]
  refs_dir <- file.path(toolkit_dir, "data", "references")

  # Path to processed RDS
  rds_path <- file.path(refs_dir, db_info$dir, "processed", species, db_info$rds_file)

  # If processed file exists and no rebuild requested, load it

  if (file.exists(rds_path) && !rebuild) {
    message(sprintf("Loading %s (%s): %s", db_info$label, species, basename(rds_path)))
    result <- readRDS(rds_path)

    # Ensure standard structure
    if (!all(c("T2G", "T2N") %in% names(result))) {
      stop("Corrupted RDS file — missing T2G/T2N. Try rebuild = TRUE.")
    }

    return(result)
  }

  # Need to rebuild from raw files
  if (rebuild || !file.exists(rds_path)) {
    message(sprintf("Building %s for %s from raw files...", db_info$label, species))
    result <- .rebuild_database(database, species, toolkit_dir)

    # Save the result
    out_dir <- dirname(rds_path)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    saveRDS(result, rds_path)
    message(sprintf("  Saved: %s", rds_path))

    return(result)
  }

  stop(sprintf("Processed file not found for %s/%s and rebuild failed.", database, species))
}


#' Rebuild a Database from Raw Files
#' @keywords internal
.rebuild_database <- function(database, species, toolkit_dir) {

  registry <- .db_registry()
  db_info <- registry[[database]]
  refs_dir <- file.path(toolkit_dir, "data", "references")

  # Source parsers
  parser_dir <- file.path(toolkit_dir, "scripts", "GSEA", "GSEA_processing")
  source(file.path(parser_dir, "parse_external_genesets.R"), local = TRUE)

  # Handle composite databases
  if (isTRUE(db_info$is_composite)) {
    message("Building composite database from components...")
    components <- lapply(db_info$components, function(comp) {
      .rebuild_database(comp, species, toolkit_dir)
    })

    # Merge T2G and T2N
    T2G <- do.call(rbind, lapply(components, `[[`, "T2G"))
    T2N <- do.call(rbind, lapply(components, `[[`, "T2N"))

    # Deduplicate by Jaccard similarity
    T2G <- .deduplicate_genesets(T2G, T2N, threshold = 0.99)

    return(list(
      T2G = T2G$T2G,
      T2N = T2G$T2N,
      source = "Unified_Mitochondrial_Pathways",
      created = Sys.time()
    ))
  }

  # Regular database
  raw_path <- file.path(refs_dir, db_info$dir, db_info$raw_file)
  if (!file.exists(raw_path)) {
    stop(sprintf("Raw file not found: %s", raw_path))
  }

  # Call the appropriate parser
  result <- switch(db_info$parser,
    parse_gmx = parse_gmx(raw_path, prefix = db_info$prefix),
    parse_mitoxplorer = parse_mitoxplorer(raw_path, prefix = db_info$prefix),
    parse_transportdb = {
      org_db <- NULL
      if (grepl("musculus", species, ignore.case = TRUE)) {
        if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
          org_db <- org.Mm.eg.db::org.Mm.eg.db
        }
      } else if (grepl("sapiens", species, ignore.case = TRUE)) {
        if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
          org_db <- org.Hs.eg.db::org.Hs.eg.db
        }
      }
      parse_transportdb(raw_path, prefix = db_info$prefix, org_db = org_db)
    },
    stop(sprintf("Unknown parser: %s", db_info$parser))
  )

  # Apply species conversion if needed
  if (db_info$convert && grepl("musculus", species, ignore.case = TRUE)) {
    message("  Converting human symbols to mouse orthologs...")
    result$T2G <- convert_human_to_mouse(result$T2G)
  }

  return(result)
}


#' Deduplicate Gene Sets by Jaccard Similarity
#' @keywords internal
.deduplicate_genesets <- function(T2G, T2N, threshold = 0.99) {

  # Build gene set lists
  gs_list <- split(T2G$gene_symbol, T2G$gs_name)
  gs_names <- names(gs_list)
  n <- length(gs_names)

  if (n <= 1) {
    return(list(T2G = T2G, T2N = T2N))
  }

  # Find redundant pairs
  to_remove <- character()
  for (i in seq_len(n - 1)) {
    if (gs_names[i] %in% to_remove) next
    for (j in (i + 1):n) {
      if (gs_names[j] %in% to_remove) next

      set_i <- gs_list[[gs_names[i]]]
      set_j <- gs_list[[gs_names[j]]]
      jaccard <- length(intersect(set_i, set_j)) / length(union(set_i, set_j))

      if (jaccard >= threshold) {
        # Prefer MITOPATHWAYS over MITOXPLORER
        if (grepl("^MITOPATHWAYS_", gs_names[j])) {
          to_remove <- c(to_remove, gs_names[i])
        } else {
          to_remove <- c(to_remove, gs_names[j])
        }
      }
    }
  }

  if (length(to_remove) > 0) {
    message(sprintf("  Removed %d redundant gene sets (Jaccard >= %.2f)",
                    length(to_remove), threshold))
    T2G <- T2G[!T2G$gs_name %in% to_remove, ]
    T2N <- T2N[!T2N$gs_name %in% to_remove, ]
  }

  return(list(T2G = T2G, T2N = T2N))
}


#' List Available Reference Databases
#'
#' Returns a data frame of all reference databases available in the toolkit,
#' both bundled and external (GATOM).
#'
#' @param toolkit_dir Character. Path to RNAseq-toolkit root.
#'   Auto-detected if NULL.
#'
#' @return Data frame with columns: database, name, bundled, species, description
#'
#' @examples
#' \dontrun{
#' list_reference_dbs()
#' }
#'
#' @export
list_reference_dbs <- function(toolkit_dir = NULL) {

  toolkit_dir <- .resolve_toolkit_dir(toolkit_dir)
  metadata_path <- file.path(toolkit_dir, "data", "references", "METADATA.yaml")

  if (!file.exists(metadata_path)) {
    stop("METADATA.yaml not found at: ", metadata_path)
  }

  if (!requireNamespace("yaml", quietly = TRUE)) {
    # Fallback: use built-in registry
    registry <- .db_registry()
    return(data.frame(
      database = names(registry),
      name = sapply(registry, `[[`, "label"),
      bundled = TRUE,
      stringsAsFactors = FALSE
    ))
  }

  meta <- yaml::read_yaml(metadata_path)

  rows <- lapply(names(meta$databases), function(db_name) {
    db <- meta$databases[[db_name]]
    data.frame(
      database    = db_name,
      name        = db$name %||% db_name,
      bundled     = isTRUE(db$bundled),
      description = db$description %||% "",
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  # List available species for each bundled database
  refs_dir <- file.path(toolkit_dir, "data", "references")
  result$species <- sapply(seq_len(nrow(result)), function(i) {
    if (!result$bundled[i]) return("(not bundled)")
    db <- meta$databases[[result$database[i]]]
    proc_dir <- file.path(refs_dir, db$directory %||% result$database[i], "processed")
    if (dir.exists(proc_dir)) {
      paste(list.dirs(proc_dir, full.names = FALSE, recursive = FALSE), collapse = ", ")
    } else {
      "(none)"
    }
  })

  return(result)
}


#' Get Reference Database Metadata
#'
#' Returns full metadata for a single reference database, including
#' citation information.
#'
#' @param database Character. Database name.
#' @param toolkit_dir Character. Path to RNAseq-toolkit root.
#'
#' @return List with metadata fields from METADATA.yaml
#'
#' @export
get_reference_db_info <- function(database, toolkit_dir = NULL) {

  toolkit_dir <- .resolve_toolkit_dir(toolkit_dir)
  metadata_path <- file.path(toolkit_dir, "data", "references", "METADATA.yaml")

  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' required. Install with: install.packages('yaml')")
  }

  meta <- yaml::read_yaml(metadata_path)

  if (!database %in% names(meta$databases)) {
    stop(sprintf("Unknown database: '%s'. Available: %s",
                 database, paste(names(meta$databases), collapse = ", ")))
  }

  db_info <- meta$databases[[database]]

  # Add citation file path
  db_dir <- db_info$directory %||% database
  bib_path <- file.path(toolkit_dir, "data", "references", db_dir, "CITATIONS.bib")
  if (file.exists(bib_path)) {
    db_info$citations_path <- bib_path
    db_info$citations_text <- readLines(bib_path, warn = FALSE)
  }

  return(db_info)
}


#' Download GATOM Reference Files
#'
#' Downloads GATOM network files from the Artyomov Lab server.
#' These files are too large to bundle with the toolkit (~24MB).
#'
#' @param dest_dir Character. Destination directory (default: "00_data/references/gatom")
#' @param species Character. "Mus_musculus" or "Homo_sapiens" (default: "Mus_musculus")
#' @param networks Character vector. Networks to download (default: c("kegg", "combined"))
#' @param overwrite Logical. Overwrite existing files (default: FALSE)
#'
#' @return Character vector of downloaded file paths (invisibly)
#'
#' @examples
#' \dontrun{
#' download_gatom_references(dest_dir = "00_data/references/gatom")
#' }
#'
#' @export
download_gatom_references <- function(
    dest_dir = "00_data/references/gatom",
    species = "Mus_musculus",
    networks = c("kegg", "combined"),
    overwrite = FALSE
) {

  base_url <- "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM"

  # Species codes
  species_code <- switch(species,
    Mus_musculus  = list(short = "Mm", ncbi = "mmu"),
    Homo_sapiens  = list(short = "Hs", ncbi = "hsa"),
    stop("Unsupported species: ", species, ". Use 'Mus_musculus' or 'Homo_sapiens'.")
  )

  # Build file list
  files_to_download <- c()

  for (net in networks) {
    files_to_download <- c(files_to_download,
      sprintf("network.%s.rds", net),
      sprintf("met.%s.db.rds", net),
      sprintf("gene2reaction.%s.%s.eg.tsv", net, species_code$ncbi)
    )
  }

  # Always download the organism annotation
  files_to_download <- c(files_to_download,
    sprintf("org.%s.eg.gatom.anno.rds", species_code$short)
  )

  files_to_download <- unique(files_to_download)

  # Create destination directory
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
    message("Created directory: ", dest_dir)
  }

  # Download each file
  downloaded <- character()
  for (fname in files_to_download) {
    dest_file <- file.path(dest_dir, fname)

    if (file.exists(dest_file) && !overwrite) {
      message(sprintf("  [skip] %s (exists, use overwrite=TRUE to replace)", fname))
      downloaded <- c(downloaded, dest_file)
      next
    }

    url <- file.path(base_url, fname)
    message(sprintf("  Downloading %s ...", fname))

    tryCatch({
      download.file(url, dest_file, mode = "wb", quiet = TRUE)
      downloaded <- c(downloaded, dest_file)
      message(sprintf("  [ok] %s (%.1f MB)",
                      fname, file.info(dest_file)$size / 1e6))
    }, error = function(e) {
      warning(sprintf("  [FAIL] %s: %s", fname, conditionMessage(e)))
    })
  }

  message(sprintf("\nDownloaded %d / %d files to: %s",
                  length(downloaded), length(files_to_download), dest_dir))

  invisible(downloaded)
}


# Null-coalesce operator (if not already defined)
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}
