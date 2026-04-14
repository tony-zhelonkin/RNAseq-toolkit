#!/usr/bin/env Rscript
#' Build Reference Databases from Raw Files
#'
#' Parses raw gene set files bundled with the RNAseq-toolkit and generates
#' processed RDS files for a given species. Idempotent — safe to re-run.
#'
#' Usage:
#'   Rscript build_reference_databases.R [species]
#'
#' Arguments:
#'   species  - "Mus_musculus" (default) or "Homo_sapiens"
#'
#' Output:
#'   data/references/{db}/processed/{species}/*.rds
#'
#' @keywords internal

# ---- Parse command-line args ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
species <- if (length(args) >= 1) args[1] else "Mus_musculus"
message(sprintf("=== Building reference databases for: %s ===\n", species))

# ---- Resolve paths ----------------------------------------------------------
# This script lives at: scripts/GSEA/GSEA_processing/build_reference_databases.R
script_dir <- dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE))
if (script_dir == ".") {
  # Fallback: try to detect from command line
  script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
}
toolkit_dir <- normalizePath(file.path(script_dir, "..", "..", ".."))
refs_dir <- file.path(toolkit_dir, "data", "references")
parser_file <- file.path(script_dir, "parse_external_genesets.R")

message("Toolkit root: ", toolkit_dir)
message("References:   ", refs_dir)
message("Parsers:      ", parser_file)

if (!file.exists(parser_file)) {
  stop("Cannot find parse_external_genesets.R at: ", parser_file)
}

# ---- Source parsers ---------------------------------------------------------
source(parser_file)

# Null-coalesce
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Helper: ensure output directory ----------------------------------------
ensure_dir <- function(path) {
  dir <- dirname(path)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# ---- Helper: Jaccard deduplication ------------------------------------------
deduplicate_genesets <- function(T2G, T2N, threshold = 0.99) {
  gs_list <- split(T2G$gene_symbol, T2G$gs_name)
  gs_names <- names(gs_list)
  n <- length(gs_names)

  if (n <= 1) return(list(T2G = T2G, T2N = T2N))

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

  list(T2G = T2G, T2N = T2N)
}


# ---- Helper: filter by size -------------------------------------------------
filter_by_size <- function(result, min_size = 5, max_size = 500) {
  gs_sizes <- table(result$T2G$gs_name)
  keep <- names(gs_sizes[gs_sizes >= min_size & gs_sizes <= max_size])
  n_before <- length(unique(result$T2G$gs_name))

  result$T2G <- result$T2G[result$T2G$gs_name %in% keep, ]
  result$T2N <- result$T2N[result$T2N$gs_name %in% keep, ]

  n_after <- length(keep)
  if (n_before != n_after) {
    message(sprintf("  Size filter (%d-%d genes): %d -> %d gene sets",
                    min_size, max_size, n_before, n_after))
  }

  result
}


# =============================================================================
# Build each database
# =============================================================================

results <- list()

# ---- 1. MitoPathways 3.0 ---------------------------------------------------
message("\n--- MitoPathways 3.0 ---")
mp_raw <- file.path(refs_dir, "mitocarta3.0", "raw", "MitoPathways3.0.gmx")

if (file.exists(mp_raw)) {
  mp_result <- parse_gmx(mp_raw, prefix = "MITOPATHWAYS")

  # Species conversion
  if (grepl("musculus", species, ignore.case = TRUE)) {
    message("  Converting human -> mouse via homologene...")
    mp_result$T2G <- convert_human_to_mouse(mp_result$T2G)
  }
  # If human, gene symbols are already human — no conversion needed

  mp_result <- filter_by_size(mp_result)
  mp_result$source <- "MitoPathways3.0"
  mp_result$created <- Sys.time()

  out_path <- file.path(refs_dir, "mitocarta3.0", "processed", species, "mito_mitopathways.rds")
  ensure_dir(out_path)
  saveRDS(mp_result, out_path)
  message(sprintf("  Saved: %s (%d gene sets, %d genes)",
                  basename(out_path),
                  length(unique(mp_result$T2G$gs_name)),
                  length(unique(mp_result$T2G$gene_symbol))))
  results$mitopathways <- mp_result
} else {
  message("  [SKIP] Raw file not found: ", mp_raw)
}


# ---- 2. mitoXplorer 3.0 ----------------------------------------------------
message("\n--- mitoXplorer 3.0 ---")
mx_raw <- file.path(refs_dir, "mitoxplorer3.0", "raw", "mouse_gene_function.txt")

if (grepl("sapiens", species, ignore.case = TRUE)) {
  message("  [SKIP] mitoXplorer is mouse-only; no human data available")
} else if (file.exists(mx_raw)) {
  mx_result <- parse_mitoxplorer(mx_raw, prefix = "MITOXPLORER")
  mx_result <- filter_by_size(mx_result)
  mx_result$source <- "mitoXplorer3.0"
  mx_result$created <- Sys.time()

  out_path <- file.path(refs_dir, "mitoxplorer3.0", "processed", species, "mito_mitoxplorer.rds")
  ensure_dir(out_path)
  saveRDS(mx_result, out_path)
  message(sprintf("  Saved: %s (%d gene sets, %d genes)",
                  basename(out_path),
                  length(unique(mx_result$T2G$gs_name)),
                  length(unique(mx_result$T2G$gene_symbol))))
  results$mitoxplorer <- mx_result
} else {
  message("  [SKIP] Raw file not found: ", mx_raw)
}


# ---- 3. Unified Mitochondrial Pathways --------------------------------------
message("\n--- Unified Mitochondrial Pathways ---")

if (!is.null(results$mitopathways)) {
  # Merge available components
  components <- list(results$mitopathways)
  if (!is.null(results$mitoxplorer)) {
    components <- c(components, list(results$mitoxplorer))
  }

  unified_T2G <- do.call(rbind, lapply(components, `[[`, "T2G"))
  unified_T2N <- do.call(rbind, lapply(components, `[[`, "T2N"))

  # Deduplicate
  deduped <- deduplicate_genesets(unified_T2G, unified_T2N, threshold = 0.99)

  unified_result <- list(
    T2G = deduped$T2G,
    T2N = deduped$T2N,
    source = "Unified_Mitochondrial_Pathways",
    created = Sys.time()
  )

  out_path <- file.path(refs_dir, "mitochondria_unified", "processed", species,
                        "unified_mito_pathways.rds")
  ensure_dir(out_path)
  saveRDS(unified_result, out_path)
  message(sprintf("  Saved: %s (%d gene sets, %d genes)",
                  basename(out_path),
                  length(unique(unified_result$T2G$gs_name)),
                  length(unique(unified_result$T2G$gene_symbol))))
  results$mito_unified <- unified_result
} else {
  message("  [SKIP] No mitopathways result available to merge")
}


# ---- 4. TransportDB 2.0 ----------------------------------------------------
message("\n--- TransportDB 2.0 ---")
tdb_raw <- file.path(refs_dir, "transportdb", "raw", "TransportDB2.0.csv")

if (file.exists(tdb_raw)) {
  # Load species-appropriate annotation DB
  org_db <- NULL
  if (grepl("musculus", species, ignore.case = TRUE)) {
    if (requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      org_db <- org.Mm.eg.db::org.Mm.eg.db
    } else {
      message("  [WARN] org.Mm.eg.db not available; gene ID conversion may be limited")
    }
  } else if (grepl("sapiens", species, ignore.case = TRUE)) {
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      org_db <- org.Hs.eg.db::org.Hs.eg.db
    } else {
      message("  [WARN] org.Hs.eg.db not available; gene ID conversion may be limited")
    }
  }

  tdb_result <- parse_transportdb(tdb_raw, prefix = "TRANSPORTDB", org_db = org_db)
  tdb_result <- filter_by_size(tdb_result)
  tdb_result$source <- "TransportDB2.0"
  tdb_result$created <- Sys.time()

  out_path <- file.path(refs_dir, "transportdb", "processed", species, "transportdb_genesets.rds")
  ensure_dir(out_path)
  saveRDS(tdb_result, out_path)
  message(sprintf("  Saved: %s (%d gene sets, %d genes)",
                  basename(out_path),
                  length(unique(tdb_result$T2G$gs_name)),
                  length(unique(tdb_result$T2G$gene_symbol))))
  results$transportdb <- tdb_result
} else {
  message("  [SKIP] Raw file not found: ", tdb_raw)
}


# =============================================================================
# Summary
# =============================================================================
message("\n=== Build Summary ===")
message(sprintf("Species: %s", species))
for (db_name in names(results)) {
  r <- results[[db_name]]
  message(sprintf("  %-20s: %3d gene sets, %5d unique genes",
                  db_name,
                  length(unique(r$T2G$gs_name)),
                  length(unique(r$T2G$gene_symbol))))
}
message("\nDone.")
