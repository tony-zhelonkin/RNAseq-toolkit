#!/usr/bin/env Rscript
#' Tests for load_reference_db.R
#'
#' Run with: Rscript tests/test_load_reference_db.R

# Resolve toolkit root (this file is at tests/test_load_reference_db.R)
.get_script_dir <- function() {
  # Try --file= argument (Rscript)
  args <- commandArgs(FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg[1]))))
  }
  # Fallback: working directory
  return(getwd())
}
toolkit_dir <- normalizePath(file.path(.get_script_dir(), ".."))

`%||%` <- function(x, y) if (is.null(x)) y else x

message("=== Testing load_reference_db.R ===")
message("Toolkit dir: ", toolkit_dir)

# Source the loader
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))

passed <- 0
failed <- 0

test <- function(name, expr) {
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      message(sprintf("  [PASS] %s", name))
      passed <<- passed + 1
    } else {
      message(sprintf("  [FAIL] %s — returned: %s", name, as.character(result)))
      failed <<- failed + 1
    }
  }, error = function(e) {
    message(sprintf("  [FAIL] %s — error: %s", name, conditionMessage(e)))
    failed <<- failed + 1
  })
}


# ---- Test list_reference_dbs() ----------------------------------------------
message("\n--- list_reference_dbs() ---")

dbs <- list_reference_dbs(toolkit_dir = toolkit_dir)
test("Returns data frame", is.data.frame(dbs))
test("Has expected columns", all(c("database", "name", "bundled") %in% colnames(dbs)))
test("Contains mitopathways", "mitopathways" %in% dbs$database)
test("Contains mitoxplorer", "mitoxplorer" %in% dbs$database)
test("Contains mito_unified", "mito_unified" %in% dbs$database)
test("Contains transportdb", "transportdb" %in% dbs$database)
test("Contains gatom", "gatom" %in% dbs$database)
test("GATOM is not bundled", !dbs$bundled[dbs$database == "gatom"])


# ---- Test load_reference_db() — mito_unified --------------------------------
message("\n--- load_reference_db('mito_unified') ---")

db <- load_reference_db("mito_unified", toolkit_dir = toolkit_dir)
test("Returns list", is.list(db))
test("Has T2G", "T2G" %in% names(db))
test("Has T2N", "T2N" %in% names(db))
test("Has source", "source" %in% names(db))
test("T2G is data.frame", is.data.frame(db$T2G))
test("T2G has gs_name column", "gs_name" %in% colnames(db$T2G))
test("T2G has gene_symbol column", "gene_symbol" %in% colnames(db$T2G))
test("T2N has gs_name column", "gs_name" %in% colnames(db$T2N))
test("T2N has description column", "description" %in% colnames(db$T2N))
test("Has >5 gene sets", length(unique(db$T2G$gs_name)) > 5)
test("Has >50 unique genes", length(unique(db$T2G$gene_symbol)) > 50)


# ---- Test load_reference_db() — mitopathways ---------------------------------
message("\n--- load_reference_db('mitopathways') ---")

db_mp <- load_reference_db("mitopathways", toolkit_dir = toolkit_dir)
test("Returns list", is.list(db_mp))
test("T2G has correct columns", all(c("gs_name", "gene_symbol") %in% colnames(db_mp$T2G)))
test("Gene set names have MITOPATHWAYS prefix",
     all(grepl("^MITOPATHWAYS_", unique(db_mp$T2G$gs_name))))
test("Has >10 gene sets", length(unique(db_mp$T2G$gs_name)) > 10)


# ---- Test load_reference_db() — mitoxplorer ----------------------------------
message("\n--- load_reference_db('mitoxplorer') ---")

db_mx <- load_reference_db("mitoxplorer", toolkit_dir = toolkit_dir)
test("Returns list", is.list(db_mx))
test("Gene set names have MITOXPLORER prefix",
     all(grepl("^MITOXPLORER_", unique(db_mx$T2G$gs_name))))


# ---- Test load_reference_db() — transportdb ----------------------------------
message("\n--- load_reference_db('transportdb') ---")

db_tdb <- load_reference_db("transportdb", toolkit_dir = toolkit_dir)
test("Returns list", is.list(db_tdb))
test("Gene set names have TRANSPORTDB prefix",
     all(grepl("^TRANSPORTDB_", unique(db_tdb$T2G$gs_name))))
test("Has >5 gene sets", length(unique(db_tdb$T2G$gs_name)) > 5)


# ---- Test error handling — nonexistent database ------------------------------
message("\n--- Error handling ---")

test("Unknown database errors", tryCatch({
  load_reference_db("nonexistent", toolkit_dir = toolkit_dir)
  FALSE
}, error = function(e) grepl("Unknown database", conditionMessage(e))))

test("GATOM errors with download hint", tryCatch({
  load_reference_db("gatom", toolkit_dir = toolkit_dir)
  FALSE
}, error = function(e) grepl("download_gatom_references", conditionMessage(e))))


# ---- Test get_reference_db_info() -------------------------------------------
message("\n--- get_reference_db_info() ---")

info <- get_reference_db_info("mitopathways", toolkit_dir = toolkit_dir)
test("Returns list", is.list(info))
test("Has name field", !is.null(info$name))
test("Has source_url", !is.null(info$source_url))
test("Has citations path", !is.null(info$citations_path))
test("Citations file exists", file.exists(info$citations_path))


# ---- Summary -----------------------------------------------------------------
message(sprintf("\n=== Results: %d passed, %d failed ===", passed, failed))
if (failed > 0) {
  quit(status = 1)
} else {
  message("All tests passed!")
}
