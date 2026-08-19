name_audit_r_dir <- function() {
  for (path in c("../../R", "R", testthat::test_path("..", "..", "R"))) {
    if (dir.exists(path) && length(list.files(path, "[.]R$"))) return(path)
  }
  NULL
}

name_audit_exports <- function() {
  r_dir <- name_audit_r_dir()
  if (is.null(r_dir)) return(NULL)
  namespace <- file.path(dirname(r_dir), "NAMESPACE")
  if (!file.exists(namespace)) return(NULL)
  lines <- readLines(namespace, warn = FALSE)
  sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", lines, value = TRUE))
}

test_that("every live export has a layer prefix or a reasoned exception", {
  exports <- name_audit_exports()
  skip_if(is.null(exports), "NAMESPACE not reachable from the test dir")
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])

  top_level <- c(
    # Original I/O verbs describe the action more clearly than an io_ prefix.
    read_counts_matrix = "specific file-reading verb",
    read_metadata = "specific file-reading verb",
    ensure_dir = "frozen utility used directly by two unmigrated consumers",
    write_session_provenance = "specific provenance-writing verb",
    # These construct the DE object or its annotation before the de_ renderers.
    build_dge = "frozen constructor predating the de_ renderer layer",
    annotate_genes = "annotation verb used before differential expression",
    # Identifier conversion and filtering are cross-layer gene operations.
    gene_to_entrez = "identifier conversion is not owned by one layer",
    entrez_to_gene = "identifier conversion is not owned by one layer",
    filter_confounder_genes = "cross-layer gene filtering verb",
    # Presentation names follow R/ggplot vocabulary and remain intentionally short.
    format_pathway_name = "general label formatter",
    theme_bulki = "house theme parallel to ggplot2 theme names",
    # Package metadata already carries the package name as its namespace.
    bulkirna_api = "package metadata",
    bulkirna_check_deps = "package metadata",
    bulkirna_stochastic = "package metadata"
  )
  prefixed <- grepl("^(gsdb_|gs_|de_|gatom_|coresh_)", live)
  observed <- sort(live[!prefixed])

  expect_true(all(nzchar(top_level)))
  expect_identical(observed, sort(names(top_level)))
})

test_that("live export formals use the audited argument vocabulary", {
  exports <- name_audit_exports()
  skip_if(is.null(exports), "NAMESPACE not reachable from the test dir")
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])

  forbidden <- list(
    species = c("organism", "organism_name", "species_name"),
    db = c("gene_sets", "genesets", "gene_set_db"),
    contrast = c("comparison", "contrast_name"),
    seed = c("random_seed", "rng_seed"),
    quiet = "silent",
    verbose = "verbosity",
    path = c("file", "filepath", "file_path", "filename"),
    min_size = c("minSize", "minsize"),
    max_size = c("maxSize", "maxsize"),
    n_cores = c("cores", "ncores", "n_workers", "workers"),
    dir = c("directory", "dest_dir", "out_dir", "output_dir")
  )

  offenders <- unlist(lapply(live, function(name) {
    formal_names <- names(formals(getExportedValue("bulkiRNA", name)))
    bad <- intersect(formal_names, unlist(forbidden, use.names = FALSE))
    if (!length(bad)) return(character())
    paste0(name, "(", bad, ")")
  }), use.names = FALSE)

  expect_identical(
    sort(offenders),
    character(0L),
    info = paste(
      "Use species, db, contrast, seed, quiet, verbose, path, min_size,",
      "max_size, n_cores, or dir according to the audited concept."
    )
  )
})
