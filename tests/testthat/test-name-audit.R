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
  if (is.null(exports)) {
    skip(paste(
      "Source-tree-only name audit: R/ and NAMESPACE are unavailable",
      "when tests run against the installed package."
    ))
  }
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
    # The bulkirna_ prefix marks package-wide metadata and dependency checks.
    bulkirna_api = "package metadata",
    bulkirna_check_deps = "package-wide dependency inspection",
    bulkirna_stochastic = "package metadata"
  )
  prefixed <- grepl("^(gsdb_|gs_|de_|gatom_|coresh_)", live)
  observed <- sort(live[!prefixed])

  expect_identical(observed, sort(names(top_level)))
})

test_that("live export formals use the complete audited vocabulary", {
  exports <- name_audit_exports()
  if (is.null(exports)) {
    skip(paste(
      "Source-tree-only argument audit: R/ and NAMESPACE are unavailable",
      "when tests run against the installed package."
    ))
  }
  api <- bulkirna_api(quiet = TRUE)
  live <- intersect(exports, api$name[api$lifecycle != "deprecated"])

  allowed <- c(
    "...", "B_cutoff", "aes_x", "annotate_counts", "baseMean", "base_family",
    "base_size", "base_theme", "biomart_host", "biomart_version", "by",
    "by_contrast", "by_direction", "cache", "caption", "center", "chunk_dir",
    "chunk_path", "coef", "collapse", "collection", "color_palette",
    "colour_by", "colours", "compare", "contrast", "count_mat", "data",
    "database", "database_label", "database_labels", "db", "db_species", "de",
    "de_results", "decision_by", "df", "dge", "dir", "direction", "download",
    "dpi", "drop", "drop_empty", "ens_ids", "ensembl_version", "entity_type",
    "entrez", "eps", "error", "expr", "facet", "fc_cutoff", "fdr_cutoff",
    "features", "fit", "fixed_p_boundary", "formats", "gene2reaction_extra",
    "genes", "genes_df", "genome_build", "gpl", "grid", "gse_id", "gsea_param",
    "height", "highlight", "highlight_gene", "id", "input_gene_name",
    "jaccard_threshold", "k_gene", "k_met", "kcdf", "keep_first_caption",
    "label", "label_method", "label_size", "labels", "legend_pos",
    "legend_position", "lifecycle", "limits", "log2FC", "m", "max.overlaps",
    "max_genes", "max_name_length", "max_size", "met_de", "method", "metric",
    "metric_label", "min_genes", "min_queries", "min_size", "multi_vals", "n",
    "n_cores", "n_top", "name", "network", "networks", "norm_method", "obj",
    "orientation", "overview", "overwrite", "p_cutoff", "p_value", "padj",
    "padj_max", "palette", "panel_heights", "path", "pathway_id",
    "pathway_names", "pathways", "pattern", "per", "plot", "plots",
    "point_size", "prefix", "prune", "pval", "pvalues", "queries", "query",
    "quiet", "ranking", "ranks", "rebuild", "refs", "required_cols", "res",
    "round_nonint", "sample_col_candidates", "sample_data", "sample_size",
    "samples_df", "scale", "schema_version", "seed", "sets", "shape_by",
    "show_grid", "show_quadrant_counts", "size_range", "solver", "sort_by",
    "species", "stat", "stat_as_nes", "strip_prefix", "subcollection",
    "subtitle", "symbol_by", "symbols", "table", "text", "title", "top",
    "top_hits", "top_n", "unique_genes", "universe", "use_biomart",
    "use_formatting", "verbose", "width", "wrap_width", "x", "x_breaks",
    "xlim_abs", "y_padding", "ylim_abs"
  )
  required <- c(
    "species", "db", "contrast", "seed", "quiet", "verbose", "path",
    "min_size", "max_size", "n_cores", "dir"
  )
  observed <- sort(unique(unlist(lapply(live, function(name) {
    names(formals(getExportedValue("bulkiRNA", name)))
  }), use.names = FALSE)))
  unexpected <- setdiff(observed, allowed)
  unused <- setdiff(allowed, observed)

  expect_identical(
    unexpected,
    character(0L),
    info = "Every new formal spelling needs an explicit name-audit decision."
  )
  expect_identical(
    unused,
    character(0L),
    info = "Remove formal spellings from the vocabulary when they leave the API."
  )
  expect_true(
    all(required %in% observed),
    info = "The canonical cross-layer formals must remain present."
  )
})
