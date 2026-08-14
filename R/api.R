#' Public API stability registry
#'
#' `bulkirna_api()` is the machine-readable contract for the package's public
#' surface. Stability and the historical signature freeze are separate axes:
#'
#' * `stable` functions follow semantic versioning. Incompatible changes need
#'   a major release, apart from removals that completed a documented
#'   deprecation cycle.
#' * `experimental` functions may change or be removed without a major release.
#' * `deprecated` functions remain callable and warn, but are scheduled for
#'   removal in the version recorded in `removed_in`.
#'
#' `frozen` records the 24 signatures carried forward from the script-library
#' API. It does not mean that a function is recommended: all 20 deprecated
#' shims are frozen until their documented removal. Conversely, most stable
#' functions were introduced after the freeze and therefore have
#' `frozen = FALSE`.
#'
#' The `superseded_by` field may name a technique or a sequence of calls rather
#' than a single function, and may say that no replacement exists when part of
#' the deprecated behaviour was deliberately retired.
#'
#' @param lifecycle `"all"`, or one or more of `"stable"`, `"experimental"`,
#'   and `"deprecated"`.
#' @param quiet Logical. If `FALSE`, print the registry and return it
#'   invisibly. If `TRUE`, do not print and return it visibly.
#' @return A tibble with one row per selected exported function and columns `name`,
#'   `layer`, `lifecycle`, `frozen`, `superseded_by`, and `removed_in`.
#' @examples
#' experimental <- bulkirna_api("experimental", quiet = TRUE)
#' @export
bulkirna_api <- function(lifecycle = "all", quiet = FALSE) {
  lifecycle_requested <- match.arg(
    lifecycle,
    choices = c("all", "stable", "experimental", "deprecated"),
    several.ok = TRUE
  )
  if (!is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) {
    stop("`quiet` must be a single non-missing logical value.", call. = FALSE)
  }

  out <- .bulkirna_api_registry()
  if (!"all" %in% lifecycle_requested) {
    out <- out[out$lifecycle %in% lifecycle_requested, , drop = FALSE]
  }
  out <- tibble::as_tibble(out)

  # n = Inf: a truncated report is worse than none, since the rows it hides
  # are exactly the ones someone is checking for.
  if (!quiet) print(out, n = Inf)
  if (!quiet) return(invisible(out))
  out
}

#' Registry of public API metadata
#'
#' Keeps the layer, lifecycle, signature-freeze, and deprecation metadata used
#' by [bulkirna_api()] in one place.
#'
#' @return A data frame with one row per exported function.
#' @keywords internal
.bulkirna_api_registry <- function() {
  stable <- list(
    gsdb = c(
      "gsdb_from_file", "gsdb_info", "gsdb_list", "gsdb_load",
      "gsdb_msigdb", "gsdb_register"
    ),
    gs = c(
      "gs_filter", "gs_leading_edge", "gs_plot_bar", "gs_plot_dot",
      "gs_plot_heatmap", "gs_plot_running", "gs_ranks", "gs_read",
      "gs_save", "gs_score", "gs_split", "gs_stat_types", "gs_test",
      "gs_to_master", "gs_top", "gs_validate_master", "gs_write"
    ),
    de = c(
      "de_bfc_plot", "de_md_plot", "de_pca", "de_pca_3d", "de_volcano",
      "de_volcano_grid"
    ),
    gatom = c(
      "download_gatom_references", "gatom_de", "gatom_genes",
      "gatom_module", "gatom_refs", "gatom_save_html"
    ),
    `top-level` = c(
      "annotate_genes", "build_dge", "bulkirna_api",
      "bulkirna_check_deps", "ensure_dir", "format_pathway_name",
      "read_counts_matrix", "read_metadata", "theme_bulki",
      "write_session_provenance"
    )
  )

  experimental <- list(
    coresh = c(
      "coresh_chunks", "coresh_convergence", "coresh_match",
      "coresh_search", "coresh_validate"
    ),
    `top-level` = c(
      "entrez_to_gene", "filter_confounder_genes", "gene_to_entrez"
    )
  )

  deprecated <- c(
    "convert_human_to_mouse", "create_MD_plot", "create_standard_volcano",
    "custom_minimal_theme_with_grid", "empty_gsea_tibble", "filter_by_size",
    "gsea_barplot", "gsea_dotplot", "gsea_dotplot_facet",
    "gsea_running_sum_plot", "list_reference_dbs", "list_to_term2gene",
    "load_reference_db", "normalize_gsea_results", "parse_gmx",
    "parse_mitoxplorer", "plot_all_gsea_results", "run_gsea",
    "run_gsea_analysis", "save_gsea_log"
  )

  stable_names <- unlist(stable, use.names = FALSE)
  experimental_names <- unlist(experimental, use.names = FALSE)
  names_all <- c(stable_names, experimental_names, deprecated)

  layer <- sub("_.*$", "", names_all)
  known_layers <- c("gs", "gsdb", "de", "gatom", "coresh", "gsea")
  layer[!layer %in% known_layers] <- "top-level"
  lifecycle <- c(
    rep("stable", length(stable_names)),
    rep("experimental", length(experimental_names)),
    rep("deprecated", length(deprecated))
  )

  frozen_names <- c(
    "normalize_gsea_results", "run_gsea", "create_standard_volcano",
    "format_pathway_name", "gsea_running_sum_plot", "list_to_term2gene",
    "gsea_barplot", "gsea_dotplot", "load_reference_db",
    "custom_minimal_theme_with_grid", "gsea_dotplot_facet", "create_MD_plot",
    "empty_gsea_tibble", "ensure_dir", "run_gsea_analysis", "save_gsea_log",
    "plot_all_gsea_results", "convert_human_to_mouse", "parse_gmx",
    "parse_mitoxplorer", "filter_by_size", "build_dge",
    "list_reference_dbs", "download_gatom_references"
  )

  successors <- vapply(
    deprecated, .bulkirna_deprecation_target, character(1L)
  )
  message_only <- c(
    "filter_by_size", "parse_mitoxplorer", "convert_human_to_mouse",
    "empty_gsea_tibble", "plot_all_gsea_results", "save_gsea_log"
  )
  stopifnot(all(is.na(successors[message_only])))
  successors[message_only] <- c(
    "the min_size/max_size arguments on gsdb_msigdb(), gsdb_load() and gsdb_from_file()",
    "gsdb_load(\"mitoxplorer\"), or gsdb_from_file() for an arbitrary file",
    "gsdb_msigdb(species = \"Mus musculus\", db_species = \"HS\")",
    paste0(
      "gs_test(); filter its result to zero rows for an empty gs_result; ",
      "there is no exported constructor, by design"
    ),
    "gs_plot_dot(), gs_plot_bar(), gs_plot_running() and gs_save()",
    "gs_save(); the free-text log has no replacement, by design"
  )

  superseded_by <- rep(NA_character_, length(names_all))
  names(superseded_by) <- names_all
  superseded_by[deprecated] <- successors[deprecated]

  removed_in <- rep(NA_character_, length(names_all))
  removed_in[names_all %in% deprecated] <- "1.0.0"

  out <- data.frame(
    name = names_all,
    layer = unname(layer),
    lifecycle = lifecycle,
    frozen = names_all %in% frozen_names,
    superseded_by = unname(superseded_by),
    removed_in = removed_in,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$name), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Read a successor from a deprecation shim
#'
#' The first expression in every legacy shim is its `.Deprecated()` call. For
#' shims that use `.Deprecated(new =)` (including an unnamed first argument),
#' evaluate that argument so the API registry cannot drift from the warning.
#' Message-only shims deliberately return `NA_character_`; their prose-level
#' successor is supplied by `bulkirna_api()`.
#'
#' @param name Name of a deprecated export.
#' @return A character scalar or `NA_character_`.
#' @keywords internal
.bulkirna_deprecation_target <- function(name) {
  namespace <- environment(.bulkirna_deprecation_target)
  fun <- get(name, envir = namespace, mode = "function", inherits = FALSE)
  expr <- body(fun)
  if (!is.call(expr) || !identical(expr[[1L]], as.name("{")) ||
      length(expr) < 2L) {
    stop("Deprecated shim `", name, "` has an unexpected body.", call. = FALSE)
  }

  call <- expr[[2L]]
  if (!is.call(call) || !identical(call[[1L]], as.name(".Deprecated"))) {
    stop("Deprecated shim `", name, "` must call .Deprecated() first.",
         call. = FALSE)
  }

  args <- as.list(call)[-1L]
  arg_names <- names(args)
  if (!is.null(arg_names) && "msg" %in% arg_names &&
      !"new" %in% arg_names) {
    return(NA_character_)
  }

  new <- if (!is.null(arg_names) && "new" %in% arg_names) {
    args[[which(arg_names == "new")[[1L]]]]
  } else {
    args[[1L]]
  }
  value <- eval(new, envir = environment(fun))
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("Deprecated shim `", name, "` has an invalid successor.", call. = FALSE)
  }
  value
}
