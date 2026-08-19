api_expect_set_allowlist <- function(observed, expected, info = NULL) {
  normalise <- function(x) sort(unique(as.character(x)))
  expect_identical(normalise(observed), normalise(expected), info = info)
}

api_expect_named_allowlist <- function(observed, expected, info = NULL) {
  normalise <- function(x) {
    if (!length(x)) {
      return(stats::setNames(character(0L), character(0L)))
    }
    x
  }
  expect_identical(normalise(observed), normalise(expected), info = info)
}

api_expect_machine_readable_targets <- function(api, deprecated, message_only) {
  machine_readable <- setdiff(
    as.character(deprecated),
    as.character(message_only)
  )
  checked <- character(0L)
  for (name in machine_readable) {
    expect_identical(
      api$superseded_by[api$name == name],
      bulkiRNA:::.bulkirna_deprecation_target(name)
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, machine_readable)
}

test_that("the API registry covers the complete namespace exactly once", {
  api <- bulkirna_api(quiet = TRUE)
  exports <- getNamespaceExports("bulkiRNA")

  expect_s3_class(api, "tbl_df")
  expect_identical(
    names(api),
    c("name", "layer", "lifecycle", "frozen", "stochastic", "superseded_by",
      "removed_in")
  )
  expect_equal(nrow(api), length(exports))
  expect_equal(anyDuplicated(api$name), 0L)
  expect_setequal(api$name, exports)
  expect_identical(api$name, sort(api$name))
})

test_that("each export has one stability lifecycle and one layer", {
  api <- bulkirna_api(quiet = TRUE)

  expect_false(anyNA(api$lifecycle))
  expect_false(anyNA(api$layer))
  expect_true(all(nzchar(api$layer)))
  expect_true(all(api$lifecycle %in% c("stable", "experimental", "deprecated")))
  expect_identical(
    api$layer[api$name == "download_gatom_references"],
    "gatom"
  )
  expect_false(any(
    api$lifecycle == "deprecated" & api$layer == "deprecated"
  ))
})

test_that("experimental status covers CoReSh, coregulation and gene-id helpers", {
  api <- bulkirna_api(quiet = TRUE)
  experimental <- c(
    "coresh_chunks", "coresh_convergence", "coresh_loadings", "coresh_match",
    "coresh_search", "coresh_sets", "coresh_validate", "gs_coregulation",
    "entrez_to_gene", "filter_confounder_genes", "gene_to_entrez"
  )

  expect_setequal(api$name[api$lifecycle == "experimental"], experimental)
  expect_true(all(api$layer[grepl("^coresh_", api$name)] == "coresh"))
  expect_identical(api$layer[api$name == "gs_coregulation"], "gs")
})

test_that("the historical signature freeze remains an independent axis", {
  api <- bulkirna_api(quiet = TRUE)
  frozen <- c(
    "normalize_gsea_results", "run_gsea", "create_standard_volcano",
    "format_pathway_name", "gsea_running_sum_plot", "list_to_term2gene",
    "gsea_barplot", "gsea_dotplot", "load_reference_db",
    "custom_minimal_theme_with_grid", "gsea_dotplot_facet", "create_MD_plot",
    "empty_gsea_tibble", "ensure_dir", "run_gsea_analysis", "save_gsea_log",
    "plot_all_gsea_results", "convert_human_to_mouse", "parse_gmx",
    "parse_mitoxplorer", "filter_by_size", "build_dge",
    "list_reference_dbs", "download_gatom_references"
  )

  expect_equal(sum(api$frozen), length(frozen))
  expect_setequal(api$name[api$frozen], frozen)
  expect_true(all(api$frozen[api$lifecycle == "deprecated"]))
  api_expect_set_allowlist(
    api$name[api$frozen & api$lifecycle == "stable"],
    c("build_dge", "ensure_dir", "format_pathway_name")
  )
  api_expect_set_allowlist(
    character(0L),
    NULL,
    info = "The frozen-and-stable allowlist must be shape-stable when empty."
  )
})

test_that("deprecated metadata names a successor and removal version", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api[api$lifecycle == "deprecated", ]
  maintained <- api[api$lifecycle != "deprecated", ]

  expect_true(all(!is.na(deprecated$superseded_by)))
  expect_true(all(nzchar(deprecated$superseded_by)))
  expect_true(all(deprecated$removed_in == "1.0.0"))
  expect_true(all(is.na(maintained$superseded_by)))
  expect_true(all(is.na(maintained$removed_in)))

  no_deprecated <- deprecated[0L, , drop = FALSE]
  expect_true(all(!is.na(no_deprecated$superseded_by)))
  expect_true(all(nzchar(no_deprecated$superseded_by)))
  expect_true(all(no_deprecated$removed_in == "1.0.0"))
})

test_that("machine-readable successors come from the deprecation calls", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]
  # Message-only shims carry no machine-readable target; every other target is
  # parsed out of the call itself.
  message_only <- c(
    "filter_by_size", "parse_mitoxplorer", "convert_human_to_mouse",
    "plot_all_gsea_results", "save_gsea_log", "empty_gsea_tibble"
  )
  # One literal anchor, so the loop below cannot pass with a helper that reads
  # the wrong argument and still agrees with itself.
  if (length(deprecated)) {
    expect_true("run_gsea" %in% deprecated)
    expect_identical(
      api$superseded_by[api$name == "run_gsea"],
      "gs_ranks() and gs_test()"
    )
  }

  api_expect_machine_readable_targets(api, deprecated, message_only)
  api_expect_machine_readable_targets(
    api[0L, , drop = FALSE],
    character(0L),
    character(0L)
  )
})

test_that("message-only shims describe the complete migration path", {
  api <- bulkirna_api(quiet = TRUE)
  expected <- c(
    filter_by_size = paste0(
      "the min_size/max_size arguments on gsdb_msigdb(), gsdb_load() and ",
      "gsdb_from_file()"
    ),
    parse_mitoxplorer = paste0(
      "gsdb_load(\"mitoxplorer\"), or gsdb_from_file() for an arbitrary ",
      "file"
    ),
    convert_human_to_mouse =
      "gsdb_msigdb(species = \"Mus musculus\", db_species = \"HS\")",
    plot_all_gsea_results = paste0(
      "gs_plot_dot(), gs_plot_bar(), gs_plot_running() and gs_save()"
    ),
    save_gsea_log =
      "gs_save(); the free-text log has no replacement, by design"
  )

  observed <- api$superseded_by[match(names(expected), api$name)]
  api_expect_named_allowlist(
    stats::setNames(observed, names(expected)),
    expected
  )
  api_expect_named_allowlist(
    character(0L),
    NULL,
    info = "The message-only shim allowlist must be shape-stable when empty."
  )
})

test_that("bulkirna_api is itself a stable top-level export", {
  row <- bulkirna_api(quiet = TRUE)
  row <- row[row$name == "bulkirna_api", ]

  expect_equal(nrow(row), 1L)
  expect_identical(row$lifecycle, "stable")
  expect_identical(row$layer, "top-level")
  expect_false(row$frozen)
  expect_true(is.na(row$superseded_by))
  expect_true(is.na(row$removed_in))
})

test_that("G4 widens a stable vocabulary with an experimental verb", {
  api <- bulkirna_api(quiet = TRUE)
  stat_types <- api[api$name == "gs_stat_types", ]
  coregulation <- api[api$name == "gs_coregulation", ]

  expect_identical(stat_types$lifecycle, "stable")
  expect_false(stat_types$frozen)
  expect_identical(coregulation$lifecycle, "experimental")
  expect_identical(coregulation$layer, "gs")
  expect_false(coregulation$frozen)
})

api_expect_deprecated_first <- function(deprecated) {
  checked <- character(0L)
  for (name in deprecated) {
    fun <- getExportedValue("bulkiRNA", name)
    first <- body(fun)[[2L]]
    expect_true(
      is.call(first) && identical(first[[1L]], as.name(".Deprecated")),
      info = name
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, deprecated)
}

api_expect_deprecated_warnings <- function(deprecated) {
  checked <- character(0L)
  for (name in deprecated) {
    fun <- getExportedValue("bulkiRNA", name)
    # The test above establishes that this is the shim's first statement.
    # Install a body-truncated copy under the real shim name. This preserves
    # .Deprecated()'s call-derived `old` value without allowing delegated work
    # in this or a future shim to run.
    deprecation_call <- body(fun)[[2L]]
    warning_fun <- fun
    body(warning_fun) <- as.call(list(as.name("{"), deprecation_call))
    shim_env <- new.env(parent = environment(fun))
    assign(name, warning_fun, envir = shim_env)
    # Both shim forms name themselves, in two spellings: `.Deprecated("target")`
    # produces `'name' is deprecated` from the call, and the six message-only
    # shims write `` `name()` is deprecated `` in their own prose. Matching the
    # bare name covers both, and still fails if the truncated body loses the
    # call-derived `old` value and the warning reads `'.Deprecated' is
    # deprecated` instead.
    expect_warning(
      eval(call(name), envir = shim_env),
      name,
      fixed = TRUE,
      class = "deprecatedWarning"
    )
    checked <- c(checked, name)
  }
  expect_identical(checked, deprecated)
}

test_that("every deprecated export calls .Deprecated first", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]

  api_expect_deprecated_first(deprecated)
  api_expect_deprecated_first(character(0L))
})

test_that("every deprecated export warns before doing any work", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]

  api_expect_deprecated_warnings(deprecated)
  api_expect_deprecated_warnings(character(0L))
})

test_that("a name's layer is invariant to lifecycle selection", {
  api <- bulkirna_api(quiet = TRUE)
  for (lifecycle in c("stable", "experimental", "deprecated")) {
    selected <- bulkirna_api(lifecycle, quiet = TRUE)
    expect_identical(
      selected$layer,
      api$layer[match(selected$name, api$name)]
    )
  }
  selected <- api[0L, , drop = FALSE]
  expect_identical(
    selected$layer,
    api$layer[match(selected$name, api$name)]
  )
})

test_that("lifecycle selection follows the dependency-report interface", {
  api <- bulkirna_api(quiet = TRUE)
  stable <- bulkirna_api("stable", quiet = TRUE)
  selected <- bulkirna_api(c("experimental", "deprecated"), quiet = TRUE)

  expect_true(all(stable$lifecycle == "stable"))
  api_expect_set_allowlist(
    unique(selected$lifecycle),
    intersect(c("experimental", "deprecated"), unique(api$lifecycle))
  )
  api_expect_set_allowlist(
    character(0L),
    NULL,
    info = "Lifecycle selection must be shape-stable when a tier is empty."
  )
  expect_error(
    bulkirna_api("not-a-lifecycle", quiet = TRUE),
    "should be one of"
  )
  expect_error(
    bulkirna_api(quiet = NA),
    "`quiet` must be a single non-missing logical value"
  )
  expect_invisible(bulkirna_api("experimental", quiet = FALSE))
})
