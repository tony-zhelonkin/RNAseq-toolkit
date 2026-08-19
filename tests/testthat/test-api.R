test_that("the API registry covers the complete namespace exactly once", {
  api <- bulkirna_api(quiet = TRUE)
  exports <- getNamespaceExports("bulkiRNA")

  expect_s3_class(api, "tbl_df")
  expect_identical(
    names(api),
    c("name", "layer", "lifecycle", "frozen", "stochastic", "superseded_by",
      "removed_in")
  )
  expect_equal(nrow(api), 75L)
  expect_equal(anyDuplicated(api$name), 0L)
  expect_setequal(api$name, exports)
  expect_identical(api$name, sort(api$name))
})

test_that("each export has one stability lifecycle and one layer", {
  api <- bulkirna_api(quiet = TRUE)

  expect_equal(sum(api$lifecycle == "deprecated"), 21L)
  expect_equal(sum(api$lifecycle == "experimental"), 8L)
  expect_equal(sum(api$lifecycle == "stable"), 46L)
  expect_false(anyNA(api$lifecycle))
  expect_false(anyNA(api$layer))
  expect_true(all(nzchar(api$layer)))
  expect_true(all(api$lifecycle %in% c("stable", "experimental", "deprecated")))
  expect_equal(nrow(api), 75L)
  expect_equal(sum(api$layer == "gs"), 17L)
  expect_equal(sum(api$layer == "gsdb"), 6L)
  expect_equal(sum(api$layer == "de"), 6L)
  expect_equal(sum(api$layer == "gatom"), 7L)
  expect_equal(sum(api$layer == "coresh"), 5L)
  expect_equal(sum(api$layer == "gsea"), 4L)
  expect_equal(sum(api$layer == "top-level"), 30L)
  expect_identical(
    api$layer[api$name == "download_gatom_references"],
    "gatom"
  )
  expect_false(any(
    api$lifecycle == "deprecated" & api$layer == "deprecated"
  ))
})

test_that("experimental status is limited to CoReSh and gene-id helpers", {
  api <- bulkirna_api(quiet = TRUE)
  experimental <- c(
    "coresh_chunks", "coresh_convergence", "coresh_match", "coresh_search",
    "coresh_validate", "entrez_to_gene", "filter_confounder_genes",
    "gene_to_entrez"
  )

  expect_setequal(api$name[api$lifecycle == "experimental"], experimental)
  expect_true(all(api$layer[grepl("^coresh_", api$name)] == "coresh"))
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

  expect_equal(sum(api$frozen), 24L)
  expect_setequal(api$name[api$frozen], frozen)
  expect_true(all(api$frozen[api$lifecycle == "deprecated"]))
  expect_setequal(
    api$name[api$frozen & api$lifecycle == "stable"],
    c("build_dge", "ensure_dir", "format_pathway_name")
  )
})

test_that("deprecated metadata names a successor and removal version", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api[api$lifecycle == "deprecated", ]
  maintained <- api[api$lifecycle != "deprecated", ]

  expect_equal(nrow(deprecated), 21L)
  expect_true(all(!is.na(deprecated$superseded_by)))
  expect_true(all(nzchar(deprecated$superseded_by)))
  expect_true(all(deprecated$removed_in == "1.0.0"))
  expect_true(all(is.na(maintained$superseded_by)))
  expect_true(all(is.na(maintained$removed_in)))
})

test_that("machine-readable successors come from the deprecation calls", {
  api <- bulkirna_api(quiet = TRUE)
  # Six shims use `.Deprecated(msg = )` and so carry no machine-readable
  # target; the other fifteen are parsed out of the call itself.
  message_only <- c(
    "filter_by_size", "parse_mitoxplorer", "convert_human_to_mouse",
    "plot_all_gsea_results", "save_gsea_log", "empty_gsea_tibble"
  )
  machine_readable <- setdiff(api$name[api$lifecycle == "deprecated"], message_only)
  expect_length(machine_readable, 15L)

  # One literal anchor, so the loop below cannot pass with a helper that reads
  # the wrong argument and still agrees with itself.
  expect_identical(
    api$superseded_by[api$name == "run_gsea"],
    "gs_ranks() and gs_test()"
  )

  for (name in machine_readable) {
    expect_identical(
      api$superseded_by[api$name == name],
      bulkiRNA:::.bulkirna_deprecation_target(name)
    )
  }
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
  expect_identical(stats::setNames(observed, names(expected)), expected)
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

test_that("every deprecated export calls .Deprecated first", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]

  for (name in deprecated) {
    fun <- getExportedValue("bulkiRNA", name)
    first <- body(fun)[[2L]]
    expect_true(
      is.call(first) && identical(first[[1L]], as.name(".Deprecated")),
      info = name
    )
  }
})

test_that("every deprecated export warns before doing any work", {
  api <- bulkirna_api(quiet = TRUE)
  deprecated <- api$name[api$lifecycle == "deprecated"]

  for (name in deprecated) {
    fun <- getExportedValue("bulkiRNA", name)
    # The test above establishes that this is the shim's first statement.
    # Evaluate only that statement: this proves the real warning fires without
    # allowing any delegated work in this or a future shim to run.
    deprecation_call <- body(fun)[[2L]]
    expect_warning(
      eval(deprecation_call, envir = environment(fun)),
      class = "deprecatedWarning"
    )
  }
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
})

test_that("lifecycle selection follows the dependency-report interface", {
  stable <- bulkirna_api("stable", quiet = TRUE)
  selected <- bulkirna_api(c("experimental", "deprecated"), quiet = TRUE)

  expect_true(all(stable$lifecycle == "stable"))
  expect_setequal(unique(selected$lifecycle), c("experimental", "deprecated"))
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
