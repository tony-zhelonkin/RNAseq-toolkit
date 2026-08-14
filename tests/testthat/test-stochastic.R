test_that("every exported seed formal is declared stochastic", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  exports <- getNamespaceExports("bulkiRNA")
  has_seed <- vapply(exports, function(name) {
    "seed" %in% names(formals(getExportedValue("bulkiRNA", name)))
  }, logical(1L))

  expect_s3_class(registry, "tbl_df")
  expect_identical(
    names(registry),
    c("name", "seed_arg", "seed_default", "source_of_randomness", "note")
  )
  expect_setequal(exports[has_seed], registry$name[!is.na(registry$seed_arg)])
})

test_that("declared seed arguments and defaults match their formals", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  declared <- registry[!is.na(registry$seed_arg), ]

  for (i in seq_len(nrow(declared))) {
    fun <- getExportedValue("bulkiRNA", declared$name[[i]])
    formal_names <- names(formals(fun))
    seed_arg <- declared$seed_arg[[i]]
    seed_default <- paste(deparse(formals(fun)[[seed_arg]]), collapse = "")

    expect_true(seed_arg %in% formal_names, info = declared$name[[i]])
    expect_identical(
      declared$seed_default[[i]], seed_default,
      info = declared$name[[i]]
    )
  }
})

test_that("gs_test exposes its seed only through dots and uses it", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  through_dots <- registry[is.na(registry$seed_arg), ]

  expect_identical(through_dots$name, "gs_test")
  expect_true("..." %in% names(formals(gs_test)))
  expect_false("seed" %in% names(formals(gs_test)))
  adapter_default <- paste(
    deparse(formals(bulkiRNA:::.gs_fgsea)$seed), collapse = ""
  )
  expect_identical(through_dots$seed_default, adapter_default)

  first <- gs_test(
    fake_ranks(), fake_gs_db(), min_size = 5L, max_size = 50L,
    n_perm_simple = 1000L, seed = 101L
  )
  second <- gs_test(
    fake_ranks(), fake_gs_db(), min_size = 5L, max_size = 50L,
    n_perm_simple = 1000L, seed = 202L
  )

  expect_false(identical(first, second))
})

test_that("the API stochastic flag is derived from the registry", {
  registry <- bulkirna_stochastic(quiet = TRUE)
  api <- bulkirna_api(quiet = TRUE)

  expect_setequal(api$name[api$stochastic], registry$name)
  expect_false(anyNA(api$stochastic))
})

test_that("session provenance records RNG structure and stochastic names", {
  .clear_ref_resolutions()
  path <- tempfile(fileext = ".txt")
  write_session_provenance(path)
  lines <- readLines(path)
  registry <- bulkirna_stochastic(quiet = TRUE)

  rng_line <- grep("^RNGkind\\(\\):", lines, value = TRUE)
  expect_length(rng_line, 1L)
  expect_true(all(c("kind =", "normal.kind =", "sample.kind =") %in%
                    regmatches(
                      rng_line,
                      gregexpr("(?:normal\\.|sample\\.)?kind =", rng_line,
                               perl = TRUE)
                    )[[1L]]))
  for (name in registry$name) {
    expect_true(any(grepl(paste0("^  ", name, " = "), lines)))
  }
})
