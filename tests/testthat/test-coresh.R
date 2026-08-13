fake_coresh_object <- function(gse = "GSE1", total_var = 10) {
  list(
    gseId = gse,
    gplId = "GPL1",
    E1024 = matrix(
      c(1024L, 0L, 9999L, 9999L, 0L, 2048L),
      nrow = 3L,
      byrow = TRUE
    ),
    rownames = c(10L, 10L, 20L),
    totalVar = total_var
  )
}

test_that("coresh_match computes pct_var from stored totalVar", {
  obj <- fake_coresh_object()
  out <- coresh_match(obj, c(10L, 20L, 99L), pvalues = FALSE)

  # match() selects rows 1 and 3; the duplicate Entrez row 2 is deliberately
  # not selected. The scaled query profile is therefore c(1, 2).
  expected <- sum(c(1, 2)^2) / 2 / 10 * 100
  expect_identical(names(out),
                   c("gse", "gpl", "pct_var", "p_value", "size"))
  expect_equal(out$pct_var, expected)
  expect_identical(out$p_value, NA_real_)
  expect_identical(out$size, 2L)
})

test_that("coresh_match returns the zero row for an absent query", {
  out <- coresh_match(fake_coresh_object(), 999L)

  expect_equal(out$pct_var, 0)
  expect_identical(out$p_value, NA_real_)
  expect_identical(out$size, 0L)
})

test_that("the p-value path stops, before any chunk is read", {
  expect_error(
    coresh_match(NULL, integer(), pvalues = TRUE),
    "not wired up yet",
    fixed = TRUE
  )
  # A non-existent chunk_dir: the stop must come from `pvalues`, not from
  # resolving a directory, so the argument is rejected before any I/O.
  expect_error(
    coresh_search(list(), chunk_dir = tempfile(), pvalues = TRUE),
    "not wired up yet",
    fixed = TRUE
  )
  err <- tryCatch(
    coresh_match(NULL, integer(), pvalues = TRUE),
    error = conditionMessage
  )
  expect_match(err, "fgsea::geseca()", fixed = TRUE)
  expect_match(err, "pvalues = FALSE", fixed = TRUE)
  # The message must not repeat the withdrawn claim about `sampleSize`.
  expect_false(grepl("does not honour", err, fixed = TRUE))
  expect_false(grepl("gesecaCpp", err, fixed = TRUE))
})

test_that("coresh_match validates its supported arguments", {
  obj <- fake_coresh_object()
  expect_error(coresh_match(obj, c(10, 20)), "integer vector")
  expect_error(coresh_match(obj, integer()), "non-empty")
  expect_error(coresh_match(obj, 10L, pvalues = NA), "`pvalues`")
  expect_error(coresh_match(obj, 10L, sample_size = 0), "`sample_size`")
  expect_error(coresh_match(obj, 10L, seed = Inf), "`seed`")
  expect_error(coresh_match(obj, 10L, eps = 0), "`eps`")

  bad <- obj
  bad$totalVar <- NULL
  expect_error(coresh_match(bad, 10L), "totalVar")
  bad <- obj
  bad$rownames <- as.character(bad$rownames)
  expect_error(coresh_match(bad, 10L), "integer Entrez")
})

test_that("coresh_chunks errors clearly for an empty directory", {
  empty <- withr::local_tempdir()
  expect_error(
    coresh_chunks(empty, cache = FALSE),
    "references/synapse-data-setup.md",
    fixed = TRUE
  )
  expect_error(coresh_chunks(empty, species = "rat"), "`species`")
  expect_error(coresh_chunks(empty, cache = NA), "`cache`")
  expect_error(coresh_chunks(1), "`chunk_dir`")
})

test_that("coresh_chunks maps species aliases and caches by resolved path", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  root <- withr::local_tempdir()
  hsa <- file.path(root, "coresh", "current", "preprocessed_chunks", "hsa")
  mmu <- file.path(root, "coresh", "current", "preprocessed_chunks", "mmu")
  dir.create(hsa, recursive = TRUE)
  dir.create(mmu, recursive = TRUE)
  file.create(file.path(hsa, "001_full_objects.qs2"))
  file.create(file.path(mmu, "001_full_objects.qs2"))
  withr::local_envvar(c(REFCACHE_ROOT = root))

  reads <- 0L
  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      reads <<- reads + 1L
      gse <- if (basename(dirname(path)) == "hsa") "GSE_HSA" else "GSE_MMU"
      list(fake_coresh_object(gse))
    },
    .package = "bulkiRNA"
  )

  human <- coresh_chunks(species = "human", cache = FALSE)
  expect_identical(human$gse, "GSE_HSA")
  expect_identical(attr(human, "provenance")$species, "hsa")
  expect_identical(coresh_chunks(species = "hsa")$gse, "GSE_HSA")
  expect_identical(reads, 1L)

  mouse <- coresh_chunks(species = "mouse", cache = FALSE)
  expect_identical(mouse$gse, "GSE_MMU")
  expect_identical(attr(mouse, "provenance")$species, "mmu")
  expect_identical(coresh_chunks(species = "mmu")$gse, "GSE_MMU")
  expect_identical(reads, 2L)
})

test_that("coresh_search ranks pct_var descending without optional engines", {
  skip_if_not(exists("local_mocked_bindings", asNamespace("testthat")))
  chunks <- withr::local_tempdir()
  file.create(file.path(chunks, "001_full_objects.qs2"))
  file.create(file.path(chunks, "002_full_objects.qs2"))

  testthat::local_mocked_bindings(
    .coresh_read_chunk = function(path) {
      if (grepl("001_", basename(path), fixed = TRUE)) {
        list(fake_coresh_object("GSE_HIGH", total_var = 10))
      } else {
        list(fake_coresh_object("GSE_LOW", total_var = 20))
      }
    },
    .require_pkg = function(...) {
      stop("an optional engine was requested", call. = FALSE)
    },
    .package = "bulkiRNA"
  )

  out <- coresh_search(
    list(query_a = c(10L, 20L, 99L)),
    chunk_dir = chunks,
    n_cores = 1L,
    pvalues = FALSE
  )
  expect_identical(names(out), c(
    "query_name", "gse", "gpl", "pct_var", "p_value", "size", "rank"
  ))
  expect_identical(out$gse, c("GSE_HIGH", "GSE_LOW"))
  expect_identical(out$rank, 1:2)
  expect_true(all(is.na(out$p_value)))
  expect_identical(attr(out, "provenance")$species, "hsa")
})

test_that("coresh_search names invalid queries and validates controls", {
  chunks <- withr::local_tempdir()
  expect_error(
    coresh_search(list(short_query = c(1L, 2L)), chunks, n_cores = 1L),
    "short_query"
  )
  expect_error(coresh_search(list(c(1L, 2L, 3L)), chunks), "named list")
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, n_cores = 0),
    "`n_cores`"
  )
  expect_error(
    coresh_search(list(q = c(1L, 2L, 3L)), chunks, species = "rat"),
    "`species`"
  )
})

test_that("coresh_convergence summarizes independent top-query support", {
  ranking <- tibble::tibble(
    query_name = c("q1", "q2", "q3", "q1", "q2", "q3", "q1"),
    gse = c("GSE_A", "GSE_A", "GSE_A", "GSE_B", "GSE_B", "GSE_C", "GSE_A"),
    pct_var = c(10, 8, 6, 9, 7, 12, 1),
    rank = c(1L, 2L, 12L, 2L, 4L, 1L, 3L)
  )

  out <- coresh_convergence(ranking, top_n = 10L, min_queries = 2L)
  expect_identical(out$gse, c("GSE_A", "GSE_B"))
  expect_identical(out$n_queries, c(2L, 2L))
  expect_identical(out$queries, c("q1, q2", "q1, q2"))
  expect_identical(out$best_rank, c(1L, 2L))
  # The duplicate q1/GSE_A platform row is not counted or averaged twice.
  expect_equal(out$mean_pct_var, c(9, 8))
})

test_that("coresh_convergence validates its inputs", {
  expect_error(coresh_convergence(list()), "data frame")
  expect_error(coresh_convergence(data.frame(gse = "GSE1")), "missing column")
  ranking <- tibble::tibble(
    query_name = "q", gse = "GSE1", pct_var = 1, rank = 1L
  )
  expect_error(coresh_convergence(ranking, top_n = 0), "`top_n`")
  expect_error(coresh_convergence(ranking, min_queries = NA), "`min_queries`")
})

test_that("coresh_validate reports every preflight failure without stopping", {
  empty <- withr::local_tempdir()
  visible <- NULL
  output <- capture.output(
    visible <- withVisible(coresh_validate(empty, species = "hsa"))
  )
  out <- visible$value

  expect_s3_class(out, "tbl_df")
  expect_false(visible$visible)
  expect_identical(names(out), c("check", "ok", "detail"))
  expect_true(length(output) > 0L)
  expect_true(all(paste0("package: ", c(
    "qs2", "coresh", "BiocParallel", "org.Hs.eg.db", "org.Mm.eg.db"
  )) %in% out$check))
  coresh_row <- out[out$check == "package: coresh", ]
  expect_identical(coresh_row$ok, requireNamespace("coresh", quietly = TRUE))
  expect_match(coresh_row$detail, "no R code or callable functions", fixed = TRUE)
  expect_false(out$ok[out$check == "chunk files"])
  expect_match(
    out$detail[out$check == "chunk files"],
    "references/synapse-data-setup.md",
    fixed = TRUE
  )
  expect_error(coresh_validate(empty, species = "rat"), "`species`")
  expect_error(coresh_validate(1), "`chunk_dir`")
})

test_that("a real CoReSh chunk index is readable when the refcache is mounted", {
  resolved <- tryCatch(
    .ref_path("coresh", "preprocessed_chunks", "hsa"),
    error = function(e) NULL
  )
  skip_if_not(!is.null(resolved) && dir.exists(resolved),
              "CoReSh refcache is not mounted")
  skip_if_not_installed("qs2")

  paths <- sort(list.files(
    resolved,
    pattern = "_full_objects\\.qs2$",
    full.names = TRUE
  ))
  skip_if_not(length(paths) > 0L, "CoReSh refcache has no human chunks")
  one_chunk <- withr::local_tempdir()
  linked <- suppressWarnings(file.symlink(
    paths[[1L]],
    file.path(one_chunk, basename(paths[[1L]]))
  ))
  skip_if_not(isTRUE(linked), "This filesystem cannot link a test chunk")

  # Exercise one real ~120 MB chunk, not the entire ~20 GB species tree.
  out <- coresh_chunks(chunk_dir = one_chunk, species = "human", cache = FALSE)
  expect_gt(nrow(out), 0L)
  expect_identical(names(out), c("gse", "gpl", "chunk"))
  expect_true(all(file.exists(unique(out$chunk))))
})
