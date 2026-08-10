test_that("a named numeric vector dispatches to fgsea and returns a gs_result", {
  res <- gs_test(fake_ranks(), fake_gs_db(), contrast = "KO-WT",
                 min_size = 5, max_size = 50)
  expect_s3_class(res, "gs_result")
  expect_true(all(res$stat_type == "NES"))
  expect_true(all(res$method == "fgsea"))
  expect_true(all(res$database == "testdb"))
  expect_true(all(res$contrast == "KO-WT"))
  expect_true("leading_edge" %in% names(res))
  expect_true(is.list(res$leading_edge))
  # names come from the gs_db pathway_names attribute
  expect_true(all(grepl("^Pretty ", res$pathway_name)))
})

test_that("fgsea recovers the planted direction", {
  res <- gs_test(fake_ranks(), fake_gs_db(), min_size = 5, max_size = 50)
  expect_equal(res$direction[res$pathway_id == "SET_UP"], "up")
  expect_equal(res$direction[res$pathway_id == "SET_DOWN"], "down")
})

test_that("a character vector dispatches to ORA", {
  db <- fake_gs_db()
  res <- gs_test(paste0("G", 1:20), db,
                 universe = paste0("G", 1:100),
                 min_size = 5, max_size = 50)
  expect_s3_class(res, "gs_result")
  expect_true(all(res$method == "ora"))
  expect_true(all(res$stat_type == "log2_fold_enrichment"))
  expect_true(all(c("fold_enrichment", "overlap") %in% names(res)))
  up <- res[res$pathway_id == "SET_UP", ]
  expect_equal(up$overlap, 20L)
  expect_gt(up$stat, 0)
})

test_that("ORA rejects query genes outside the universe", {
  expect_error(
    gs_test("ZZZ9", fake_gs_db(), universe = paste0("G", 1:100),
            min_size = 1),
    "universe"
  )
})

test_that("multiple databases become multiple `database` values", {
  db <- list(
    alpha = fake_gs_db(database = "alpha"),
    beta = fake_gs_db(list(BETA_1 = paste0("G", 1:15)), database = "beta")
  )
  res <- gs_test(fake_ranks(), db, min_size = 5, max_size = 50)
  expect_setequal(unique(res$database), c("alpha", "beta"))
})

test_that("pooling contrasts is rbind", {
  db <- fake_gs_db()
  a <- gs_test(fake_ranks(), db, contrast = "A", min_size = 5, max_size = 50)
  b <- gs_test(rev(fake_ranks()), db, contrast = "B",
               min_size = 5, max_size = 50)
  pooled <- rbind(a, b)
  expect_s3_class(pooled, "gs_result")
  expect_equal(nrow(pooled), nrow(a) + nrow(b))
  expect_setequal(unique(pooled$contrast), c("A", "B"))
})

test_that("gs_test errors informatively on an unsupported input", {
  expect_error(gs_test(list(1), fake_gs_db()), "no method")
  expect_error(gs_test(unname(fake_ranks()), fake_gs_db()), "named")
  expect_error(gs_test(fake_ranks(), list(a = 1)), "gs_db")
  expect_error(gs_test(fake_ranks(), fake_gs_db(), method = "ora"), "fgsea")
})

test_that("an empty result is a zero-row gs_result, not NULL", {
  db <- fake_gs_db(list(TOO_BIG = paste0("G", 1:100)))
  res <- gs_test(fake_ranks(), db, min_size = 5, max_size = 10)
  expect_s3_class(res, "gs_result")
  expect_equal(nrow(res), 0L)
})

test_that("gs_test on a gs_matrix runs limma and reports t statistics", {
  skip_if_not_installed("limma")
  m <- matrix(
    c(rnorm(4, 2), rnorm(4, -2), rnorm(8, 0)), nrow = 2, byrow = TRUE,
    dimnames = list(c("SET_UP", "SET_MID"), paste0("s", 1:8))
  )
  gm <- bulkiRNA:::gs_matrix(
    m, database = "testdb", method = "gsva",
    sample_data = fake_sample_data()
  )
  res <- gs_test(gm, design = ~ 0 + group, contrast = "groupKO-groupWT")
  expect_s3_class(res, "gs_result")
  expect_true(all(res$stat_type == "t"))
  expect_true(all(res$method == "gsva"))
  expect_equal(res$contrast[1], "groupKO-groupWT")
  expect_true(all(is.na(res$n_genes)))
})

test_that("gs_test on a gs_matrix needs a design", {
  gm <- bulkiRNA:::gs_matrix(
    matrix(1:8, nrow = 2, dimnames = list(c("a", "b"), paste0("s", 1:4))),
    database = "d", method = "gsva"
  )
  expect_error(gs_test(gm), "design")
})
