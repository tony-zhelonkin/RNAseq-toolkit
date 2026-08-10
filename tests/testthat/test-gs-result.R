test_that("gs_result fills scalar metadata and derives direction", {
  res <- fake_gs_result(3L)
  expect_s3_class(res, "gs_result")
  expect_true(all(c("database", "contrast", "method", "stat_type") %in% names(res)))
  expect_identical(unique(res$database), "testdb")
  expect_identical(res$direction, c("down", "ns", "up"))
})

test_that("core columns are present, typed, and canonically ordered", {
  res <- fake_gs_result(2L)
  core <- bulkiRNA:::.gs_core_cols
  expect_identical(names(res)[seq_along(core)], core)
  expect_type(res$pathway_id, "character")
  expect_type(res$n_genes, "integer")
  expect_type(res$stat, "double")
})

test_that("a missing core column is an error naming the column", {
  df <- fake_gs_df()
  df$padj <- NULL
  expect_error(
    bulkiRNA:::gs_result(df, database = "d", contrast = "c",
                         method = "fgsea", stat_type = "NES"),
    "padj"
  )
})

test_that("stat_type is a closed vocabulary", {
  expect_error(
    bulkiRNA:::gs_result(fake_gs_df(), database = "d", contrast = "c",
                         method = "fgsea", stat_type = "odds_ratio"),
    "Unknown stat_type"
  )
  expect_true(all(
    c("NES", "t", "log2_fold_enrichment") %in% names(gs_stat_types())
  ))
})

test_that("validate_gs_result rejects a bad direction", {
  res <- fake_gs_result()
  res$direction <- "UP"
  expect_error(bulkiRNA:::validate_gs_result(res), "direction")
})

test_that("optional columns survive and stay after the core block", {
  df <- fake_gs_df(2L)
  df$es <- c(0.4, -0.4)
  res <- bulkiRNA:::gs_result(df, database = "d", contrast = "c",
                              method = "fgsea", stat_type = "NES")
  res$leading_edge <- list(c("A", "B"), "C")
  res <- bulkiRNA:::validate_gs_result(res)
  expect_true("es" %in% names(res))
  expect_gt(match("es", names(res)), length(bulkiRNA:::.gs_core_cols))
})

test_that("rbind pools contrasts and keeps the class", {
  pooled <- rbind(fake_gs_result(2L, "A-B"), fake_gs_result(2L, "C-D"))
  expect_s3_class(pooled, "gs_result")
  expect_identical(nrow(pooled), 4L)
  expect_setequal(unique(pooled$contrast), c("A-B", "C-D"))
})

test_that("rbind fills optional columns absent from one input", {
  a <- fake_gs_result(1L, "A-B")
  b <- fake_gs_result(1L, "C-D")
  b$es <- 0.5
  pooled <- rbind(a, b)
  expect_true("es" %in% names(pooled))
  expect_true(is.na(pooled$es[1]))
})

test_that("dplyr verbs keep the class, and drop it when a core column goes", {
  res <- fake_gs_result(3L)
  expect_s3_class(dplyr::filter(res, stat > 0), "gs_result")
  expect_s3_class(dplyr::arrange(res, padj), "gs_result")
  thin <- dplyr::select(res, "pathway_id", "stat")
  expect_false(inherits(thin, "gs_result"))
  expect_s3_class(thin, "tbl_df")
})

test_that("as_tibble drops the subclass without losing data", {
  res <- fake_gs_result(3L)
  tb <- tibble::as_tibble(res)
  expect_false(inherits(tb, "gs_result"))
  expect_identical(tb$stat, res$stat)
})

test_that("summary counts pathways and significant hits per database x contrast", {
  s <- summary(rbind(fake_gs_result(3L, "A-B"), fake_gs_result(3L, "C-D")),
               padj_cutoff = 0.05)
  expect_identical(nrow(s), 2L)
  expect_true(all(s$n_pathways == 3L))
  expect_true(all(s$n_sig >= s$n_up))
})

test_that("gs_stat_label reads the axis label off stat_type", {
  expect_identical(bulkiRNA:::gs_stat_label(fake_gs_result()), "NES")
  mixed <- fake_gs_result(2L)
  mixed$stat_type <- c("NES", "t")
  expect_identical(bulkiRNA:::gs_stat_label(mixed), "statistic")
})

test_that("print shows the gs_result header", {
  expect_output(print(fake_gs_result()), "gs_result")
})
