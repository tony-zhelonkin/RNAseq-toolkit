test_that("gs_plot_bar draws one bar per selected pathway", {
  res <- fake_plot_result(n = 8L)
  p <- gs_plot_bar(res, top = 4)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(layer_data_for(p)), 4L)
})

test_that("the value axis is labelled from stat_type", {
  expect_equal(gs_plot_bar(fake_plot_result())$labels$y, "NES")
  expect_equal(
    gs_plot_bar(fake_plot_result(stat_type = "t"))$labels$y, "t statistic"
  )
})

test_that("non-significant pathways are shown, not dropped", {
  res <- fake_plot_result(n = 8L)
  src <- attr(gs_plot_bar(res, top = 8), "gs_source")
  expect_equal(nrow(src), 8L)
  expect_true(any(!src$significant))
})

test_that("significant bars are outlined and others are transparent", {
  src <- attr(gs_plot_bar(fake_plot_result(n = 8L), top = 8), "gs_source")
  expect_setequal(unique(src$outline), c("black", "transparent"))
  expect_true(all(src$outline[src$significant] == "black"))
})

test_that("padj_max applies a hard filter before selection", {
  res <- fake_plot_result(n = 8L)
  src <- attr(gs_plot_bar(res, top = 8, padj_max = 0.05), "gs_source")
  expect_true(all(src$padj < 0.05))
  expect_lt(nrow(src), 8L)
})

test_that("selection by absolute statistic keeps both tails", {
  res <- fake_plot_result(n = 10L)
  src <- attr(gs_plot_bar(res, top = 4, sort_by = "stat"), "gs_source")
  expect_true(any(src$stat > 0) && any(src$stat < 0))
})

test_that("facets read database_label", {
  res <- fake_plot_result(n = 3L, databases = c("msigdb_H", "mitopathways"))
  src <- attr(gs_plot_bar(res, top = 2, facet = "database"), "gs_source")
  expect_true(all(grepl("^Label for ", as.character(src$.facet_row))))
})

test_that("an empty result gives an empty plot", {
  p <- gs_plot_bar(fake_plot_result(n = 3L)[0, ])
  expect_equal(p$labels$subtitle, "No pathways to plot")
})

test_that("bad input is rejected", {
  expect_error(gs_plot_bar(list()), "`x` must be a gs_result")
  expect_error(gs_plot_bar(fake_plot_result(), padj_max = "low"),
               "`padj_max`")
})
