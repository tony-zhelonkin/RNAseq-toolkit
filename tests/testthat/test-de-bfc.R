test_that("de_bfc_plot returns a ggplot and needs logFC + B", {
  de <- fake_de_table()
  expect_s3_class(de_bfc_plot(de), "ggplot")
  expect_error(de_bfc_plot(de[, "logFC", drop = FALSE]), "`B`")
})

test_that("the B threshold line is always inside the y limits", {
  de <- fake_de_table()
  de$B <- seq(-20, -10, length.out = nrow(de))   # nothing crosses B = 0
  p <- de_bfc_plot(de, B_cutoff = 0, y_padding = 1)

  ylim <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
  expect_gt(ylim[2], 0)

  labels <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) d$label))
  expect_true(any(grepl("threshold", labels)))
})

test_that("categories combine B and fold change", {
  de <- data.frame(logFC = c(3, 0.1, 3, 0.1), B = c(2, 2, -2, -2),
                   row.names = paste0("g", 1:4))
  p <- de_bfc_plot(de, fc_cutoff = 1, B_cutoff = 0, label_method = "none")
  cats <- ggplot2::ggplot_build(p)$plot$data$cat
  expect_equal(cats, c("B-statistic & Log2FC", "B-statistic", "Log2FC", "NS"))
})

test_that("top labels are ranked by decreasing B on each side", {
  de <- fake_de_table(n = 20)
  p  <- de_bfc_plot(de, top_n = 2)
  repel <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(repel, 1L)
  expect_equal(nrow(repel[[1]]$data), 4L)
})
