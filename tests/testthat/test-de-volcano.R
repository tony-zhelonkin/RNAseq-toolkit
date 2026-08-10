test_that("de_volcano returns a ggplot and validates its input", {
  de <- fake_de_table()
  expect_s3_class(de_volcano(de), "ggplot")
  expect_error(de_volcano(de[, c("logFC", "B")]), "missing column")
})

test_that("the FDR dashed line sits at the raw p of the FDR boundary", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1)

  boundary <- max(de$P.Value[de$adj.P.Val <= 0.05])
  expect_equal(hline_yintercepts(p), -log10(boundary))
  # ... and NOT at -log10(p_cutoff), which is the bug this guards against.
  expect_false(isTRUE(all.equal(hline_yintercepts(p), -log10(0.05))))
})

test_that("decision_by = 'p' puts the line at -log10(p_cutoff)", {
  de <- fake_de_table()
  p  <- de_volcano(de, decision_by = "p", p_cutoff = 0.05, fc_cutoff = 1)
  expect_equal(hline_yintercepts(p), -log10(0.05))
})

test_that("fixed_p_boundary overrides the derived boundary", {
  de <- fake_de_table()
  p  <- de_volcano(de, fixed_p_boundary = 1e-3)
  expect_equal(hline_yintercepts(p), 3)
})

test_that("no significant genes means no line and an annotation instead", {
  de <- fake_de_table()
  de$adj.P.Val <- 0.9
  p <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05)

  expect_length(hline_yintercepts(p), 0L)
  labels <- unlist(lapply(ggplot2::ggplot_build(p)$data, function(d) d$label))
  expect_true(any(grepl("No genes pass FDR", labels)))
  expect_match(p$labels$caption, "No genes pass FDR")
})

test_that("the caption reports the realised p, not the FDR cutoff", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, decision_by = "fdr", p_cutoff = 0.05, fc_cutoff = 1)
  boundary <- signif(max(de$P.Value[de$adj.P.Val <= 0.05]), 2)
  expect_match(p$labels$caption, sprintf("p ≤ %.2g", boundary), fixed = TRUE)
})

test_that("annotate_counts appends counts to the populated legend entry", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, fc_cutoff = 1, annotate_counts = TRUE)
  labs <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")$labels
  expect_true(any(grepl("↑", labs)))
})

test_that(".volcano_sig_counts walks the priority order", {
  df <- data.frame(logFC = c(2, -2, 0.1), sig_dec = c(TRUE, TRUE, FALSE),
                   sig_fc = c(TRUE, TRUE, FALSE))
  expect_equal(.volcano_sig_counts(df)$category, "both")
  expect_equal(.volcano_sig_counts(df)$n_up, 1L)

  df$sig_fc <- FALSE
  expect_equal(.volcano_sig_counts(df)$category, "fdr")

  df$sig_dec <- FALSE
  none <- .volcano_sig_counts(df)
  expect_equal(none$category, "none")
  expect_true(is.na(none$legend_key))
})

test_that("vertical orientation swaps the axes and the boundary line", {
  de <- fake_de_table(n = 20, n_sig = 5)
  p  <- de_volcano(de, fc_cutoff = 1, orientation = "vertical")
  boundary <- max(de$P.Value[de$adj.P.Val <= 0.05])

  expect_equal(vline_xintercepts(p), -log10(boundary))
  expect_setequal(hline_yintercepts(p), c(-1, 1))
  expect_equal(p$labels$y, "log2(FC)")
})

test_that("highlight_gene gets its own always-on label layer", {
  de <- fake_de_table()
  p  <- de_volcano(de, highlight_gene = "Gene12")
  repel <- Filter(function(l) inherits(l$geom, "GeomTextRepel"), p$layers)
  expect_length(repel, 2L)
  bold <- Filter(function(l) identical(l$aes_params$fontface, "bold"), repel)
  expect_length(bold, 1L)
  # the highlight layer carries only the highlighted gene, and never drops it
  expect_equal(rownames(bold[[1]]$data), "Gene12")
  expect_identical(bold[[1]]$geom_params$max.overlaps, Inf)
})

test_that("de_volcano_grid combines panels on shared limits", {
  skip_if_not_installed("patchwork")
  de <- fake_de_table()
  p  <- de_volcano(de, fc_cutoff = 1, orientation = "vertical")
  g  <- de_volcano_grid(list(A = p, B = p))
  expect_s3_class(g, "patchwork")
  expect_error(de_volcano_grid(list()), "non-empty list")
})
