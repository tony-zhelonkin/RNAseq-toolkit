# Contract tests for gs_plot_running(). Plot assertions read
# ggplot_build(p)$data / $layout, never pixels.

rs_ranks <- function(n = 60L) {
  stats::setNames(seq(3, -3, length.out = n), paste0("G", seq_len(n)))
}

rs_sets <- function() {
  list(
    SET_A = paste0("G", c(1:6, 30)),
    SET_B = paste0("G", c(55:60, 20)),
    SET_C = paste0("G", seq(2, 60, by = 6))
  )
}

rs_db <- function() {
  sets <- rs_sets()
  structure(
    sets,
    pathway_names = c(SET_A = "Alpha response", SET_B = "Beta response",
                      SET_C = "Gamma response"),
    database = "testdb", species = "Homo sapiens", gene_id_type = "symbol",
    class = "gs_db"
  )
}

rs_result <- function() {
  bulkiRNA:::gs_result(
    data.frame(
      pathway_id = c("SET_A", "SET_B", "SET_C"),
      pathway_name = c("Alpha response", "Beta response", "Gamma response"),
      n_genes = c(7L, 7L, 10L),
      n_genes_tested = c(7L, 7L, 10L),
      stat = c(2.1, -1.8, 0.4),
      p_value = c(0.001, 0.01, 0.4),
      padj = c(0.01, 0.05, 0.6),
      stringsAsFactors = FALSE
    ),
    database = "testdb", contrast = "KO-WT", method = "fgsea",
    stat_type = "NES"
  )
}

test_that("returns a plain ggplot with three panels sharing one x scale", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))
  b <- ggplot2::ggplot_build(p)
  expect_equal(nrow(b$layout$layout), 3L)
  expect_equal(as.character(b$layout$layout$panel), c("es", "ticks", "stats"))
  xr <- lapply(b$layout$panel_params, function(pp) pp$x$get_limits())
  expect_equal(xr[[1]], xr[[2]])
  expect_equal(xr[[1]], xr[[3]])
})

test_that("panel heights honour panel_heights exactly", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       panel_heights = c(3, 1, 2))
  g <- ggplot2::ggplotGrob(p)
  h <- as.numeric(g$heights[grepl("null", as.character(g$heights))])
  expect_length(h, 3L)
  expect_equal(h / h[2], c(3, 1, 2), tolerance = 1e-4)
})

test_that("panel_heights is validated", {
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               panel_heights = c(1, 2)),
               "three positive finite numbers")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               panel_heights = c(1, 0, 2)),
               "three positive finite numbers")
})

test_that("colours are keyed by pathway id, not by position or label", {
  # Palette given in reverse-alphabetical order, and labels that would sort
  # the other way: the old positional bug would swap these.
  pal <- c(SET_C = "#000001", SET_B = "#000002", SET_A = "#000003")
  p <- gs_plot_running(rs_db(), ranks = rs_ranks(),
                       pathways = c("SET_A", "SET_B", "SET_C"),
                       palette = pal,
                       labels = c(SET_A = "zzz", SET_C = "aaa"))
  gd <- ggplot2::get_guide_data(p, "colour")
  expect_equal(as.character(gd$.value), c("SET_A", "SET_B", "SET_C"))
  expect_equal(gd$colour, c("#000003", "#000002", "#000001"))
  expect_equal(gd$.label, c("zzz", "Beta response", "aaa"))
})

test_that("an unnamed palette zips to the declared pathway order", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       pathways = c("SET_C", "SET_A"),
                       palette = c("#111111", "#222222"))
  gd <- ggplot2::get_guide_data(p, "colour")
  expect_equal(gd$colour[as.character(gd$.value) == "SET_C"], "#111111")
  expect_equal(gd$colour[as.character(gd$.value) == "SET_A"], "#222222")
})

test_that("a partially matching named palette warns and falls back", {
  expect_warning(
    gs_plot_running(rs_sets(), ranks = rs_ranks(),
                    pathways = c("SET_A", "SET_B"),
                    palette = c(SET_A = "#111111")),
    "do not cover every plotted pathway"
  )
})

test_that("the ES curve comes from fgsea, unaltered up to the panel window", {
  ranks <- rs_ranks()
  sets <- rs_sets()
  p <- gs_plot_running(sets, ranks = ranks, pathways = "SET_A")
  b <- ggplot2::ggplot_build(p)
  # the geom_line layer is the last one added
  drawn <- b$data[[length(b$data)]]
  ref <- fgsea::plotEnrichmentData(pathway = sets$SET_A, stats = ranks)$curve
  expect_equal(nrow(drawn), nrow(ref))
  expect_equal(drawn$x, as.numeric(ref$rank))
  # window mapping is linear, so correlation is exactly 1
  expect_equal(stats::cor(drawn$y, ref$ES), 1, tolerance = 1e-8)
})

test_that("y labels read in original units and the tick panel has none", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks())
  b <- ggplot2::ggplot_build(p)
  es_labs <- as.numeric(b$layout$panel_params[[1]]$y$get_labels())
  expect_true(all(abs(es_labs) <= 1.2))
  expect_length(b$layout$panel_params[[2]]$y$get_breaks(), 0L)
  met_labs <- as.numeric(b$layout$panel_params[[3]]$y$get_labels())
  expect_true(max(met_labs) >= 2 && min(met_labs) <= -2)
})

test_that("a gs_result selects its top pathways by abs(stat) and needs a db", {
  res <- rs_result()
  ranks <- rs_ranks()
  expect_error(gs_plot_running(res, ranks = ranks), "`db` is required")
  p <- gs_plot_running(res, ranks = ranks, db = rs_db(), top = 2)
  gd <- ggplot2::get_guide_data(p, "colour")
  expect_equal(sort(as.character(gd$.value)), c("SET_A", "SET_B"))
  expect_equal(gd$.label[as.character(gd$.value) == "SET_A"], "Alpha response")
})

test_that("integer pathways index gs_result rows, but not a bare set list", {
  res <- rs_result()
  p <- gs_plot_running(res, ranks = rs_ranks(), db = rs_db(),
                       pathways = c(3L, 1L))
  gd <- ggplot2::get_guide_data(p, "colour")
  expect_equal(as.character(gd$.value), c("SET_C", "SET_A"))
  expect_error(gs_plot_running(res, ranks = rs_ranks(), db = rs_db(),
                               pathways = 99L),
               "index outside the 3 rows")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(), pathways = 1L),
               "index the rows of a `gs_result`")
})

test_that("ranks and gene sets can arrive as attributes of x", {
  x <- rs_result()
  attr(x, "ranks") <- rs_ranks()
  attr(x, "gene_sets") <- rs_sets()
  expect_s3_class(gs_plot_running(x, top = 1), "ggplot")
})

test_that("missing or malformed ranks error clearly", {
  expect_error(gs_plot_running(rs_sets()), "`ranks` is required")
  expect_error(gs_plot_running(rs_sets(), ranks = 1:5),
               "must be a \\*named\\* numeric vector")
})

test_that("unknown pathways and empty overlaps error", {
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               pathways = "NOPE"),
               "not found in `db`")
  expect_error(gs_plot_running(list(SET_X = c("zz1", "zz2")),
                              ranks = rs_ranks()),
               "no genes in `ranks`")
})

test_that("an empty gs_result is an error, not an empty plot", {
  res <- rs_result()[0, , drop = FALSE]
  expect_error(gs_plot_running(res, ranks = rs_ranks(), db = rs_db()),
               "no rows")
})

test_that("legend_position and metric_label reach the plot", {
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(),
                       legend_position = "none",
                       metric_label = "t statistic")
  b <- ggplot2::ggplot_build(p)
  expect_equal(b$plot$theme$legend.position, "none")
  labeller <- b$plot$facet$params$labeller
  expect_equal(unname(unlist(labeller("stats"))), "t statistic")
  expect_error(gs_plot_running(rs_sets(), ranks = rs_ranks(),
                               metric_label = c("a", "b")),
               "single string")
})

test_that("long labels are wrapped, never truncated", {
  long <- paste(rep("verylongword", 6), collapse = " ")
  p <- gs_plot_running(rs_sets(), ranks = rs_ranks(), pathways = "SET_A",
                       labels = c(SET_A = long), max_name_length = 20)
  lab <- ggplot2::get_guide_data(p, "colour")$.label[[1]]
  expect_true(grepl("\n", lab))
  expect_equal(gsub("\n", " ", lab), long)
})
