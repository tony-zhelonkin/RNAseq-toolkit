# The frozen dotplot/barplot/running-sum names take a bare `gseaResult` S4
# object. Neither clusterProfiler nor DOSE is a dependency of this package
# (see the handback report), so we define a minimal stand-in class locally,
# with just the three slots `.dep_gsea_to_gs_result()` reads.
methods::setClass(
  "gseaResult",
  representation(
    result = "data.frame",
    geneList = "numeric",
    geneSets = "list"
  )
)

fake_gsea_result <- function(n = 6L) {
  nes <- seq(-2.5, 2.5, length.out = n)
  df <- data.frame(
    ID = paste0("SET_", seq_len(n)),
    Description = paste0("HALLMARK_SET_", seq_len(n)),
    setSize = 20L + seq_len(n),
    enrichmentScore = nes / 2,
    NES = nes,
    pvalue = seq(1e-6, 0.4, length.out = n),
    p.adjust = seq(1e-4, 0.6, length.out = n),
    qvalues = seq(1e-4, 0.6, length.out = n),
    rank = seq_len(n),
    core_enrichment = vapply(seq_len(n), function(i) {
      paste(paste0("G", seq_len(1 + (i %% 4))), collapse = "/")
    }, character(1)),
    stringsAsFactors = FALSE
  )
  rownames(df) <- df$ID
  genes <- paste0("G", 1:10)
  gene_list <- stats::setNames(seq(3, -3, length.out = length(genes)), genes)
  gene_sets <- stats::setNames(
    lapply(seq_len(n), function(i) paste0("G", seq_len(1 + (i %% 4)))),
    df$ID
  )
  methods::new("gseaResult", result = df, geneList = gene_list,
              geneSets = gene_sets)
}

test_that(".dep_gsea_to_gs_result builds a valid gs_result", {
  g <- fake_gsea_result()
  x <- .dep_gsea_to_gs_result(g)
  expect_s3_class(x, "gs_result")
  expect_equal(nrow(x), 6L)
  expect_equal(unique(x$stat_type), "NES")
  expect_equal(unique(x$database), "gsea")
  expect_identical(attr(x, "ranks"), g@geneList)
  expect_identical(attr(x, "gene_sets"), g@geneSets)
})

test_that("gsea_dotplot warns and forwards to gs_plot_dot", {
  g <- fake_gsea_result()
  expect_warning(p <- gsea_dotplot(g, showCategory = 4), "deprecated")
  expect_s3_class(p, "ggplot")
})

test_that("gsea_dotplot filterBy maps direction correctly", {
  g <- fake_gsea_result(n = 8L)
  p_pos <- suppressWarnings(gsea_dotplot(g, filterBy = "NES_positive",
                                         showCategory = 20))
  expect_true(all(attr(p_pos, "gs_source")$stat > 0))

  p_neg <- suppressWarnings(gsea_dotplot(g, filterBy = "NES_negative",
                                         showCategory = 20))
  expect_true(all(attr(p_neg, "gs_source")$stat < 0))
})

test_that("gsea_dotplot_facet warns and forwards with facet = direction", {
  g <- fake_gsea_result(n = 8L)
  expect_warning(
    p <- gsea_dotplot_facet(g, showCategory = 3),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  src <- attr(p, "gs_source")
  expect_setequal(as.character(unique(src$.facet_row)), c("Up", "Down"))
})

test_that("gsea_barplot warns, hard-filters by padj_cutoff, forwards", {
  g <- fake_gsea_result(n = 8L)
  expect_warning(
    p <- gsea_barplot(g, padj_cutoff = 0.6, top_n = 20),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  expect_true(all(attr(p, "gs_source")$padj < 0.6))
})

test_that("gsea_running_sum_plot warns and forwards to gs_plot_running", {
  skip_if_not_installed("enrichplot")
  skip_if_not_installed("patchwork")
  g <- fake_gsea_result()
  expect_warning(
    p <- gsea_running_sum_plot(g, gene_set_ids = c("SET_1", "SET_2")),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
})

test_that("custom_minimal_theme_with_grid warns and forwards to theme_bulki", {
  expect_warning(th <- custom_minimal_theme_with_grid(), "deprecated")
  expect_s3_class(th, "theme")
  # The 14pt floor in theme_bulki() means the old default (12) does not
  # round-trip -- this is the documented divergence, not a bug.
  expect_equal(th$text$size, 14)

  th2 <- suppressWarnings(custom_minimal_theme_with_grid(base_size = 20))
  expect_equal(th2$text$size, 20)
})

test_that("plot_all_gsea_results warns, writes files, forwards to .gs_plot_all", {
  g1 <- fake_gsea_result(n = 6L)
  g2 <- fake_gsea_result(n = 6L)
  out <- file.path(tempdir(), "pagr")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_warning(
    written <- plot_all_gsea_results(
      list(dbA = g1, dbB = g2), analysis_name = "an1", out_root = out,
      n_pathways = 5
    ),
    "deprecated"
  )
  expect_true(length(written) > 0L)
  expect_true(all(file.exists(written)))
})

test_that("plot_all_gsea_results returns empty invisibly for an empty list", {
  expect_warning(
    written <- plot_all_gsea_results(list(), "an1", tempdir()),
    "deprecated"
  )
  expect_equal(written, character(0))
})

test_that("save_gsea_log warns, writes a log, forwards to .gs_write_log", {
  g <- fake_gsea_result(n = 6L)
  d <- file.path(tempdir(), "sgl")
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  expect_warning(
    path <- save_gsea_log(g, "log.txt", padj_cutoff = 0.6, dir = d),
    "deprecated"
  )
  expect_true(file.exists(path))
  expect_true(any(grepl("Gene-set results log", readLines(path))))
})

test_that("create_standard_volcano warns and forwards to de_volcano", {
  de <- fake_de_table()
  expect_warning(p <- create_standard_volcano(de), "deprecated")
  expect_s3_class(p, "ggplot")

  p_ref <- de_volcano(de)
  expect_equal(p$labels, p_ref$labels)
})

test_that("create_MD_plot warns and forwards to de_md_plot", {
  de <- fake_de_table()
  fit <- list(Amean = de$AveExpr)
  expect_warning(
    p <- create_MD_plot(fit, coef = "KO_vs_WT", de_results = de),
    "deprecated"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "MD plot: KO_vs_WT")
})
