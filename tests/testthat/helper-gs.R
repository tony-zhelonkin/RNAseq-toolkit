# Small builders shared by the unit tests. Not a fixture — the committed
# fixture in tests/fixtures/ belongs to the golden harness and is off limits.

fake_gs_df <- function(n = 3L, stat = NULL) {
  data.frame(
    pathway_id = paste0("SET_", seq_len(n)),
    pathway_name = paste("Set", seq_len(n)),
    n_genes = seq_len(n) * 10L,
    n_genes_tested = seq_len(n) * 8L,
    stat = stat %||% seq(-1, 1, length.out = n),
    p_value = seq(0.001, 0.2, length.out = n),
    padj = seq(0.01, 0.4, length.out = n),
    stringsAsFactors = FALSE
  )
}

fake_gs_result <- function(n = 3L, contrast = "KO-WT", ...) {
  bulkiRNA:::gs_result(
    fake_gs_df(n, ...),
    database = "testdb", contrast = contrast,
    method = "fgsea", stat_type = "NES"
  )
}

fake_gs_matrix <- function(n_path = 3L, n_samp = 4L) {
  m <- matrix(
    seq_len(n_path * n_samp) / 10,
    nrow = n_path,
    dimnames = list(paste0("SET_", seq_len(n_path)),
                    paste0("s", seq_len(n_samp)))
  )
  bulkiRNA:::gs_matrix(
    m,
    database = "testdb", method = "gsva",
    sample_data = data.frame(
      group = rep(c("WT", "KO"), length.out = n_samp),
      row.names = colnames(m)
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# `tests/fixtures/` is in .Rbuildignore, so it ships with the source tree but not
# with the built package. Tests that need it run under devtools::test() and the
# golden harness, and skip under R CMD check rather than erroring there.
coresh_micro_fixture <- function() {
  path <- testthat::test_path("..", "fixtures", "coresh-chunk-micro.rds")
  testthat::skip_if_not(file.exists(path),
                        "coresh-chunk-micro.rds is not in the built package")
  readRDS(path)
}
