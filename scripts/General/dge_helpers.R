suppressPackageStartupMessages({
  library(edgeR)
})

# Build DGEList with integer counts (rounded if needed)
build_dge <- function(count_mat, samples_df, genes_df, round_nonint = TRUE, norm_method = "TMM") {
  cm <- as.matrix(count_mat)
  storage.mode(cm) <- "numeric"
  nonint <- any(abs(cm - round(cm)) > .Machine$double.eps^0.5)
  if (nonint) {
    if (round_nonint) {
      message("[info] rounding non-integer counts for DGEList (edgeR).")
      cm <- round(cm)
    } else {
      stop("Counts must be integers for edgeR; set round_nonint=TRUE.")
    }
  }
  storage.mode(cm) <- "integer"

  # ensure sample order matches
  stopifnot(identical(colnames(cm), rownames(samples_df)))
  # ensure gene order matches
  stopifnot(nrow(cm) == nrow(genes_df))

  dge <- DGEList(counts = cm, genes = genes_df, samples = samples_df)
  dge <- calcNormFactors(dge, method = norm_method)
  dge
}