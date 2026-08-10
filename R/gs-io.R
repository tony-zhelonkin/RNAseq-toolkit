#' Write a `gs_result` to the standard table layout
#'
#' `gs_write()` and [gs_read()] are the only functions in the compute layer
#' that touch disk. The layout under `dir` is the project convention:
#'
#' ```
#' <dir>/
#'   by_contrast/<contrast>/<name>_<database>.tsv
#'   _overview/<name>_all.tsv
#'   _overview/<name>_summary.tsv
#' ```
#'
#' `dir` is normally `03_results/<stage>/tables`. List columns
#' (`leading_edge`) are collapsed to `"/"`-separated strings on the way out and
#' split again by [gs_read()], so the round trip is lossless.
#'
#' @param x A [gs_result].
#' @param dir Directory to write into; created if missing.
#' @param name Character stem for the file names. Default `"gsea"`.
#' @param overview Logical. Also write the pooled table and the
#'   [summary.gs_result()] table under `_overview/`.
#' @param by_contrast Logical. Also write one file per contrast x database.
#' @return The directory, invisibly, with the written paths in the `"files"`
#'   attribute.
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' d <- file.path(tempdir(), "tables")
#' gs_write(res, d)
#' gs_read(d)
#' @export
gs_write <- function(x, dir, name = "gsea", overview = TRUE,
                     by_contrast = TRUE) {
  .gs_check_result(x)
  if (!is.character(dir) || length(dir) != 1L || !nzchar(dir)) {
    stop("`dir` must be a single non-empty directory path.", call. = FALSE)
  }
  ensure_dir(dir)
  files <- character(0)

  if (by_contrast) {
    parts <- gs_split(x, by = c("contrast", "database"))
    for (nm in names(parts)) {
      part <- parts[[nm]]
      if (!nrow(part)) next
      cdir <- file.path(dir, "by_contrast", .gs_slug(part$contrast[1L]))
      ensure_dir(cdir)
      f <- file.path(
        cdir, paste0(name, "_", .gs_slug(part$database[1L]), ".tsv")
      )
      .gs_write_tsv(part, f)
      files <- c(files, f)
    }
  }

  if (overview) {
    odir <- file.path(dir, "_overview")
    ensure_dir(odir)
    f_all <- file.path(odir, paste0(name, "_all.tsv"))
    .gs_write_tsv(x, f_all)
    f_sum <- file.path(odir, paste0(name, "_summary.tsv"))
    utils::write.table(
      summary(x), f_sum,
      sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    files <- c(files, f_all, f_sum)
  }

  structure(dir, files = files)
}

#' Read a `gs_result` back from the standard table layout
#'
#' Reads `_overview/<name>_all.tsv` when it exists, otherwise every
#' `by_contrast/*/<name>_*.tsv` and row-binds them. Also accepts a path to a
#' single `.tsv` file.
#'
#' @param dir Directory written by [gs_write()], or a single `.tsv` path.
#' @param name Character stem used when writing.
#' @return A [gs_result].
#' @examples
#' db <- structure(
#'   list(SET_A = c("A", "B", "C", "D"), SET_B = c("C", "D", "E", "F")),
#'   pathway_names = c(SET_A = "Set A", SET_B = "Set B"),
#'   database = "demo", species = "Homo sapiens", gene_id_type = "symbol",
#'   class = "gs_db"
#' )
#' res <- gs_test(stats::setNames(c(3, 2, 1, -1, -2, -3), LETTERS[1:6]),
#'                db, min_size = 1, max_size = 10)
#' d <- file.path(tempdir(), "tables2")
#' gs_write(res, d)
#' gs_read(d)
#' @export
gs_read <- function(dir, name = "gsea") {
  if (!is.character(dir) || length(dir) != 1L) {
    stop("`dir` must be a single path.", call. = FALSE)
  }
  paths <- if (grepl("\\.tsv$", dir)) {
    dir
  } else {
    all_f <- file.path(dir, "_overview", paste0(name, "_all.tsv"))
    if (file.exists(all_f)) {
      all_f
    } else {
      list.files(
        file.path(dir, "by_contrast"),
        pattern = paste0("^", name, "_.*\\.tsv$"),
        recursive = TRUE, full.names = TRUE
      )
    }
  }
  missing <- paths[!file.exists(paths)]
  if (!length(paths) || length(missing)) {
    stop("No gs_result tables found at ", sQuote(dir),
         ". Expected `_overview/", name, "_all.tsv` or ",
         "`by_contrast/*/", name, "_*.tsv`.", call. = FALSE)
  }
  parts <- lapply(paths, .gs_read_tsv)
  gs_result(dplyr::bind_rows(parts))
}

# ---- internals --------------------------------------------------------------

#' Write one `gs_result` as a TSV, flattening list columns
#'
#' @param x A [gs_result].
#' @param path File path.
#' @return `path`, invisibly.
#' @keywords internal
.gs_write_tsv <- function(x, path) {
  df <- as.data.frame(.as_plain_tibble(x))
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(
        df[[nm]], function(z) paste(as.character(unlist(z)), collapse = "/"),
        character(1L)
      )
    }
  }
  ensure_parent_dir(path)
  utils::write.table(
    df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  invisible(path)
}

#' Read one TSV back, restoring the `leading_edge` list column
#'
#' @param path File path.
#' @return A data frame.
#' @keywords internal
.gs_read_tsv <- function(path) {
  df <- utils::read.delim(
    path, sep = "\t", header = TRUE, quote = "",
    stringsAsFactors = FALSE, na.strings = c("", "NA"),
    check.names = FALSE
  )
  if (!is.null(df[["leading_edge"]])) {
    df[["leading_edge"]] <- lapply(df[["leading_edge"]], function(s) {
      if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, "/")[[1L]]
    })
  }
  df
}

#' Make a path-safe slug from a label
#'
#' @param x Character vector.
#' @return Character vector with path-hostile characters replaced by `_`.
#' @keywords internal
.gs_slug <- function(x) {
  x <- gsub("[^A-Za-z0-9._+-]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}
