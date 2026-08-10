#' Ensure a directory exists
#'
#' Creates each directory in `path` if it is missing, including parents.
#' Idempotent. `path` is a **directory** path, never a file path — if you have a
#' file path, pass `dirname(path)` (or use the internal [ensure_parent_dir()]).
#'
#' @param path Character vector of directory paths.
#' @return `path`, invisibly, so the call can be inlined.
#' @examples
#' d <- file.path(tempdir(), "figures")
#' ensure_dir(d)
#' dir.exists(d)
#' @export
ensure_dir <- function(path) {
  if (is.null(path) || length(path) == 0L) {
    stop("`path` must be a non-empty character vector of directory paths.", call. = FALSE)
  }
  if (!is.character(path) || anyNA(path) || any(!nzchar(path))) {
    stop("`path` must contain non-missing, non-empty directory paths.", call. = FALSE)
  }
  for (p in path) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Ensure the parent directory of a file path exists
#'
#' File-path counterpart to [ensure_dir()]. Internal; kept separate so the two
#' semantics never live under one name again.
#'
#' @param path Character vector of **file** paths.
#' @return `path`, invisibly.
#' @keywords internal
ensure_parent_dir <- function(path) {
  ensure_dir(unique(dirname(path)))
  invisible(path)
}

#' Default value for NULL
#'
#' @param x,y `x` unless it is `NULL`, in which case `y`.
#' @return `x` or `y`.
#' @name null-default
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
