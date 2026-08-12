# Package-local record of reference data resolved during this session.
.ref_resolution_state <- new.env(parent = emptyenv())

#' Resolve a path in the shared reference cache
#'
#' An explicit `path` takes precedence over the shared cache. Otherwise the
#' path is resolved below `$REFCACHE_ROOT/<source>/current`, and the snapshot
#' selected by `current` is recorded for session provenance.
#'
#' @param source Character scalar naming one refcache source. It must be a
#'   single path segment.
#' @param ... Character scalar path components below the source's `current`
#'   snapshot. Ignored when `path` is supplied.
#' @param path Optional caller-supplied path. Used as given.
#' @param root Optional refcache root. `NULL` reads `REFCACHE_ROOT`.
#' @param must_exist Logical scalar; whether the returned path must exist.
#' @return A character scalar with `snapshot`, `source`, `root`, and
#'   `caller_supplied` attributes. `snapshot` and `root` are `NA_character_`
#'   for a caller-supplied path.
#' @keywords internal
.ref_path <- function(source, ..., path = NULL, root = NULL,
                      must_exist = TRUE) {
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(source) || source == "." ||
      grepl("[/\\\\]|\\.\\.", source)) {
    stop("`source` must be a single path segment -- no `/`, `\\`, `..` or ",
         "`.`; got ", sQuote(paste0(source, collapse = ", ")), ".",
         call. = FALSE)
  }
  source <- unname(source)
  if (!is.logical(must_exist) || length(must_exist) != 1L ||
      is.na(must_exist)) {
    stop("`must_exist` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  caller_supplied <- !is.null(path)
  if (caller_supplied) {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
      stop("`path` must be a non-empty character scalar or `NULL`.",
           call. = FALSE)
    }
    if (must_exist && !file.exists(path)) {
      stop("Caller-supplied reference `path` does not exist: ", path, ".",
           call. = FALSE)
    }
    resolved <- unname(as.character(path))
    snapshot <- NA_character_
    resolved_root <- NA_character_
  } else {
    components <- unname(list(...))
    valid_components <- vapply(
      components,
      function(x) is.character(x) && length(x) == 1L && !is.na(x),
      logical(1L)
    )
    if (length(valid_components) && !all(valid_components)) {
      stop("Each path component in `...` must be a character scalar.",
           call. = FALSE)
    }

    if (is.null(root)) {
      root <- Sys.getenv("REFCACHE_ROOT", unset = "")
    }
    if (!is.character(root) || length(root) != 1L || is.na(root)) {
      stop("`root` must be a character scalar or `NULL`.", call. = FALSE)
    }
    if (!nzchar(root)) {
      stop(
        "Reference data \"", source, "\" is not available. ",
        "`REFCACHE_ROOT` is unset. Bind the shared refcache into the ",
        "container and set `REFCACHE_ROOT=/refcache`, or pass `path=` to ",
        "point at a local copy.",
        call. = FALSE
      )
    }

    source_dir <- file.path(root, source)
    if (!dir.exists(source_dir)) {
      stop(
        "Reference data \"", source, "\" is not available. The source ",
        "directory is missing under `REFCACHE_ROOT`: ", source_dir, ". ",
        "Bind the shared refcache into the container and set ",
        "`REFCACHE_ROOT=/refcache`, or pass `path=` to point at a local copy.",
        call. = FALSE
      )
    }

    current <- file.path(source_dir, "current")
    if (!dir.exists(current)) {
      stop(
        "Reference data \"", source, "\" is not available. `current` is ",
        "missing or dangling under the source directory: ", current, ". ",
        "The refcache refresh may have failed. Repair the refcache snapshot, ",
        "or pass `path=` to point at a local copy.",
        call. = FALSE
      )
    }

    link_target <- Sys.readlink(current)
    if (length(link_target) != 1L || is.na(link_target) ||
        !nzchar(link_target)) {
      link_target <- normalizePath(current, winslash = "/", mustWork = TRUE)
    }
    snapshot <- basename(link_target)
    resolved_root <- unname(root)
    resolved <- if (length(components)) {
      unname(do.call(file.path, c(list(current), components)))
    } else {
      unname(current)
    }

    if (length(resolved) != 1L) {
      stop("Path components in `...` must resolve to one path.", call. = FALSE)
    }
    if (must_exist && !file.exists(resolved)) {
      stop("Reference data path for \"", source, "\" does not exist: ",
           resolved, ".", call. = FALSE)
    }
  }

  out <- structure(
    resolved,
    snapshot = snapshot,
    source = source,
    root = resolved_root,
    caller_supplied = caller_supplied
  )
  assign(
    source,
    list(source = source, snapshot = snapshot, path = resolved),
    envir = .ref_resolution_state
  )
  out
}

#' Return reference data resolved during this session
#'
#' @return A data frame with character columns `source`, `snapshot`, and
#'   `path`, containing at most one row per source.
#' @keywords internal
.ref_resolutions <- function() {
  keys <- sort(ls(envir = .ref_resolution_state, all.names = TRUE))
  if (!length(keys)) {
    return(data.frame(
      source = character(0),
      snapshot = character(0),
      path = character(0),
      stringsAsFactors = FALSE
    ))
  }

  records <- lapply(
    keys,
    get,
    envir = .ref_resolution_state,
    inherits = FALSE
  )
  data.frame(
    source = unname(vapply(records, `[[`, character(1L), "source")),
    snapshot = unname(vapply(records, `[[`, character(1L), "snapshot")),
    path = unname(vapply(records, `[[`, character(1L), "path")),
    stringsAsFactors = FALSE
  )
}

#' Clear the session reference-data record
#'
#' @return `NULL`, invisibly.
#' @keywords internal
.clear_ref_resolutions <- function() {
  keys <- ls(envir = .ref_resolution_state, all.names = TRUE)
  if (length(keys)) {
    rm(list = keys, envir = .ref_resolution_state)
  }
  invisible(NULL)
}
