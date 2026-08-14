#' Evaluate an expression with a pinned RNG seed
#'
#' `set.seed()` alone is not enough because `BiocParallel::bplapply()` changes
#' the active generator. This helper is the only place in bulkiRNA allowed to
#' seed stochastic work.
#'
#' @param seed Whole-number RNG seed, or `NULL` to evaluate `expr` without
#'   touching RNG state.
#' @param expr An expression to evaluate exactly once after seeding.
#' @return The value of `expr`.
#' @examples
#' bulkiRNA:::.with_pinned_seed(NULL, 1 + 1)
#' @keywords internal
.with_pinned_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    # RNGkind() re-randomises the seed, so restore the kind before the vector.
    # A legacy-sampler warning describes the caller's choice, not ours.
    suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}
