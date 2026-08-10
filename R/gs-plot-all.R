#' Render the standard plot set for a gene-set result
#'
#' The body behind the deprecated `plot_all_gsea_results()`, rebuilt on the new
#' renderers: for each database in `x` it writes an up dotplot, a down dotplot,
#' a direction-faceted dotplot, a barplot and a text log under
#' `<out_dir>/<database>/`. Every figure goes through [gs_save()], so each one
#' arrives with its source table beside it.
#'
#' Internal on purpose -- the new golden path is to call the renderer you want
#' and save it. This exists so the deprecation shim has something to delegate
#' to.
#'
#' @param x A [gs_result].
#' @param out_dir Output directory; created if missing.
#' @param name Stem prepended to every file name.
#' @param top Number of pathways per plot.
#' @param padj_cutoff FDR threshold for highlighting and for the log.
#' @param width,height Figure size in inches.
#' @param verbose Logical. Report progress with [message()].
#' @return A character vector of every written path, invisibly.
#' @keywords internal
.gs_plot_all <- function(x, out_dir, name = "gsea", top = 20,
                         padj_cutoff = 0.05, width = 8, height = 6,
                         verbose = FALSE) {
  .gs_plot_check_result(x)
  if (!is.character(out_dir) || length(out_dir) != 1L || !nzchar(out_dir)) {
    stop("`out_dir` must be a single non-empty directory path.",
         call. = FALSE)
  }
  ensure_dir(out_dir)
  written <- character(0)

  for (db in unique(x[["database"]])) {
    part <- x[x[["database"]] == db, , drop = FALSE]
    if (nrow(part) == 0L) next
    db_dir <- file.path(out_dir, db)
    ensure_dir(db_dir)
    stem <- file.path(db_dir, paste0(name, "_", db))
    if (verbose) {
      message("Rendering ", db, " (", nrow(part), " pathways) into ", db_dir)
    }

    specs <- list(
      list(
        suffix = "_up_dot",
        plot = function() {
          gs_plot_dot(part, top = top, direction = "up",
                      highlight = padj_cutoff, title = paste0(db, ": up"))
        },
        height = height
      ),
      list(
        suffix = "_down_dot",
        plot = function() {
          gs_plot_dot(part, top = top, direction = "down",
                      highlight = padj_cutoff, title = paste0(db, ": down"))
        },
        height = height
      ),
      list(
        suffix = "_facet_dot",
        plot = function() {
          gs_plot_dot(part, top = top, facet = "direction",
                      highlight = padj_cutoff, title = db)
        },
        height = height * 1.4
      ),
      list(
        suffix = "_bar",
        plot = function() {
          gs_plot_bar(part, top = top, highlight = padj_cutoff, title = db)
        },
        height = height
      )
    )

    for (spec in specs) {
      written <- c(written, gs_save(
        spec$plot(), paste0(stem, spec$suffix),
        width = width, height = spec$height
      ))
    }
    written <- c(
      written,
      .gs_write_log(part, paste0(stem, "_log.txt"),
                    padj_cutoff = padj_cutoff)
    )
  }
  invisible(written)
}
