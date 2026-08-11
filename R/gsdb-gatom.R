#' Download GATOM reference network files
#'
#' GATOM's atom-transition metabolic networks are too large to bundle
#' (~24 MB), so they are fetched on demand from the Artyomov Lab server. This
#' is a downloader, not a gene-set provider: the files it writes are consumed
#' by GATOM itself, and it does not return a [gs_db()].
#'
#' Files already present are skipped unless `overwrite = TRUE`. A file that
#' fails to download warns and is left out of the return value, so a partial
#' run is visible rather than silent.
#'
#' @param dest_dir Character(1) destination directory; created if missing.
#' @param species Character(1), `"Mus_musculus"` or `"Homo_sapiens"`.
#' @param networks Character vector of networks to fetch, e.g.
#'   `c("kegg", "combined")`.
#' @param overwrite Logical(1); re-download files that already exist.
#' @return Character vector of downloaded (or already present) file paths,
#'   invisibly.
#' @examples
#' \dontrun{
#' download_gatom_references(dest_dir = "00_data/references/gatom")
#' }
#' @export
download_gatom_references <- function(
    dest_dir = "00_data/references/gatom",
    species = "Mus_musculus",
    networks = c("kegg", "combined"),
    overwrite = FALSE
) {
  base_url <-
    "http://artyomovlab.wustl.edu/publications/supp_materials/GATOM"

  code <- switch(
    gsub(" ", "_", species, fixed = TRUE),
    Mus_musculus = list(short = "Mm", ncbi = "mmu"),
    Homo_sapiens = list(short = "Hs", ncbi = "hsa"),
    stop("`species` must be \"Mus_musculus\" or \"Homo_sapiens\"; got ",
         paste0("\"", species, "\""), ".", call. = FALSE)
  )
  if (!is.character(networks) || !length(networks) || anyNA(networks) ||
        any(!nzchar(networks))) {
    stop("`networks` must be a non-empty character vector, e.g. ",
         "c(\"kegg\", \"combined\").", call. = FALSE)
  }

  files <- unique(c(
    unlist(lapply(networks, function(net) {
      c(sprintf("network.%s.rds", net),
        sprintf("met.%s.db.rds", net),
        sprintf("gene2reaction.%s.%s.eg.tsv", net, code$ncbi))
    })),
    sprintf("org.%s.eg.gatom.anno.rds", code$short)
  ))

  ensure_dir(dest_dir)

  downloaded <- character()
  for (fname in files) {
    dest_file <- file.path(dest_dir, fname)

    # An interrupted transfer used to leave a partial file in place that every
    # later run reported as `[skip] ... (exists)`, so `gatom_refs()` then died
    # in `readRDS()` naming neither the file nor the fix. A zero-byte file is
    # treated as absent, and the transfer itself goes to `<dest>.part` and is
    # renamed only on success, so an interruption cannot produce one.
    if (file.exists(dest_file) && !isTRUE(overwrite)) {
      if (isTRUE(file.info(dest_file)$size > 0)) {
        message(sprintf("  [skip] %s (exists; use overwrite = TRUE to replace)",
                        fname))
        downloaded <- c(downloaded, dest_file)
        next
      }
      message(sprintf("  [redo] %s (present but empty -- refetching)", fname))
    }

    message(sprintf("  Downloading %s ...", fname))
    part_file <- paste0(dest_file, ".part")
    ok <- tryCatch({
      utils::download.file(file.path(base_url, fname), part_file,
                           mode = "wb", quiet = TRUE)
      if (!isTRUE(file.info(part_file)$size > 0)) {
        stop("the download produced an empty file", call. = FALSE)
      }
      if (!file.rename(part_file, dest_file)) {
        stop("could not move the download into place at ", dest_file,
             call. = FALSE)
      }
      TRUE
    }, error = function(e) {
      unlink(part_file)
      warning(sprintf("  [FAIL] %s: %s", fname, conditionMessage(e)),
              call. = FALSE)
      FALSE
    })
    if (ok) {
      downloaded <- c(downloaded, dest_file)
      message(sprintf("  [ok] %s (%.1f MB)", fname,
                      file.info(dest_file)$size / 1e6))
    }
  }

  message(sprintf("Downloaded %d / %d files to: %s",
                  length(downloaded), length(files), dest_dir))
  invisible(downloaded)
}
