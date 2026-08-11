#' Load an MSigDB collection
#'
#' Thin provider over [msigdbr::msigdbr()]. MSigDB is downloaded at first use
#' and cached under `$HOME/.cache/R/msigdbr`, so a cold cache needs network
#' access.
#'
#' @param species Character(1) target species, e.g. `"Mus musculus"` or
#'   `"Homo sapiens"`. Human sets are mapped to orthologs when `db_species`
#'   is `"HS"`.
#' @param collection Character(1) MSigDB collection, e.g. `"H"` for the
#'   hallmarks or `"C2"`. See [msigdbr::msigdbr_collections()].
#' @param subcollection Character(1) sub-collection such as `"CP:REACTOME"`,
#'   or `NULL` for the whole collection.
#' @param db_species Character(1) source database, `"HS"` (human MSigDB, the
#'   default, with ortholog mapping) or `"MM"` (mouse-native sets).
#' @param min_size,max_size Integer(1) or `NULL`; drop sets outside these
#'   bounds.
#' @param verbose Logical(1); message what was loaded.
#' @return A [gs_db()] whose `database` key is derived from the collection --
#'   `"msigdb_H"`, `"msigdb_C2_CP_KEGG"` -- with `"MSigDB H"` kept as the
#'   `database_label` attribute for display.
#' @examplesIf requireNamespace("msigdbr", quietly = TRUE) && interactive()
#' db <- gsdb_msigdb("Mus musculus", collection = "H")
#' summary(db)[1:3, ]
#' @export
gsdb_msigdb <- function(species = "Mus musculus",
                        collection = "H",
                        subcollection = NULL,
                        db_species = c("HS", "MM"),
                        min_size = NULL,
                        max_size = NULL,
                        verbose = FALSE) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("`gsdb_msigdb()` requires the msigdbr package. Install it with ",
         "install.packages(\"msigdbr\").", call. = FALSE)
  }
  species <- .gsdb_species_label(species)
  if (!is.character(collection) || length(collection) != 1L ||
        is.na(collection) || !nzchar(collection)) {
    stop("`collection` must be a single MSigDB collection name such as ",
         "\"H\" or \"C2\"; see `msigdbr::msigdbr_collections()`.",
         call. = FALSE)
  }
  db_species <- match.arg(db_species)

  args <- list(db_species = db_species, species = species,
               collection = collection)
  if (!is.null(subcollection)) {
    args$subcollection <- subcollection
  }
  tbl <- do.call(msigdbr::msigdbr, args)

  if (!nrow(tbl)) {
    stop("MSigDB returned no gene sets for `collection = \"", collection,
         "\"`", if (is.null(subcollection)) "" else
           paste0(", `subcollection = \"", subcollection, "\"`"),
         " and `species = \"", species, "\"`. Check the names against ",
         "`msigdbr::msigdbr_collections()`.", call. = FALSE)
  }

  sets <- split(as.character(tbl$gene_symbol), as.character(tbl$gs_name))

  # `pathway_names` is the axis/legend label (see `gs_result`), so it is the
  # formatted set name -- "Myc Targets V1" -- not MSigDB's `gs_description`,
  # which is a full sentence ("A subgroup of genes regulated by MYC - version 1
  # (v1).") and unusable as a label. The old toolkit formatted the id at every
  # call site instead; doing it once here keeps every renderer consistent.
  # The description is genuinely useful, so it is kept alongside rather than
  # dropped.
  labels <- stats::setNames(format_pathway_name(names(sets)), names(sets))
  descriptions <- NULL
  if ("gs_description" %in% names(tbl)) {
    uniq <- tbl[!duplicated(tbl$gs_name), c("gs_name", "gs_description")]
    descriptions <- stats::setNames(as.character(uniq$gs_description),
                                    as.character(uniq$gs_name))
  }

  label <- paste(c("MSigDB", collection, subcollection), collapse = " ")
  key <- paste(c("msigdb", collection,
                 if (is.null(subcollection)) NULL
                 else strsplit(subcollection, "[:.]")[[1]]),
               collapse = "_")
  db <- gs_db(sets, database = key, species = species,
              pathway_names = labels, database_label = label)
  db <- .gs_filter_size(db, min_size, max_size, verbose = verbose)
  # After the filter, not before: `[.gs_db` subsets the attributes it knows
  # about and would leave this one describing dropped sets.
  if (!is.null(descriptions)) {
    attr(db, "pathway_descriptions") <- descriptions[names(db)]
  }
  if (verbose) {
    message(sprintf("Loaded %s (%s): %d sets.", label, species, length(db)))
  }
  db
}
