#' Warn about gene identifiers that did not map
#'
#' @param unmapped Character vector of unmapped identifiers.
#' @param total Integer count of input identifiers.
#' @param label Description of the identifier type for the warning.
#' @return `NULL`, invisibly.
#' @keywords internal
.warn_unmapped_gene_ids <- function(unmapped, total, label) {
  if (!length(unmapped)) return(invisible(NULL))

  shown <- utils::head(unmapped, 20L)
  warning(
    length(unmapped), "/", total, " ", label, " failed to map: ",
    paste(shown, collapse = ", "),
    if (length(unmapped) > length(shown)) ", ..." else "",
    call. = FALSE
  )
  invisible(NULL)
}

#' Convert gene symbols to Entrez identifiers
#'
#' Maps human or mouse gene symbols through the corresponding Bioconductor
#' organism annotation package. Unmapped symbols produce a warning that names
#' at most the first 20 failures and are dropped from the result. Mapped values
#' retain their input order, including repeated symbols.
#'
#' @param symbols Character vector of gene symbols.
#' @param species Character(1): `"human"`, `"Homo sapiens"`, or `"hsa"`;
#'   alternatively `"mouse"`, `"Mus musculus"`, or `"mmu"`.
#' @param multi_vals How ambiguous mappings are handled. One of `"first"`,
#'   `"filter"`, or `"asNA"`, passed to [AnnotationDbi::mapIds()]. Strategies
#'   that return list columns are not supported because this function always
#'   returns an integer vector.
#' @return An unnamed integer vector of mapped Entrez identifiers. Unmapped
#'   inputs are omitted and the relative order of mapped inputs is preserved.
#' @examples
#' if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   gene_to_entrez(c("TP53", "EGFR"))
#' }
#' @export
gene_to_entrez <- function(symbols, species = "human", multi_vals = "first") {
  sp <- .species(species)
  if (!is.character(symbols) || anyNA(symbols) || any(!nzchar(symbols))) {
    stop("`symbols` must be a character vector of non-missing, non-empty ",
         "gene symbols.", call. = FALSE)
  }
  valid_multi_vals <- c("first", "filter", "asNA")
  if (!is.character(multi_vals) || length(multi_vals) != 1L ||
        is.na(multi_vals) || !multi_vals %in% valid_multi_vals) {
    stop("`multi_vals` must be one of \"first\", \"filter\", or \"asNA\" ",
         "so `gene_to_entrez()` can return an integer vector.", call. = FALSE)
  }
  if (!length(symbols)) return(integer(0L))

  .require_pkg(
    "AnnotationDbi", "`gene_to_entrez()`",
    'BiocManager::install("AnnotationDbi")'
  )
  .require_pkg(
    sp$orgdb,
    sprintf("`gene_to_entrez(species = \"%s\")`", sp$common),
    sprintf('BiocManager::install("%s")', sp$orgdb)
  )
  db <- getExportedValue(sp$orgdb, sp$orgdb)

  keys <- unique(symbols)
  mapped_keys <- suppressMessages(AnnotationDbi::mapIds(
    db,
    keys = keys,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = multi_vals
  ))
  mapped <- as.character(mapped_keys[match(symbols, names(mapped_keys))])
  missing <- is.na(mapped) | !nzchar(mapped)
  .warn_unmapped_gene_ids(symbols[missing], length(symbols), "symbols")

  as.integer(mapped[!missing])
}

#' Convert Entrez identifiers to gene symbols
#'
#' Maps human or mouse Entrez identifiers through the corresponding
#' Bioconductor organism annotation package. Unmapped identifiers produce a
#' warning that names at most the first 20 failures and are dropped. Mapped
#' values retain their input order and are named by the input Entrez identifier
#' so they can be aligned back to another result.
#'
#' @param entrez Integer, numeric, or character vector of Entrez identifiers.
#' @param species Character(1): `"human"`, `"Homo sapiens"`, or `"hsa"`;
#'   alternatively `"mouse"`, `"Mus musculus"`, or `"mmu"`.
#' @return A character vector of mapped gene symbols, named by the input Entrez
#'   identifiers. Unmapped inputs are omitted and the relative input order is
#'   preserved.
#' @examples
#' if (requireNamespace("AnnotationDbi", quietly = TRUE) &&
#'     requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   entrez_to_gene(c(7157L, 1956L))
#' }
#' @export
entrez_to_gene <- function(entrez, species = "human") {
  sp <- .species(species)
  if ((!is.character(entrez) && !is.numeric(entrez)) || anyNA(entrez)) {
    stop("`entrez` must be an integer, numeric, or character vector of ",
         "non-missing Entrez identifiers.", call. = FALSE)
  }
  if (is.numeric(entrez) && any(!is.finite(entrez))) {
    stop("`entrez` must contain only finite Entrez identifiers.", call. = FALSE)
  }
  if (is.numeric(entrez) && any(entrez != floor(entrez))) {
    stop("Numeric `entrez` values must be whole Entrez identifiers.",
         call. = FALSE)
  }
  entrez <- if (is.numeric(entrez)) {
    format(entrez, scientific = FALSE, trim = TRUE)
  } else {
    as.character(entrez)
  }
  if (any(!nzchar(entrez))) {
    stop("`entrez` must not contain empty Entrez identifiers.", call. = FALSE)
  }
  if (!length(entrez)) return(character(0L))

  .require_pkg(
    "AnnotationDbi", "`entrez_to_gene()`",
    'BiocManager::install("AnnotationDbi")'
  )
  .require_pkg(
    sp$orgdb,
    sprintf("`entrez_to_gene(species = \"%s\")`", sp$common),
    sprintf('BiocManager::install("%s")', sp$orgdb)
  )
  db <- getExportedValue(sp$orgdb, sp$orgdb)

  keys <- unique(entrez)
  mapped_keys <- suppressMessages(AnnotationDbi::mapIds(
    db,
    keys = keys,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  ))
  mapped <- as.character(mapped_keys[match(entrez, names(mapped_keys))])
  missing <- is.na(mapped) | !nzchar(mapped)
  .warn_unmapped_gene_ids(entrez[missing], length(entrez), "Entrez IDs")

  stats::setNames(mapped[!missing], entrez[!missing])
}

#' Confounder-gene matching rules
#'
#' This is the single source for both the accepted `drop` categories and their
#' matching rules in [filter_confounder_genes()], so the public choices and the
#' filter implementation cannot drift apart. Regexes use human capitalization;
#' the caller applies them case-insensitively so they also match mouse symbols.
#'
#' @return A named list. Each category contains either a `pattern` regular
#'   expression or an explicit `symbols` character vector.
#' @keywords internal
.confounder_gene_patterns <- function() {
  list(
    ribosomal = list(pattern = "^RP[SL][0-9]+"),
    mito = list(pattern = "^MT-"),
    # The rule as documented in prose was `^HB[ABDEG][0-9]*`, which matches
    # HBEGF -- a growth factor, not a globin. Anchoring the digits, and
    # requiring a hyphen for the mouse spelling, excludes it. The character
    # class also gains M, Q and Z so that mu, theta and zeta globin stop
    # escaping a filter named after them.
    hemoglobin = list(pattern = "^HB[ABDEGMQZ][0-9]*$|^HB[ABDEGMQZ]-"),
    cell_cycle = list(symbols = c(
      "MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "CDC20", "BIRC5", "TYMS"
    )),
    sex = list(symbols = c("XIST", "RPS4Y1", "DDX3Y", "KDM5D"))
  )
}

#' Remove common confounder genes from a symbol vector
#'
#' Removes selected categories of genes that commonly dominate expression
#' signatures for technical or broadly shared biological reasons. Pattern
#' matching and explicit-symbol matching are both case-insensitive, so one rule
#' set covers human uppercase symbols such as `MT-ND1` and mouse title-case
#' symbols such as `mt-Nd1`. That is why there is no `species` argument: an
#' argument that only validated its own value would suggest the rules changed
#' with it. A message reports the number of input elements matched by every
#' active category whenever at least one element is removed.
#'
#' @param symbols Character vector of human or mouse gene symbols.
#' @param drop Character vector selecting any of `"ribosomal"`, `"mito"`,
#'   `"hemoglobin"`, `"cell_cycle"`, and `"sex"`.
#' @return A character vector containing the unfiltered symbols in their
#'   original order.
#' @examples
#' filter_confounder_genes(
#'   c("TP53", "RPS3", "MT-ND1", "MKI67"),
#'   drop = c("ribosomal", "mito", "cell_cycle")
#' )
#' @export
filter_confounder_genes <- function(
    symbols,
    drop = c("ribosomal", "mito", "hemoglobin", "cell_cycle", "sex")) {
  if (!is.character(symbols) || anyNA(symbols) || any(!nzchar(symbols))) {
    stop("`symbols` must be a character vector of non-missing, non-empty ",
         "gene symbols.", call. = FALSE)
  }

  rules <- .confounder_gene_patterns()
  if (!is.character(drop) || anyNA(drop) || any(!nzchar(drop))) {
    stop("`drop` must be a character vector containing zero or more of: ",
         paste0("\"", names(rules), "\"", collapse = ", "), ".",
         call. = FALSE)
  }
  unknown <- setdiff(unique(drop), names(rules))
  if (length(unknown)) {
    stop("Unknown `drop` categor", if (length(unknown) == 1L) "y" else "ies",
         ": ", paste0("\"", unknown, "\"", collapse = ", "), ". Valid ",
         "categories are: ",
         paste0("\"", names(rules), "\"", collapse = ", "), ".",
         call. = FALSE)
  }

  active <- unique(drop)
  counts <- stats::setNames(integer(length(active)), active)
  remove <- rep(FALSE, length(symbols))
  for (category in active) {
    rule <- rules[[category]]
    if (!is.null(rule$pattern)) {
      matched <- grepl(rule$pattern, symbols, ignore.case = TRUE, perl = TRUE)
    } else {
      matched <- toupper(symbols) %in% toupper(rule$symbols)
    }
    counts[[category]] <- sum(matched)
    remove <- remove | matched
  }

  if (any(remove)) {
    for (category in active) {
      message("filter_confounder_genes(): ", category, " removed ",
              counts[[category]], " gene(s).")
    }
  }

  symbols[!remove]
}
