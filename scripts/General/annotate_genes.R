annotate_genes_from_ensembl <- function(ens_ids, try_biomart = TRUE,
                                        input_gene_name = NULL,
                                        biomart_version = NULL,
                                        biomart_host    = NULL) {
  strip_version <- function(x) sub("\\..*$","", x)
  stable <- strip_version(ens_ids)

  SYMBOL   <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = stable, column = "SYMBOL",   keytype = "ENSEMBL", multiVals = "first")
  ENTREZID <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = stable, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

  ann <- tibble::tibble(
    Symbol       = ifelse(!is.na(SYMBOL) & SYMBOL != "", unname(SYMBOL), stable),
    Ensembl      = stable,
    mgi_symbol   = ifelse(!is.na(SYMBOL) & SYMBOL != "", unname(SYMBOL), stable),  # fallback; overwritten below if biomaRt succeeds
    ENTREZID     = unname(ENTREZID),
    gene_biotype = NA_character_
  )

  if (try_biomart) {
    bt <- try({
      suppressPackageStartupMessages(library(biomaRt))
      # Build useEnsembl() call — pass version/host only when non-NULL
      biomart_args <- list("genes", dataset = "mmusculus_gene_ensembl")
      if (!is.null(biomart_version)) biomart_args[["version"]] <- biomart_version
      if (!is.null(biomart_host))    biomart_args[["host"]]    <- biomart_host
      ensembl <- do.call(biomaRt::useEnsembl, biomart_args)
      biomaRt::getBM(attributes = c("ensembl_gene_id","gene_biotype","mgi_symbol"),
                     filters = "ensembl_gene_id", values = unique(stable), mart = ensembl) %>%
        dplyr::as_tibble()
    }, silent = TRUE)

    if (!inherits(bt, "try-error")) {
      ann <- ann %>%
        dplyr::left_join(
          bt %>% dplyr::distinct(ensembl_gene_id, gene_biotype, mgi_symbol),
          by = c("Ensembl" = "ensembl_gene_id")
        ) %>%
        dplyr::mutate(
          # prefer biomaRt biotype when available
          gene_biotype = dplyr::coalesce(.data$gene_biotype.y, .data$gene_biotype.x),
          # keep biomaRt mgi_symbol; fall back to org.db Symbol when missing
          mgi_symbol   = dplyr::coalesce(.data$mgi_symbol.y, .data$mgi_symbol.x),
          # Symbol column: prefer mgi_symbol from biomaRt, fall back to org.db Symbol
          Symbol       = dplyr::coalesce(.data$mgi_symbol.y, .data$Symbol)
        ) %>%
        dplyr::select(Symbol, Ensembl, mgi_symbol, ENTREZID, gene_biotype)
    }
  }

  # Thread input_gene_name: join by stripped stable id, NA when arg is NULL
  if (!is.null(input_gene_name)) {
    # input_gene_name must be a named vector (names = stripped stable ids) or
    # an unnamed vector parallel to ens_ids; normalise to named vector keyed on stable
    if (is.null(names(input_gene_name))) {
      ign_named <- stats::setNames(input_gene_name, stable)
    } else {
      ign_named <- input_gene_name
    }
    ann <- ann %>%
      dplyr::mutate(input_gene_name = unname(ign_named[.data$Ensembl]))
  } else {
    ann <- ann %>%
      dplyr::mutate(input_gene_name = NA_character_)
  }

  # Ensure final column order: Symbol, Ensembl, mgi_symbol, ENTREZID, gene_biotype, input_gene_name
  ann %>% dplyr::select(Symbol, Ensembl, mgi_symbol, ENTREZID, gene_biotype, input_gene_name)
}