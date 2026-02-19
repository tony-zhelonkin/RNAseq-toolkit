annotate_genes_from_ensembl <- function(ens_ids, try_biomart = TRUE) {
  strip_version <- function(x) sub("\\..*$","", x)
  stable <- strip_version(ens_ids)

  SYMBOL   <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = stable, column = "SYMBOL",   keytype = "ENSEMBL", multiVals = "first")
  ENTREZID <- AnnotationDbi::mapIds(org.Mm.eg.db, keys = stable, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

  ann <- tibble::tibble(
    Symbol       = ifelse(!is.na(SYMBOL) & SYMBOL != "", unname(SYMBOL), stable),
    Ensembl      = stable,
    ENTREZID     = unname(ENTREZID),
    gene_biotype = NA_character_
  )

  if (try_biomart) {
    bt <- try({
      suppressPackageStartupMessages(library(biomaRt))
      ensembl <- biomaRt::useEnsembl("genes", dataset = "mmusculus_gene_ensembl")
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
          # prefer mgi_symbol over SYMBOL if present
          Symbol       = dplyr::coalesce(.data$mgi_symbol, .data$Symbol)
        ) %>%
        dplyr::select(Symbol, Ensembl, ENTREZID, gene_biotype)
    }
  }
  ann
}