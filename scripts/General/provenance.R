# provenance.R — gene-index writer and session provenance recorder
#
# Functions:
#   write_gene_index(ann_df, outfile)         — G1: rename Symbol/Ensembl cols and write TSV
#   write_session_provenance(outfile, ...)    — G6: sessionInfo + genome build + Ensembl release

suppressPackageStartupMessages(library(data.table))


# write_gene_index -------------------------------------------------------------
# Renames `Ensembl → ensembl_gene_id` and `Symbol → SYMBOL` from the output of
# annotate_genes_from_ensembl(), then writes a tab-separated gene dictionary with
# exactly these columns in order:
#   ensembl_gene_id, SYMBOL, mgi_symbol, ENTREZID, gene_biotype, input_gene_name
#
# Parameters:
#   ann_df  — data.frame/tibble with columns Symbol, Ensembl, mgi_symbol,
#             ENTREZID, gene_biotype, input_gene_name  (as returned by
#             annotate_genes_from_ensembl())
#   outfile — output path; written as a tab-delimited file (.tsv)
write_gene_index <- function(ann_df, outfile) {
  required <- c("Symbol", "Ensembl", "mgi_symbol", "ENTREZID", "gene_biotype", "input_gene_name")
  missing  <- setdiff(required, names(ann_df))
  if (length(missing)) {
    stop("write_gene_index: ann_df is missing column(s): ", paste(missing, collapse = ", "))
  }

  out <- data.table::as.data.table(ann_df)
  data.table::setnames(out,
                       old = c("Ensembl", "Symbol"),
                       new = c("ensembl_gene_id", "SYMBOL"))

  # Enforce canonical column order
  col_order <- c("ensembl_gene_id", "SYMBOL", "mgi_symbol", "ENTREZID",
                 "gene_biotype", "input_gene_name")
  out <- out[, col_order, with = FALSE]

  data.table::fwrite(out, file = outfile, sep = "\t", na = "NA")
  invisible(outfile)
}


# write_session_provenance -----------------------------------------------------
# Writes a plain-text provenance file containing:
#   - genome build (when supplied)
#   - Ensembl annotation release (when supplied)
#   - resolved biomaRt archive release (when a biomaRt session is available at call time)
#   - full sessionInfo() block
#
# Parameters:
#   outfile         — path to write the provenance text
#   genome_build    — character; e.g. "mm10" (from config project.genome_build)
#   ensembl_version — character/integer; Ensembl release passed to useEnsembl() (NULL = floating)
write_session_provenance <- function(outfile, genome_build = NULL, ensembl_version = NULL) {
  lines <- character(0)

  lines <- c(lines, paste0("Provenance recorded: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))

  if (!is.null(genome_build)) {
    lines <- c(lines, paste0("Genome build: ", genome_build))
  }

  if (!is.null(ensembl_version)) {
    lines <- c(lines, paste0("Ensembl version (requested): ", ensembl_version))
  }

  # If biomaRt is loaded, attempt to record the resolved archive release.
  # sessionInfo() only records the biomaRt *package* version, not the live remote
  # archive release — so we query listEnsemblArchives() to capture the actual release.
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    arch_release <- tryCatch({
      archives <- biomaRt::listEnsemblArchives()
      # The row with current_release == TRUE (or "1") is the live release
      current_row <- archives[!is.na(archives$current_release) &
                                archives$current_release %in% c(TRUE, "1", "true"), ]
      if (nrow(current_row) > 0) {
        paste0(current_row$version[1], " (", current_row$date[1], ")")
      } else {
        # Fall back to first row (most recent) when flag is absent
        paste0(archives$version[1], " (", archives$date[1], ")")
      }
    }, error = function(e) {
      paste0("unavailable — ", conditionMessage(e))
    })
    lines <- c(lines, paste0("Ensembl archive release (resolved): ", arch_release))
  }

  lines <- c(lines, "", "--- sessionInfo ---")
  si <- utils::capture.output(utils::sessionInfo())
  lines <- c(lines, si)

  writeLines(lines, con = outfile)
  invisible(outfile)
}
