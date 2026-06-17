suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(stringr)
  library(dplyr)
})

# 1) Read featureCounts-like wide table (Geneid + samples...) or generic (gene_id).
#    Also handles Salmon gene-level output (gene_id + gene_name + sample columns).
#    Returns a numeric matrix with an attribute:
#      attr(mat, "input_gene_name") — named character vector (name = stripped stable ID,
#        value = gene_name) for Salmon gene-level input; NA for featureCounts/generic.
read_counts_matrix <- function(fp) {
  dt <- fread(fp)
  # Accept common shapes:
  # - featureCounts-processed: Geneid + sample columns
  # - featureCounts raw: columns 1..6 metadata, samples from 7..N
  # - Salmon gene-level: gene_id + gene_name + sample columns
  # - generic: gene_id + sample columns (no gene_name)
  if ("Geneid" %in% names(dt) && !"Chr" %in% names(dt)) {
    gene_col <- "Geneid"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, setdiff(names(dt), gene_col), with = FALSE])
    input_gene_name <- NA
  } else if (all(c("Geneid","Chr","Start","End","Strand","Length") %in% names(dt))) {
    gene_col <- "Geneid"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, (7):ncol(dt), with = FALSE])
    input_gene_name <- NA
  } else if (all(c("gene_id","gene_name") %in% names(dt))) {
    # Salmon gene-level: gene_id + gene_name + sample columns
    gene_col <- "gene_id"
    rn <- dt[[gene_col]]
    sample_cols <- setdiff(names(dt), c("gene_id", "gene_name"))
    mat <- as.matrix(dt[, sample_cols, with = FALSE])
    # Build input_gene_name: first gene_name per stripped stable id (handles versioned ids)
    stripped <- sub("\\..*$", "", rn)
    gn_dt <- data.table(stripped_id = stripped, gene_name = dt[["gene_name"]])
    first_gn <- gn_dt[, .(gene_name = gene_name[1L]), by = stripped_id]
    input_gene_name <- setNames(first_gn$gene_name, first_gn$stripped_id)
  } else if ("gene_id" %in% names(dt)) {
    # Generic: gene_id without gene_name
    gene_col <- "gene_id"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, setdiff(names(dt), gene_col), with = FALSE])
    input_gene_name <- NA
  } else {
    stop("Unsupported counts file format for: ", fp)
  }
  storage.mode(mat) <- "numeric"
  # clean column names: strip paths and extensions
  colnames(mat) <- basename(colnames(mat))
  colnames(mat) <- str_remove(colnames(mat), "\\.bam$|\\.sam$|\\.sorted$|\\.markdup$|\\.txt$")
  rownames(mat) <- rn
  attr(mat, "input_gene_name") <- input_gene_name
  mat
}

# 2) Read metadata (Excel or CSV) and standardize the Sample_ID column.
#
#   fp                   — path to an .xlsx or .csv file
#   sample_col_candidates — column name(s) tried in order for the sample-ID column
#   required_cols        — opt-in character vector of additional column names that MUST
#                          be present (stop() if any are missing); empty by default so
#                          non-project-specific sheets load without error
read_metadata <- function(fp,
                          sample_col_candidates = c("Sample_ID", "Sample ID"),
                          required_cols = character(0)) {
  # Dispatch on extension
  ext <- tolower(tools::file_ext(fp))
  if (ext %in% c("xlsx", "xls")) {
    md <- readxl::read_xlsx(fp)
  } else {
    # .csv and any other plain-text format
    md <- as.data.frame(data.table::fread(fp), check.names = FALSE)
  }
  nm <- names(md)

  # Resolve sample-ID column
  sample_col <- NULL
  for (cand in sample_col_candidates) {
    if (cand %in% nm) { sample_col <- cand; break }
  }
  if (is.null(sample_col)) {
    stop("Metadata must have one of: ", paste(sample_col_candidates, collapse = ", "))
  }
  names(md)[names(md) == sample_col] <- "Sample_ID"

  # Enforce only the caller-specified required columns
  if (length(required_cols)) {
    missing <- setdiff(required_cols, names(md))
    if (length(missing)) stop("Metadata missing required columns: ", paste(missing, collapse = ", "))
  }

  as.data.frame(md, check.names = FALSE)
}

# 3) Align metadata rows to counts columns; subset to shared samples
align_metadata_to_counts <- function(md, counts_cols) {
  extra_md   <- setdiff(md$Sample_ID, counts_cols)
  missing_md <- setdiff(counts_cols, md$Sample_ID)
  if (length(missing_md)) stop("Samples in counts missing in metadata: ", paste(missing_md, collapse=", "))
  if (length(extra_md))  message("[info] Metadata has extra samples (ignored): ", paste(extra_md, collapse=", "))

  md2 <- md[match(counts_cols, md$Sample_ID), , drop = FALSE]
  stopifnot(identical(md2$Sample_ID, counts_cols))
  md2
}

# 4) Write “annotated wide” matrix with top metadata rows
#    add_cols: a data.frame of annotation columns (same row order as mat rows)
write_annotated_matrix <- function(mat, md, add_cols, outfile) {
  stopifnot(is.matrix(mat))
  samples <- colnames(mat)
  stopifnot(identical(samples, md$Sample_ID))

  # Map metadata rows (order matches your Excel header names)
  top_map <- list(
    treat1                = md[["Treatment 1"]],
    treat1_exposure_ante  = md[["Duration of Treatment 1 before Treatment 2"]],
    treat1_exposure_total = md[["Total Duration of Treatment 1"]],
    treat2                = md[["Treatment 2"]],
    treat2_exposure       = md[["Duration of Treatment 2"]],
    bio_replicate         = md[["Biological Replicate (mouse)"]],
    batch                 = md[["Batch"]]
  )

  # Body = annotation columns + counts
  stopifnot(nrow(add_cols) == nrow(mat))
  body <- cbind(add_cols, as.data.frame(mat, check.names = FALSE))

  # Build top rows with the SAME annotation columns as 'add_cols'
  annot_cols <- colnames(add_cols)

  top_block <- do.call(
    rbind,
    lapply(names(top_map), function(k) {
      # blank row for all annotation columns
      annot_row <- as.list(setNames(rep("", length(annot_cols)), annot_cols))
      # put the row label into the first annotation column (assumed "Symbol" if present)
      first_col <- annot_cols[1]
      annot_row[[first_col]] <- k
      # bind metadata vector across sample columns
      df <- as.data.frame(c(annot_row, as.list(top_map[[k]])), check.names = FALSE)
      colnames(df) <- c(annot_cols, samples)
      df
    })
  )

  # Combine and write
  out <- rbind(top_block, body)
  data.table::fwrite(out, outfile, sep = "\t", quote = FALSE, na = "")
  message("Wrote: ", outfile, " [rows: ", nrow(out), ", cols: ", ncol(out), "]")
}
