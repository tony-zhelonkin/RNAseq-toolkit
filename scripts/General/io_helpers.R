suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(stringr)
  library(dplyr)
})

# 1) Read featureCounts-like wide table (Geneid + samples...) or generic (gene_id)
read_counts_matrix <- function(fp) {
  dt <- fread(fp)
  # Accept common shapes:
  # - featureCounts-processed: Geneid + sample columns
  # - featureCounts raw: columns 1..6 metadata, samples from 7..N
  if ("Geneid" %in% names(dt) && !"Chr" %in% names(dt)) {
    gene_col <- "Geneid"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, setdiff(names(dt), gene_col), with = FALSE])
  } else if (all(c("Geneid","Chr","Start","End","Strand","Length") %in% names(dt))) {
    gene_col <- "Geneid"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, (7):ncol(dt), with = FALSE])
  } else if ("gene_id" %in% names(dt)) {
    gene_col <- "gene_id"
    rn <- dt[[gene_col]]
    mat <- as.matrix(dt[, setdiff(names(dt), gene_col), with = FALSE])
  } else {
    stop("Unsupported counts file format for: ", fp)
  }
  storage.mode(mat) <- "numeric"
  # clean column names: strip paths and extensions
  colnames(mat) <- basename(colnames(mat))
  colnames(mat) <- str_remove(colnames(mat), "\\.bam$|\\.sam$|\\.sorted$|\\.markdup$|\\.txt$")
  rownames(mat) <- rn
  mat
}

# 2) Read Excel metadata and standardize key fields
read_metadata <- function(xlsx_fp) {
  md <- readxl::read_xlsx(xlsx_fp)
  nm <- names(md)

  # Sample ID column name could be "Sample_ID" or "Sample ID"
  sample_col <- if ("Sample_ID" %in% nm) "Sample_ID" else if ("Sample ID" %in% nm) "Sample ID" else stop("Metadata must have 'Sample_ID' or 'Sample ID'.")
  names(md)[names(md) == sample_col] <- "Sample_ID"

  # Normalize expected fields (present in your sheet)
  # Keep original names, but we’ll reference these exact ones:
  #  Treatment 1
  #  Duration of Treatment 1 before Treatment 2
  #  Total Duration of Treatment 1
  #  Treatment 2
  #  Duration of Treatment 2
  #  Biological Replicate (mouse)
  #  Batch
  needed <- c("Sample_ID",
              "Treatment 1","Duration of Treatment 1 before Treatment 2","Total Duration of Treatment 1",
              "Treatment 2","Duration of Treatment 2",
              "Biological Replicate (mouse)","Batch")
  missing <- setdiff(needed, names(md))
  if (length(missing)) stop("Metadata missing columns: ", paste(missing, collapse=", "))

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
