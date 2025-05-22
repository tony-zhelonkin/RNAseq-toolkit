#' Save GSEA Plot with Specified Dimensions
#'
#' Central function to save GSEA plots with consistent dimensions and scaling.
#'
#' @param plot ggplot object to save
#' @param filename Output filename (with or without path)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param base_font_size Base font size to use (scales with dimensions)
#' @param dir Output directory (optional)
#' @param dpi Resolution for raster outputs
#'
#' @return Invisibly returns the original plot object
#' @export
# -------------------------------------------------------------------
# save_gsea_plot()
#   • opens its own pdf() device
#   • ALWAYS closes that device via on.exit()
#   • returns the (possibly re-themed) plot invisibly
# -------------------------------------------------------------------
save_gsea_plot <- function(
    plot, 
    filename,
    width, 
    height,
    base_font_size = 10,
    dir = NULL,
    dpi = 300) {
  if (!is.null(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(dir, filename)
  }

  ## ----- dynamic font scaling (unchanged) --------------------------
  scale_fac <- sqrt((width * height) / (7 * 5))
  plot <- plot +
    theme(
      text = element_text(size = base_font_size * scale_fac),
      axis.title = element_text(size = base_font_size * scale_fac),
      axis.text = element_text(size = base_font_size * scale_fac * 0.9),
      plot.title = element_text(size = base_font_size * scale_fac * 1.2),
      legend.title = element_text(size = base_font_size * scale_fac * 0.9),
      legend.text = element_text(size = base_font_size * scale_fac * 0.8),
      strip.text = element_text(size = base_font_size * scale_fac)
    )

  ## ----- open -> print -> close ------------------------------------
  grDevices::pdf(filename,
    width = width, height = height,
    family = "sans", pointsize = base_font_size
  )
  on.exit(grDevices::dev.off(), add = TRUE) # <── GUARANTEED CLOSE

  print(plot) # write page
  invisible(plot)
}

#' Get Standard Plot Parameters by Database
#'
#' Returns predefined width, height, and font size parameters based on the database type.
#'
#' @param db_name Name of the database (e.g., "hallmark", "gobp", "kegg")
#' @return List with width, height, and font_size elements
#' @export
get_db_plot_params <- function(db_name) {
  # Lowercase for consistency
  db_lower <- tolower(db_name)

  # Default parameters
  default_params <- list(width = 7, height = 5, font_size = 10)

  # Database-specific parameters
  db_params <- list(
    hallmark = list(width = 7, height = 5, font_size = 10),
    gobp = list(width = 8, height = 8, font_size = 9),
    gomf = list(width = 8, height = 7, font_size = 9),
    gocc = list(width = 7, height = 6, font_size = 9),
    kegg = list(width = 8, height = 7, font_size = 9),
    reactome = list(width = 9, height = 8, font_size = 8),
    biocarta = list(width = 7, height = 6, font_size = 9),
    wiki = list(width = 8, height = 7, font_size = 9)
  )

  # Return database-specific parameters or default if not found
  return(db_params[[db_lower]] %||% default_params)
}

#' Simplified version of smart text wrapping function
#'
#' @param text Text to wrap
#' @param width Maximum character width
#' @return Wrapped text with newlines
smart_wrap <- function(text, width = 40) {
  words <- unlist(strsplit(text, " "))
  if (length(words) <= 1) {
    return(text)
  }

  total_chars <- nchar(text)

  if (total_chars > width * 1.5) {
    # Very long text: split into three parts
    third_point <- ceiling(length(words) / 3)
    two_thirds <- third_point * 2

    part1 <- paste(words[1:third_point], collapse = " ")
    part2 <- paste(words[(third_point + 1):two_thirds], collapse = " ")
    part3 <- paste(words[(two_thirds + 1):length(words)], collapse = " ")

    return(paste(part1, part2, part3, sep = "\n"))
  } else if (total_chars > width) {
    # Moderately long text: split in half
    mid_point <- ceiling(length(words) / 2)

    part1 <- paste(words[1:mid_point], collapse = " ")
    part2 <- paste(words[(mid_point + 1):length(words)], collapse = " ")

    return(paste(part1, part2, sep = "\n"))
  }

  return(text)
}

#' Save GSEA results to a text log file
#'
#' @param gsea_obj GSEA result object
#' @param filename Output filename (with or without path)
#' @param padj_cutoff Adjusted p-value cutoff
#' @param dir Output directory (optional)
#'
#' @return Invisibly returns TRUE if successful
#' @export
save_gsea_log <- function(
    gsea_obj,
    filename,
    padj_cutoff = 0.05,
    dir = NULL
) {
  if (!is.null(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(dir, filename)
  }
  
  # Extract results
  results <- as.data.frame(gsea_obj@result)
  
  # Open connection to file
  con <- file(filename, "w")
  on.exit(close(con), add = TRUE)  # Ensure file is closed when function exits
  
  # Write header information
  cat("GSEA Results Log\n", file = con)
  cat("================\n\n", file = con)
  cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", file = con)
  
  # Write summary statistics
  cat("Summary Statistics:\n", file = con)
  cat("-----------------\n", file = con)
  cat("Total pathways analyzed:", nrow(results), "\n", file = con)
  
  sig_up <- sum(results$p.adjust < padj_cutoff & results$NES > 0)
  sig_down <- sum(results$p.adjust < padj_cutoff & results$NES < 0)
  
  cat("Significant upregulated pathways (padj <", padj_cutoff, "):", sig_up, "\n", file = con)
  cat("Significant downregulated pathways (padj <", padj_cutoff, "):", sig_down, "\n", file = con)
  cat("Total significant pathways:", sig_up + sig_down, "\n\n", file = con)
  
  # If no significant results
  if (sig_up + sig_down == 0) {
    cat("No significant pathways were found at adjusted p-value cutoff of", padj_cutoff, "\n\n", file = con)
    return(invisible(TRUE))
  }
  
  # Write detailed results for significant pathways
  cat("Significant Pathways:\n", file = con)
  cat("-------------------\n\n", file = con)
  
  # Upregulated pathways
  if (sig_up > 0) {
    cat("UPREGULATED PATHWAYS:\n", file = con)
    cat("=====================\n\n", file = con)
    
    up_results <- results[results$p.adjust < padj_cutoff & results$NES > 0, ]
    up_results <- up_results[order(up_results$p.adjust), ]
    
    for (i in 1:nrow(up_results)) {
      pathway <- up_results[i, ]
      cat(i, ". ", pathway$Description, "\n", file = con)
      cat("   NES: ", sprintf("%.3f", pathway$NES), "\n", file = con)
      cat("   Adjusted p-value: ", sprintf("%.3e", pathway$p.adjust), "\n", file = con)
      cat("   Leading edge size: ", pathway$setSize, "\n", file = con)
      cat("   Leading edge genes: ", pathway$core_enrichment, "\n\n", file = con)
    }
  }
  
  # Downregulated pathways
  if (sig_down > 0) {
    cat("DOWNREGULATED PATHWAYS:\n", file = con)
    cat("=======================\n\n", file = con)
    
    down_results <- results[results$p.adjust < padj_cutoff & results$NES < 0, ]
    down_results <- down_results[order(down_results$p.adjust), ]
    
    for (i in 1:nrow(down_results)) {
      pathway <- down_results[i, ]
      cat(i, ". ", pathway$Description, "\n", file = con)
      cat("   NES: ", sprintf("%.3f", pathway$NES), "\n", file = con)
      cat("   Adjusted p-value: ", sprintf("%.3e", pathway$p.adjust), "\n", file = con)
      cat("   Leading edge size: ", pathway$setSize, "\n", file = con)
      cat("   Leading edge genes: ", pathway$core_enrichment, "\n\n", file = con)
    }
  }
  
  cat("End of GSEA Results Log\n", file = con)
  
  return(invisible(TRUE))
}
