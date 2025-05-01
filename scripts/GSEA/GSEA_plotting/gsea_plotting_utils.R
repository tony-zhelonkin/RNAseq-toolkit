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
save_gsea_plot <- function(plot, filename,
                           width, height,
                           base_font_size = 10,
                           dir  = NULL,
                           dpi  = 300)
{
  if (!is.null(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(dir, filename)
  }

  ## ----- dynamic font scaling (unchanged) --------------------------
  scale_fac <- sqrt((width * height) / (7 * 5))
  plot <- plot +
    theme(
      text        = element_text(size = base_font_size * scale_fac),
      axis.title  = element_text(size = base_font_size * scale_fac),
      axis.text   = element_text(size = base_font_size * scale_fac * 0.9),
      plot.title  = element_text(size = base_font_size * scale_fac * 1.2),
      legend.title= element_text(size = base_font_size * scale_fac * 0.9),
      legend.text = element_text(size = base_font_size * scale_fac * 0.8),
      strip.text  = element_text(size = base_font_size * scale_fac)
    )

  ## ----- open -> print -> close ------------------------------------
  grDevices::pdf(filename, width = width, height = height,
                 family = "sans", pointsize = base_font_size)
  on.exit(grDevices::dev.off(), add = TRUE)   # <── GUARANTEED CLOSE

  print(plot)                                 # write page
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
  if (length(words) <= 1) return(text)
  
  total_chars <- nchar(text)
  
  if (total_chars > width * 1.5) {
    # Very long text: split into three parts
    third_point <- ceiling(length(words) / 3)
    two_thirds <- third_point * 2
    
    part1 <- paste(words[1:third_point], collapse = " ")
    part2 <- paste(words[(third_point+1):two_thirds], collapse = " ")
    part3 <- paste(words[(two_thirds+1):length(words)], collapse = " ")
    
    return(paste(part1, part2, part3, sep = "\n"))
  } else if (total_chars > width) {
    # Moderately long text: split in half
    mid_point <- ceiling(length(words) / 2)
    
    part1 <- paste(words[1:mid_point], collapse = " ")
    part2 <- paste(words[(mid_point+1):length(words)], collapse = " ")
    
    return(paste(part1, part2, sep = "\n"))
  }
  
  return(text)
}
