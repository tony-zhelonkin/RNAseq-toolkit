###############################################################################
##  Shared Plotting Utilities for RNA-seq Analysis
##  Provides DRY helper functions to eliminate code duplication
###############################################################################

#' Ensure directory exists
#'
#' Creates directory if it doesn't exist, with recursive parent creation
#' and suppressed warnings. Idempotent operation.
#'
#' @param path Character vector of directory paths to create
#' @return NULL (invisibly)
#' @examples
#' ensure_dir("/path/to/output/plots")
#' ensure_dir(here::here("03_Results", "Analysis", "Plots"))
ensure_dir <- function(path) {
  if (length(path) == 0 || is.null(path)) {
    stop("ensure_dir: path cannot be NULL or empty")
  }

  for (p in path) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
  }

  invisible(NULL)
}

#' Save plot to PDF with consistent settings
#'
#' Wrapper around PDF device creation that handles device management,
#' error recovery, and ensures devices are properly closed.
#'
#' @param plot ggplot object or expression to evaluate for base plots
#' @param filepath Full path to output PDF file (including extension)
#' @param width Plot width in inches (default: 7)
#' @param height Plot height in inches (default: 7)
#' @param ensure_dir_exists Logical, create parent directory if needed (default: TRUE)
#' @return NULL (invisibly)
#' @examples
#' p <- ggplot(mtcars, aes(mpg, wt)) + geom_point()
#' save_plot(p, "output/myplot.pdf", width = 8, height = 6)
save_plot <- function(plot, filepath, width = 7, height = 7,
                      ensure_dir_exists = TRUE) {

  if (ensure_dir_exists) {
    parent_dir <- dirname(filepath)
    ensure_dir(parent_dir)
  }

  # Close any existing devices to prevent accumulation
  if (dev.cur() > 1) {
    tryCatch(dev.off(), error = function(e) NULL)
  }

  # Open PDF device
  pdf(filepath, width = width, height = height)

  # Ensure device is closed even if error occurs
  on.exit({
    if (dev.cur() > 1) {
      dev.off()
    }
  }, add = TRUE)

  # Print plot
  if (inherits(plot, "ggplot")) {
    print(plot)
  } else if (inherits(plot, "pheatmap")) {
    # pheatmap objects need special handling
    grid::grid.draw(plot$gtable)
  } else {
    # For base plots or other objects, just evaluate
    plot
  }

  invisible(NULL)
}

#' Load checkpoint RDS file
#'
#' Convenience wrapper for loading checkpoint files with error handling
#' and consistent messaging.
#'
#' @param name Name of checkpoint file (without .rds extension)
#' @param dir Directory containing checkpoints (default: checkpoint_dir from parent env)
#' @param required Logical, stop if file missing? (default: TRUE)
#' @return Loaded R object, or NULL if not required and missing
#' @examples
#' fit <- load_checkpoint("fit_object")
#' results <- load_checkpoint("gsea_results", dir = "checkpoints")
load_checkpoint <- function(name, dir = NULL, required = TRUE) {
  # Try to get checkpoint_dir from parent environment if not specified
  if (is.null(dir)) {
    if (exists("checkpoint_dir", envir = parent.frame())) {
      dir <- get("checkpoint_dir", envir = parent.frame())
    } else if (exists("checkpoint_dir", envir = .GlobalEnv)) {
      dir <- checkpoint_dir
    } else {
      stop("load_checkpoint: 'dir' not specified and checkpoint_dir not found in environment")
    }
  }

  # Construct filepath
  if (!grepl("\\.rds$", name)) {
    name <- paste0(name, ".rds")
  }
  filepath <- file.path(dir, name)

  # Check if file exists
  if (!file.exists(filepath)) {
    if (required) {
      stop("Required checkpoint file not found: ", filepath)
    } else {
      warning("Optional checkpoint file not found: ", filepath)
      return(NULL)
    }
  }

  # Load and return
  message("  Loading: ", basename(filepath))
  readRDS(filepath)
}

#' Process all contrasts with topTable
#'
#' Eliminates duplicate topTable() calls by processing all contrasts once
#' and returning a named list of results tables.
#'
#' @param fit limma MArrayLM fit object
#' @param contrasts Contrasts matrix (colnames used as names)
#' @param number Number of genes to return (default: Inf for all genes)
#' @param sort_by Sorting method passed to topTable (default: "none")
#' @return Named list of data frames, one per contrast
#' @examples
#' contrast_tables <- process_all_contrasts(fit, contrasts)
#' # Access individual contrast results
#' g32a_d35 <- contrast_tables[["G32A_vs_Ctrl_D35"]]
process_all_contrasts <- function(fit, contrasts, number = Inf, sort_by = "none") {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("limma package required for process_all_contrasts")
  }

  contrast_names <- colnames(contrasts)

  message("Processing ", length(contrast_names), " contrasts...")

  results <- lapply(contrast_names, function(co) {
    tbl <- limma::topTable(fit, coef = co, number = number, sort.by = sort_by)
    message("  ✓ ", co, ": ", nrow(tbl), " genes")
    tbl
  })

  names(results) <- contrast_names

  message("✓ All contrasts processed")

  results
}

#' Save contrast results tables to CSV
#'
#' Batch export of contrast results to CSV files with consistent naming.
#'
#' @param contrast_tables Named list of data frames from process_all_contrasts()
#' @param output_dir Directory to save CSV files
#' @param suffix Suffix to add to filenames (default: "_results.csv")
#' @return NULL (invisibly)
#' @examples
#' save_contrast_csvs(contrast_tables, "03_Results/02_Analysis/DE_results")
save_contrast_csvs <- function(contrast_tables, output_dir, suffix = "_results.csv") {
  ensure_dir(output_dir)

  message("Saving ", length(contrast_tables), " contrast result tables...")

  for (co in names(contrast_tables)) {
    filepath <- file.path(output_dir, paste0(co, suffix))
    write.csv(contrast_tables[[co]], file = filepath, row.names = TRUE)
  }

  message("✓ CSV files saved to: ", output_dir)

  invisible(NULL)
}

#' Get database-specific plot parameters
#'
#' Returns consistent plot dimensions and font sizes for different databases.
#' Falls back to default values if database not recognized.
#'
#' @param db_name Name of database (e.g., "hallmark", "kegg", "syngo")
#' @return List with width, height, and font_size
#' @examples
#' params <- get_plot_params("syngo")
#' save_plot(my_plot, "output.pdf", width = params$width, height = params$height)
get_plot_params <- function(db_name) {
  # Default parameters
  defaults <- list(width = 10, height = 8, font_size = 11)

  # Database-specific overrides
  params <- list(
    hallmark = list(width = 10, height = 6, font_size = 11),
    kegg = list(width = 10, height = 10, font_size = 10),
    reactome = list(width = 10, height = 12, font_size = 10),
    gobp = list(width = 10, height = 14, font_size = 9),
    gocc = list(width = 10, height = 10, font_size = 10),
    gomf = list(width = 10, height = 10, font_size = 10),
    wiki = list(width = 10, height = 12, font_size = 10),
    syngo = list(width = 10, height = 10, font_size = 10),
    cgp = list(width = 10, height = 12, font_size = 10),
    tf = list(width = 10, height = 12, font_size = 10),
    canon = list(width = 10, height = 10, font_size = 10)
  )

  if (db_name %in% names(params)) {
    return(params[[db_name]])
  } else {
    return(defaults)
  }
}

#' Format contrast name for display
#'
#' Converts underscore-separated contrast names to more readable format.
#'
#' @param contrast_name Character vector of contrast names
#' @param replace_underscore Replacement for underscores (default: " ")
#' @return Character vector of formatted names
#' @examples
#' format_contrast_name("G32A_vs_Ctrl_D35")
#' # Returns: "G32A vs Ctrl D35"
format_contrast_name <- function(contrast_name, replace_underscore = " ") {
  gsub("_", replace_underscore, contrast_name)
}

#' Close all graphics devices safely
#'
#' Closes all open graphics devices, useful for cleanup after plotting loops.
#'
#' @return NULL (invisibly)
#' @examples
#' close_all_devices()
close_all_devices <- function() {
  while (dev.cur() > 1) {
    tryCatch(dev.off(), error = function(e) NULL)
  }
  invisible(NULL)
}

#' Message with timestamp
#'
#' Print message with timestamp prefix for better log tracking.
#'
#' @param ... Arguments passed to message()
#' @param timestamp Logical, include timestamp? (default: TRUE)
#' @return NULL (invisibly)
#' @examples
#' log_message("Starting analysis...")
log_message <- function(..., timestamp = TRUE) {
  if (timestamp) {
    ts <- format(Sys.time(), "[%H:%M:%S]")
    message(ts, " ", ...)
  } else {
    message(...)
  }
  invisible(NULL)
}

###############################################################################
##  Package-level exports
###############################################################################

# If used as a package, export these functions
if (exists(".packageName")) {
  export(ensure_dir, save_plot, load_checkpoint, process_all_contrasts,
         save_contrast_csvs, get_plot_params, format_contrast_name,
         close_all_devices, log_message)
}

message("✓ utils_plotting.R loaded successfully")
