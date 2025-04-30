#' Run Comprehensive GSEA Analysis and Visualization Pipeline
#'
#' Performs GSEA using `run_gsea` for multiple specified MSigDB databases/collections
#' on a given differential expression table. Generates a standard set of plots
#' (dotplots, barplot, running sum) for each database.
#'
#' @param de_table Differential expression results table
#' @param analysis_name Name of the analysis for plot labeling
#' @param rank_metric Column to use for ranking genes (default: 't')
#' @param species Species for MSigDB gene sets (default: "Mus musculus")
#' @param n_pathways Number of top pathways to display (default: 30)
#' @param padj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @param save_plots Logical, save generated plots (default: TRUE)
#' @param output_dir Directory to save plots
#' @param databases List of databases to analyze (default: NULL, uses predefined set)
#' @param nperm Number of permutations for GSEA (default: 100000)
#' @param pvalue_cutoff P-value cutoff for GSEA (default: 0.05)
#' @param plot_width Width of plots in inches (default: 7)
#' @param plot_height Height of plots in inches (default: 5)
#' @param base_font_size Base font size for plots (default: 8)
#'
#' @return List of GSEA results for each database
#' @export
`%||%` <- function(x, y) if (is.null(x)) y else x   # null‑coalescing operator

run_gsea_analysis <- function(
    de_table,
    analysis_name,
    rank_metric   = "t",
    species       = "Mus musculus",
    n_pathways    = 30,
    padj_cutoff   = 0.05,
    save_plots    = TRUE,
    output_dir    = "./GSEA_Plots",
    databases     = NULL,
    nperm         = 100000,
    pvalue_cutoff = 0.05,
    plot_width    = 7,
    plot_height   = 5,
    base_font_size = 8
) {
    message("DEBUG: Starting run_gsea_analysis")
    ############################################################################
    ## 1. sanity checks                                                        ##
    ############################################################################
    message("DEBUG: Performing sanity checks")
    stopifnot(is.data.frame(de_table))
    stopifnot(rank_metric %in% colnames(de_table))
    message("DEBUG: Sanity checks passed")

    ############################################################################
    ## 2. helper functions                                                     ##
    ############################################################################
    # Bring in the processing + plotting helpers
    source_scripts <- function() {
        message("DEBUG: Sourcing helper scripts")
        base_dir <- "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/Kelsey_Followup"
        paths <- c(
            "R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R",
            "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
            "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
            "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R",
            "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
            "R_GSEA_visualisations/scripts/custom_minimal_theme.R")
        for (p in paths) {
            file <- file.path(base_dir, p)
            if (file.exists(file)) {
                message("DEBUG: Sourcing ", file)
                source(file)
            } else {
                message("DEBUG: File not found: ", file)
            }
        }
        message("DEBUG: Finished sourcing helper scripts")
    }
    source_scripts()

    if (save_plots && !dir.exists(output_dir)) {
        message("DEBUG: Creating output directory: ", output_dir)
        dir.create(output_dir, recursive = TRUE)
    }

    ############################################################################
    ## 3. default database list                                                ##
    ############################################################################
    message("DEBUG: Setting up database list")
    if (is.null(databases)) {
        databases <- list(
            hallmark = list(name = "Hallmark",                db_species = "MM", collection = "H",  subcollection = ""),
            gobp     = list(name = "GO Biological Process",   db_species = "MM", collection = "C5", subcollection = "GO:BP"),
            gomf     = list(name = "GO Molecular Function",   db_species = "MM", collection = "C5", subcollection = "GO:MF"),
            gocc     = list(name = "GO Cellular Component",   db_species = "MM", collection = "C5", subcollection = "GO:CC"),
            kegg     = list(name = "KEGG",                    db_species = "MM", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
            reactome = list(name = "Reactome",                db_species = "MM", collection = "C2", subcollection = "CP:REACTOME"),
            biocarta = list(name = "BioCarta",                db_species = "MM", collection = "C2", subcollection = "CP:BIOCARTA"),
            wiki     = list(name = "WikiPathways",            db_species = "MM", collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
            grtd     = list(name = "Regulatory sets TF M3:GTRD", db_species = "MM", collection = "M3", subcollection = "GTRD")
        )
        message("DEBUG: Using default database list with ", length(databases), " databases")
    } else {
        message("DEBUG: Using provided database list with ", length(databases), " databases")
    }

    ############################################################################
    ## 4. iterate through databases                                            ##
    ############################################################################
    message("DEBUG: Starting database iteration")
    gsea_results <- list()

    for (db_name in names(databases)) {
        cfg <- databases[[db_name]]
        message("▶  ", cfg$name)
        message("DEBUG: Processing database: ", db_name, " (", cfg$name, ")")

        message("DEBUG: Running GSEA for ", db_name)
        gsea_res <- tryCatch(
            run_gsea(
                DE_results     = de_table,
                rank_metric    = rank_metric,
                species        = species,
                db_species     = cfg$db_species,
                collection     = cfg$collection,
                subcollection  = cfg$subcollection %||% "",
                nperm          = nperm,
                pvalue_cutoff  = pvalue_cutoff
            ),
            error = function(e) {
                warning("GSEA failed for ", db_name, ": ", e$message)
                message("DEBUG: GSEA failed for ", db_name, ": ", e$message)
                NULL
            })
        gsea_results[[db_name]] <- gsea_res
        if (is.null(gsea_res) || nrow(gsea_res@result) == 0L) {
            message("DEBUG: No results for ", db_name, ", skipping")
            next
        }
        message("DEBUG: GSEA successful for ", db_name, ", got ", nrow(gsea_res@result), " results")

        # Helper that saves the plot with proper dimensions
        save_plot <- function(plot, file, width, height, dpi = 300) {
            ggplot2::ggsave(
                filename = file,
                plot     = plot,
                width    = width,
                height   = height,
                units    = "in",
                dpi      = dpi,
                limitsize = FALSE
            )
            message("   ↳ saved ", basename(file))
        }


        # Up dotplot
        message("DEBUG: Creating up dotplot")
        p_up <- tryCatch(gsea_dotplot(gsea_res,
                                      filterBy = "NES_positive",
                                      showCategory = n_pathways,
                                      padj_cutoff = padj_cutoff,
                                      base_font_size = base_font_size,
                                      wrap_width = 50,
                                      width = plot_width,
                                      height = plot_height,
                                      title = sprintf("%s %s Up", analysis_name, db_name)),
                         error = function(e) {
                             message("DEBUG: Error creating up dotplot: ", e$message)
                             NULL
                         })
        if (!is.null(p_up)) message("DEBUG: Up dotplot created successfully")

        save_plot(p_up, 
                  file.path(output_dir, sprintf("%s_%s_up_dot.pdf", analysis_name, db_name)),
                  width = plot_width,
                  height = plot_height)

        # Down dotplot
        message("DEBUG: Creating down dotplot")
        p_down <- tryCatch(gsea_dotplot(gsea_res,
                                        filterBy = "NES_negative",
                                        showCategory = n_pathways,
                                        padj_cutoff = padj_cutoff,
                                        base_font_size = base_font_size,
                                        wrap_width = 50,
                                        width = plot_width,
                                        height = plot_height,
                                        title = sprintf("%s %s Down", analysis_name, db_name)),
                           error = function(e) {
                               message("DEBUG: Error creating down dotplot: ", e$message)
                               NULL
                           })
        if (!is.null(p_down)) message("DEBUG: Down dotplot created successfully")

        save_plot(
            p_down, 
            file.path(output_dir, sprintf("%s_%s_down_dot.pdf", analysis_name, db_name)),
            width = plot_width,
            height = plot_height)

        # Facet dotplot
        message("DEBUG: Creating facet dotplot")
        p_facet <- tryCatch(gsea_dotplot_facet(gsea_res,
                                              showCategory = n_pathways,
                                              padj_cutoff = padj_cutoff,
                                              base_font_size = base_font_size,
                                              wrap_width = 50,
                                              width = plot_width,
                                              height = plot_height * 1.4, # Slightly taller for faceted plot
                                              title = sprintf("%s %s", analysis_name, db_name)),
                           error = function(e) {
                               message("DEBUG: Error creating facet dotplot: ", e$message)
                               NULL
                           })
        if (!is.null(p_facet)) message("DEBUG: Facet dotplot created successfully")
        save_plot(
            p_facet, 
            file.path(output_dir, sprintf("%s_%s_facet.pdf", analysis_name, db_name)),
            width = plot_width,
            height = plot_height * 1.4)

        # NES barplot
        message("DEBUG: Creating NES barplot")
        p_bar <- tryCatch(gsea_barplot(gsea_res,
                                      top_n = n_pathways,
                                      padj_cutoff = padj_cutoff,
                                      base_font_size = base_font_size,
                                      width = plot_width,
                                      height = plot_height,
                                      title = sprintf("%s %s NES", analysis_name, db_name)),
                         error = function(e) {
                             message("DEBUG: Error creating NES barplot: ", e$message)
                             NULL
                         })
        if (!is.null(p_bar)) message("DEBUG: NES barplot created successfully")
        save_plot(p_bar, file.path(output_dir, sprintf("%s_%s_nes_bar.pdf", analysis_name, db_name)))

        # Running-sum plot
        message("DEBUG: Starting running sum plot for ", db_name)
        tryCatch({
            # Get top 5 pathways by absolute NES
            idx <- order(abs(gsea_res@result$NES), decreasing = TRUE)[1:min(5, nrow(gsea_res@result))]
            message("DEBUG: Selected indices: ", paste(idx, collapse = ", "))
            message("DEBUG: Selected pathway IDs: ", paste(gsea_res@result$ID[idx], collapse = ", "))
            
            # Create running sum plot
            p_run <- gsea_running_sum_plot(
                gsea_res, 
                gene_set_ids = idx,
                title = sprintf("%s %s Running‑sum", analysis_name, db_name),
                combine = TRUE,
                width = 10,  # Wider for running sum plot
                height = 8   # Taller for running sum plot
            )
            
            if (!is.null(p_run)) {
                message("DEBUG: Running sum plot created successfully")
                save_plot(p_run, file.path(output_dir, sprintf("%s_%s_running_sum.pdf", analysis_name, db_name)))
            } else {
                message("DEBUG: Running sum plot is NULL, not saving")
            }
        }, error = function(e) {
            warning("Running‑sum plot failed for ", db_name, ": ", e$message)
            message("DEBUG: Running sum plot failed: ", e$message)
        })
    }

    message("✓  GSEA analysis finished.")
    invisible(gsea_results)
}
