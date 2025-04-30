#' Run Simplified GSEA Analysis and Visualization Pipeline
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
#'
#' @return List of GSEA results for each database
#' @export
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
    pvalue_cutoff = 0.05
) {
    ############################################################################
    ## 1. Source helper functions                                             ##
    ############################################################################
    base_dir <- "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/Kelsey_Followup"
    paths <- c(
        "R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea.R",
        "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R",
        "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
        "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
        "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_barplot.R",
        "R_GSEA_visualisations/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
        "R_GSEA_visualisations/scripts/custom_minimal_theme.R"
    )
    for (p in paths) {
        source(file.path(base_dir, p))
    }

    if (save_plots && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    ############################################################################
    ## 2. Set up database list                                                ##
    ############################################################################
    if (is.null(databases)) {
        databases <- list(
            hallmark = list(name = "Hallmark",                db_species = "MM", collection = "H",  subcollection = ""),
            gobp     = list(name = "GO Biological Process",   db_species = "MM", collection = "C5", subcollection = "GO:BP"),
            gomf     = list(name = "GO Molecular Function",   db_species = "MM", collection = "C5", subcollection = "GO:MF"),
            gocc     = list(name = "GO Cellular Component",   db_species = "MM", collection = "C5", subcollection = "GO:CC"),
            kegg     = list(name = "KEGG",                    db_species = "MM", collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
            reactome = list(name = "Reactome",                db_species = "MM", collection = "C2", subcollection = "CP:REACTOME"),
            biocarta = list(name = "BioCarta",                db_species = "MM", collection = "C2", subcollection = "CP:BIOCARTA"),
            wiki     = list(name = "WikiPathways",            db_species = "MM", collection = "C2", subcollection = "CP:WIKIPATHWAYS")
        )
    }

    ############################################################################
    ## 3. Process each database                                               ##
    ############################################################################
    gsea_results <- list()

    for (db_name in names(databases)) {
        cfg <- databases[[db_name]]
        message("▶  Processing ", cfg$name)

        # Get database-specific plot parameters
        plot_params <- get_db_plot_params(db_name)
        
        # Run GSEA for this database
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
                NULL
            })
        
        gsea_results[[db_name]] <- gsea_res
        if (is.null(gsea_res) || nrow(gsea_res@result) == 0L) {
            message("No results for ", db_name, ", skipping")
            next
        }
        
        if (save_plots) {
            # Create database-specific output directory
            db_dir <- file.path(output_dir, db_name)
            if (!dir.exists(db_dir)) {
                dir.create(db_dir, recursive = TRUE)
            }
            
            message("DEBUG: Starting up-regulated dotplot for ", db_name)
            # Generate and save up-regulated dotplot
            p_up <- gsea_dotplot(
                gsea_res,
                filterBy = "NES_positive",
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s Up", analysis_name, db_name),
                wrap_width = 50,
                pos_color = "#fc8d59", 
                neg_color = "#91bfdb"
            )
            
            save_gsea_plot(
                p_up,
                filename = sprintf("%s_%s_up_dot.pdf", analysis_name, db_name),
                width = plot_params$width,
                height = plot_params$height,
                base_font_size = plot_params$font_size,
                dir = db_dir
            )
            
            message("DEBUG: Starting down-regulated dotplot for ", db_name)
            # Generate and save down-regulated dotplot
            p_down <- gsea_dotplot(
                gsea_res,
                filterBy = "NES_negative",
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s Down", analysis_name, db_name),
                wrap_width = 50,
                pos_color = "#fc8d59", 
                neg_color = "#91bfdb"
            )
            
            save_gsea_plot(
                p_down,
                filename = sprintf("%s_%s_down_dot.pdf", analysis_name, db_name),
                width = plot_params$width,
                height = plot_params$height,
                base_font_size = plot_params$font_size,
                dir = db_dir
            )
            
            message("DEBUG: Starting facetted dotplot for ", db_name)
            # Generate and save faceted dotplot
            p_facet <- gsea_dotplot_facet(
                gsea_res,
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s", analysis_name, db_name),
                wrap_width = 50,
                pos_color = "#fc8d59", 
                neg_color = "#91bfdb"
            )
            
            save_gsea_plot(
                p_facet,
                filename = sprintf("%s_%s_facet.pdf", analysis_name, db_name),
                width = plot_params$width,
                height = plot_params$height * 1.4, # Slightly taller for faceted plot
                base_font_size = plot_params$font_size,
                dir = db_dir
            )
            
            message("DEBUG: Starting NES barplot for ", db_name)
            # Generate and save NES barplot
            p_bar <- gsea_barplot(
                gsea_res,
                padj_cutoff = padj_cutoff,
                top_n = n_pathways,
                title = sprintf("%s %s NES", analysis_name, db_name),
                pos_color = "#fc8d59", 
                neg_color = "#91bfdb"
            )
            
            save_gsea_plot(
                p_bar,
                filename = sprintf("%s_%s_nes_bar.pdf", analysis_name, db_name),
                width = plot_params$width,
                height = plot_params$height,
                base_font_size = plot_params$font_size,
                dir = db_dir
            )

            message("DEBUG: Starting running sum plots for ", db_name)
            # Generate running sum plots 
            top5   <- order(abs(gsea_res@result$NES), decreasing = TRUE)[1:5]

            p_run  <- gsea_running_sum_plot(
                        gsea_res,
                        gene_set_ids = top5,
                        base_size    = plot_params$font_size
            )

            save_gsea_plot(
                p_run,
                filename = sprintf("%s_%s_running_sum.pdf", analysis_name, db_name),
                width    = plot_params$width,
                height   = plot_params$height * 1.2,   # a bit taller
                base_font_size = plot_params$font_size,
                dir      = db_dir
            )


        }
    }

    message("✓  GSEA analysis finished.")
    invisible(gsea_results)
}
