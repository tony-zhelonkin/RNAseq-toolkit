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
#' @param pvalue_cutoff P-value cutoff for storing GSEA results (default: 1, stores ALL pathways).
#'        Set to 1 to retain all pathways in results for downstream analysis.
#'        Display filtering is controlled separately by padj_cutoff.
#'
#' @return List of GSEA results for each database
#' @export

# R_GSEA_visualisations/scripts/GSEA/GSEA_processing/run_gsea_analysis.R
# ---------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x # null-coalescing operator

run_gsea_analysis <- function(
    de_table,
    analysis_name,
    rank_metric = "t",
    species = "Mus musculus",
    n_pathways = 30,
    padj_cutoff = 0.05,
    save_plots = TRUE,
    output_dir = "./GSEA_Plots",
    databases = NULL,
    nperm = 100000,
    pvalue_cutoff = 1, # save all pathways for downstream visualisations
    sample_annotation = NULL,
    sample_order = NULL,
    helper_root = NULL # <── NEW  (default = NULL)
    ) {
    # ------------------------------------------------------------------ #
    # 1.  locate helper scripts                                          #
    # ------------------------------------------------------------------ #
    if (!is.null(helper_root)) {
        base_dir <- helper_root # <-- use caller-supplied root
    } else {
        this_file <- attr(body(run_gsea_analysis), "srcfile")$filename
        base_dir <- if (!is.null(this_file)) {
            dirname(dirname(dirname(this_file)))
        } else {
            getwd()
        }
    }

    helper_paths <- c(
        "/scripts/custom_minimal_theme.R",
        "/scripts/GSEA/GSEA_plotting/gsea_plotting_utils.R",
        "/scripts/GSEA/GSEA_plotting/format_pathway_names.R",  # Load formatter first
        "/scripts/GSEA/GSEA_plotting/gsea_dotplot.R",
        "/scripts/GSEA/GSEA_plotting/gsea_dotplot_facet.R",
        "/scripts/GSEA/GSEA_plotting/gsea_barplot.R",
        "/scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R",
        "/scripts/GSEA/GSEA_plotting/gsea_heatmap.R",
        "/scripts/GSEA/GSEA_processing/run_gsea.R",
        "/scripts/DE/volcano_helpers.R"
    )

    for (hp in helper_paths) {
        full <- file.path(base_dir, hp) # define first …
        message("[DEBUG] sourcing ", full) # … then print
        if (file.exists(full)) {
            source(full)
        } else {
            message("[run_gsea_analysis] helper not found → ", full)
        }
    }
    # ––––– 2  default DB list  ---------------------------------------------------
    # infer msigdbr species code from full species name 
    dbsp <- if (grepl("sapiens", species, ignore.case = TRUE)) "HS" else "MM"

    if (is.null(databases)) {
            databases <- list(
            hallmark = list(name = "Hallmark", db_species = dbsp, collection = "H", subcollection = ""),
            canon    = list(name = "Canonical Pathways", db_species = dbsp, collection = "C2", subcollection = "CP"),
            gobp     = list(name = "GO BP", db_species = dbsp, collection = "C5", subcollection = "GO:BP"),
            gomf     = list(name = "GO MF", db_species = dbsp, collection = "C5", subcollection = "GO:MF"),
            gocc     = list(name = "GO CC", db_species = dbsp, collection = "C5", subcollection = "GO:CC"),
            kegg     = list(name = "KEGG", db_species = dbsp, collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
            reactome = list(name = "Reactome", db_species = dbsp, collection = "C2", subcollection = "CP:REACTOME"),
            wiki     = list(name = "WikiPath", db_species = dbsp, collection = "C2", subcollection = "CP:WIKIPATHWAYS"),
            cgp     = list(name = "Chem-Genetic Perturbations", db_species = dbsp, collection = "C2", subcollection = "CGP"),
            tf      = list(name = "GTRD", db_species = dbsp, collection = "C3", subcollection = "TFT:GTRD")
        )
    }

    # ––––– 3  loop through DBs  --------------------------------------------------
    gsea_results <- list()

    for (db_name in names(databases)) {
        cfg <- databases[[db_name]]
        message("▶  ", cfg$name)

        res <- tryCatch(
            run_gsea(
                DE_results = de_table,
                rank_metric = rank_metric,
                species = species,
                db_species = cfg$db_species,
                collection = cfg$collection,
                subcollection = cfg$subcollection,
                nperm = nperm,
                pvalue_cutoff = pvalue_cutoff
            ),
            error = function(e) {
                warning(e)
                NULL
            }
        )

        gsea_results[[db_name]] <- res
        if (is.null(res) || nrow(res@result) == 0L) next

        if (!save_plots) next
        db_dir <- file.path(output_dir, db_name)
        dir.create(db_dir, showWarnings = FALSE)

        # dot-plots / barplot / running-sum  (unchanged ↓)
        plot_par <- get_db_plot_params(db_name)

        save_gsea_plot(
            gsea_dotplot(
                gsea_obj = res,
                filterBy = "NES_positive",
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s Up", analysis_name, db_name)
            ),
            filename = sprintf("%s_%s_up_dot.pdf", analysis_name, db_name),
            width = plot_par$width,
            height = plot_par$height,
            base_font_size = plot_par$font_size,
            dir = db_dir
        )

        save_gsea_plot(
            gsea_dotplot(
                gsea_obj = res,
                filterBy = "NES_negative",
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s Down", analysis_name, db_name)
            ),
            filename = sprintf("%s_%s_down_dot.pdf", analysis_name, db_name),
            width = plot_par$width,
            height = plot_par$height,
            base_font_size = plot_par$font_size,
            dir = db_dir
        )

        save_gsea_plot(
            gsea_dotplot_facet(
                gsea_obj = res,
                showCategory = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s", analysis_name, db_name)
            ),
            filename = sprintf("%s_%s_facet.pdf", analysis_name, db_name),
            width = plot_par$width,
            height = plot_par$height * 1.4,
            base_font_size = plot_par$font_size,
            dir = db_dir
        )

        save_gsea_plot(
            gsea_barplot(
                gsea_obj = res,
                top_n = n_pathways,
                padj_cutoff = padj_cutoff,
                title = sprintf("%s %s NES", analysis_name, db_name)
            ),
            filename = sprintf("%s_%s_nes_bar.pdf", analysis_name, db_name),
            width = plot_par$width,
            height = plot_par$height,
            base_font_size = plot_par$font_size,
            dir = db_dir
        )

        top5 <- order(abs(res@result$NES), decreasing = TRUE)[1:5]
        save_gsea_plot(
            gsea_running_sum_plot(
                gsea_obj = res,
                gene_set_ids = top5,
                base_size = plot_par$font_size
            ),
            filename = sprintf("%s_%s_running_sum.pdf", analysis_name, db_name),
            width = plot_par$width,
            height = plot_par$height * 1.2,
            base_font_size = plot_par$font_size,
            dir = db_dir
        )

        # ––––– NEW: per-database sample × pathway heat-map  -----------------------
        if (!is.null(sample_annotation)) {
            # 1) build a sample × pathway matrix of NES (0 if not sign.)
            top_tbl <- res@result |>
                dplyr::filter(p.adjust < as.numeric(padj_cutoff)) |>
                dplyr::arrange(p.adjust) |>
                utils::head(n_pathways)
                
            # Check if we have any significant pathways
            if (nrow(top_tbl) > 0) {
                geneset <- top_tbl$ID
                nes_vec <- top_tbl$NES
                names(nes_vec) <- geneset
                
                # replicate the vector for every sample (one contrast → same NES),
                # but keep dimension n_sample × n_pathways so pheatmap works.
                sample_ids <- rownames(sample_annotation)
                mat <- matrix(rep(nes_vec, each = length(sample_ids)),
                    nrow = length(sample_ids), byrow = TRUE,
                    dimnames = list(sample_ids, geneset)
                )

                # reorder columns (samples) if requested
                if (!is.null(sample_order)) {
                    mat <- mat[sample_order, , drop = FALSE]
                }

                tryCatch({
                    gsea_heatmap_save(mat,
                        file = file.path(
                            db_dir,
                            sprintf("%s_%s_heatmap.pdf", analysis_name, db_name)
                        ),
                        annotation_col = sample_annotation[rownames(mat), , drop = FALSE],
                        ann_colors = NULL,
                        main = sprintf("%s – %s", analysis_name, cfg$name),
                        gaps_col = NULL,
                        cluster_cols = FALSE
                    ) # keep your order
                }, error = function(e) {
                    warning("Error generating heatmap for ", db_name, ": ", e$message)
                })
            } else {
                message("No significant pathways found for ", db_name, " with padj_cutoff = ", padj_cutoff)
            }
        }
    }

    invisible(gsea_results)
}
