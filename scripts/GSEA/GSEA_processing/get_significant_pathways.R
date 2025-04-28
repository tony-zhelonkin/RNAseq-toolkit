get_significant_pathways <- function(gsea_results, q_cutoff = 0.01) {
    significant_paths <- lapply(gsea_results, function(x) {
        if(!is.null(x)) {
            x$ID[x$p.adjust < q_cutoff]
        }
    })
    unique(unlist(significant_paths))
}