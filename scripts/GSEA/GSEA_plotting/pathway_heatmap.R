# Modified heatmap function with controlled column order
plot_pathway_heatmap <- function(scores, title, annotation_col) {
    pheatmap(t(scores),  
             scale = "row",
             clustering_method = "complete",
             clustering_distance_rows = "correlation",
             cluster_cols = FALSE,  # Disable column clustering
             show_rownames = TRUE,
             show_colnames = FALSE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("navy", "white", "red"))(50),
             main = title,
             gaps_col = c(12, 24),  # Add visual separators between main groups
             colnames = sample_order)  # Set column order
}