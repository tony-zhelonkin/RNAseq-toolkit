GSEA_faceted_dotplot <- function(gsea_obj, showCategory = 10, font.size = 10, title = "GSEA Faceted Dotplot") {
    # Get both positive and negative results
    pos_data <- as.data.frame(gsea_obj@result) %>%
        filter(NES > 0, qvalue < 0.05) %>%
        arrange(desc(NES)) %>%
        head(showCategory)
    
    neg_data <- as.data.frame(gsea_obj@result) %>%
        filter(NES < 0, qvalue < 0.05) %>%
        arrange(NES) %>%
        head(showCategory)
    
    # Combine data
    plot_data <- rbind(pos_data, neg_data) %>%
        mutate(Direction = ifelse(NES > 0, "Upregulated", "Downregulated"))
    
    # Create faceted plot
    p <- ggplot(plot_data, aes(x = abs(NES), y = reorder(Description, NES))) +
        geom_point(aes(size = -log10(qvalue), color = Direction)) +
        scale_color_manual(values = c("Upregulated" = "orange", "Downregulated" = "skyblue")) +
        facet_grid(Direction ~ ., scales = "free_y", space = "free") +
        labs(x = "Normalized Enrichment Score", y = NULL) +
        theme_minimal() +
        theme(text = element_text(size = font.size))
    
    return(p)
}