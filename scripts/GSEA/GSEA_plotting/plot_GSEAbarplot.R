GSEA_barplot <- function(
    gsea_obj, 
    title = "GSEA NES Barplot", 
    q_cut = 0.05, top_n = 30,
    width = 10, height = 8,
    output_dir = "3_Results/imgs/GSEA/",
    save_plot = TRUE,
    dpi = 300) {
    # Extract significant pathways
    gsea_data <- as.data.frame(gsea_obj@result) %>%
        filter(qvalue < q_cut) %>%
        arrange(desc(abs(NES))) %>%
        head(top_n) %>%
        arrange(NES)
    
    p <- ggplot(gsea_data, aes(x = reorder(Description, NES), y = NES)) +
        geom_bar(stat = "identity", aes(fill = NES > 0)) +
        scale_fill_manual(values = c("skyblue", "orange")) +
        coord_flip() +
        labs(title = title,
             x = NULL,
             y = "Normalized Enrichment Score") +
        custom_minimal_theme_with_grid() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(lineheight = 0.8)
        )
    
    # Save the plot if save_plot is TRUE
    if (save_plot) {
        # Create directory if it doesn't exist
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        
        # Create filename from title
        filename <- file.path(output_dir, paste0(gsub(" ", "_", title), "_barplot.pdf"))
        
        ggsave(
            filename = filename,
            plot = p,
            width = width,
            height = height,
            dpi = dpi
        )
    }
    
    return(p)
}