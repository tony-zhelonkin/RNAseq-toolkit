source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/GSEA_barplot.R")

# Function to run GSEA and create plots for one DE table
run_gsea_analysis <- function(
    de_table, 
    timepoint, 
    n_cat = 15,
    output_dir = "3_Results/imgs/GSEA/") {
  
  # Define database-specific plot parameters
  db_plot_params <- list(
      HALLMARK = list(width = 10, height = 7, font.size = 10),
      REACTOME = list(width = 20, height = 11, font.size = 8),
      KEGG = list(width = 12, height = 8, font.size = 9),
      GOBP = list(width = 13, height = 8, font.size = 9)
  )
  
  # Create the output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
  }
  
  # A list to store all GSEA results by database
  result_list <- list()
  
  # Loop through each database config
  for (db_name in names(db_configs)) {
      
      params <- db_plot_params[[db_name]]
      
      # Run GSEA
      gsea_result <- runGSEA(
          DE_results = de_table,
          rank_metric = "t",
          category = db_configs[[db_name]]$category,
          subcategory = db_configs[[db_name]]$subcategory,
          padj_method = "fdr",
          nperm = 100000,
          pvalue_cutoff = 0.05
      )
      
      # Store the GSEA result in our output list
      result_list[[db_name]] <- gsea_result
      
      # Generate Upregulated (NES>0) dotplot
      GSEA_dotplot(
          gsea_result,
          filterBy = "NES_positive",
          sortBy = "GeneRatio",
          font.size = params$font.size,
          showCategory = n_cat,
          q_cut = 0.05,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          min.dotSize = 2,
          title = paste0(timepoint, " ", db_name, " Upregulated"),
          output_dir = output_dir,
          width = params$width,
          height = params$height
      )
      
      # Generate Downregulated (NES<0) dotplot
      GSEA_dotplot(
          gsea_result,
          filterBy = "NES_negative",
          sortBy = "GeneRatio",
          font.size = params$font.size,
          showCategory = n_cat,
          q_cut = 0.05,
          replace_ = TRUE,
          capitalize_1 = FALSE,
          capitalize_all = FALSE,
          min.dotSize = 2,
          title = paste0(timepoint, " ", db_name, " Downregulated"),
          output_dir = output_dir,
          width = params$width,
          height = params$height
      )
      
      # Generate NES barplot
      GSEA_barplot(
          gsea_result,
          top_n = n_cat * 3,
          title = paste0(timepoint, " ", db_name, " NES"),
          width = params$width,
          height = params$height,
          output_dir = output_dir
      )
  }
  
  # Return the list of GSEA results for this timepoint
  return(result_list)
}