##############################
### 1. Prepare environment ###
#############################

# Load libraries and functions
source("1_Scripts/Load_libraries.R")
source("1_Scripts/process_rnaseq_data.R")
source("1_Scripts/PCA.R")
source("1_Scripts/standard_volcano.R")
source("1_Scripts/GSEA_dotplot.R")
source("1_Scripts/runGSEA.R")
source("1_Scripts/runningSumGSEAplot.R")
source("1_Scripts/combined_volcano.R")
source("1_Scripts/runGSEA_pool.R")

#####################
### 2. Read data ###
#####################

# Prepare data
## Read count data
counts <- read.delim("0_Data/STAR_Data/counts.txt", row.names = 1)

## Read sample information
sample_info <- read.delim("0_Data/STAR_Data/seq_reference.txt", row.names = 1)
colnames(sample_info) <- c("Group", "Time", "Treatment", "RANKL")
sample_info$Sample <- rownames(sample_info)

## Create DGEList object with sample information
DGErankl <- process_rnaseq_data(counts, sample_info)
# Display sample information
#DT::datatable(DGErankl$samples)

# PCA
create_pca_plot(DGErankl, title = "PCA Plot")


#####################
### 3. Pool GSEA ###
#####################
# Set design matrix
## Set combinatorial design
design <- model.matrix(~ 0 + group, data = DGErankl$samples)
## Rename the columns to remove "group" prefix
colnames(design) <- gsub("group", "", colnames(design))

## Set contrast matrix
### Pool contrast design
contrasts.p <- makeContrasts(
    # STm effect without RANKL at 4h
    STm_4h_0 = t4h_STm_0 - t4h_mock_0,
    # STm effect with RANKL at 4h
    STm_4h_100 = t4h_STm_100 - t4h_mock_100,
    # RANKL effect without infection at 4h
    RANKL_4h_mock = t4h_mock_100 - t4h_mock_0,
    # RANKL effect with infection at 4h
    RANKL_4h_STm = t4h_STm_100 - t4h_STm_0,
    # STm effect without RANKL at 24h
    STm_24h_0 = t24h_STm_0 - t24h_mock_0,
    # STm effect with RANKL at 24h
    STm_24h_100 = t24h_STm_100 - t24h_mock_100,
    # RANKL effect without infection at 24h
    RANKL_24h_mock = t24h_mock_100 - t24h_mock_0,
    # RANKL effect with infection at 24h
    RANKL_24h_STm = t24h_STm_100 - t24h_STm_0,
    # Design
    levels = design
)

### Select contrast desing
contrasts.m <- makeContrasts(
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    # RANKL time-dependent effect
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) - ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    # Design
    levels = design
)

## Fit the model
fit <- edgeR::voomLmFit(
    counts = DGErankl,
    design = design,
    sample.weights = TRUE
)
## Compute contrasts from the linear model fit
fit <- contrasts.fit(fit, contrasts.p)
## Compute Bayes moderated t-statistics and p-values
fit <- eBayes(fit, robust = TRUE)

##############
#### GSEA ####
# Run pooled GSEA analysis
pooled_gsea_results <- run_pooled_gsea(fit, contrasts.p, DGErankl)

sample_order <- with(DGErankl$samples, {
  # We see 'rankl' is numeric: 0 or 100
  # 'treatment' is "mock" or "STm"
  # 'timepoint' is numeric: 4 or 24

  rankl_levels <- c(0, 100)
  treatment_levels <- c("mock", "STm")
  time_levels <- c(4, 24)

  ordered_samples <- character()
  for(r in rankl_levels) {
    for(tr in treatment_levels) {
      for(t in time_levels) {
        # Create a logical index to find rows that match r, tr, t
        idx <- (rankl == r & treatment == tr & timepoint == t)
        these_cells <- rownames(DGErankl$samples)[idx]
        ordered_samples <- c(ordered_samples, these_cells)
      }
    }
  }
  ordered_samples
})

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

# Calculate normalized expression scores
norm_counts <- cpm(DGErankl, log = TRUE)

# Add sample annotations
annotation_col <- data.frame(
    Time = DGErankl$samples$timepoint,
    Treatment = DGErankl$samples$treatment,
    RANKL = DGErankl$samples$rankl,
    row.names = colnames(norm_counts)
)

# Create heatmaps with annotations
# Define annotation colors
ann_colors <- list(
    timepoint = c("4h" = "#1f77b4", "24h" = "#ff7f0e"),
    treatment = c("mock" = "#2ca02c", "STm" = "#d62728"),
    rankl = c("0" = "#9467bd", "100" = "#8c564b")
)

# Create and save heatmaps
pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/A.Pair-wise comparisons/hallmark_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$hallmark, "Hallmark Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/A.Pair-wise comparisons/kegg_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$kegg, "KEGG Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/A.Pair-wise comparisons/gobp_heatmap.pdf", 
    width = 12, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$gobp, "GO Biological Processes", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/A.Pair-wise comparisons/reactome_heatmap.pdf", 
    width = 18, 
    height = 8)
plot_pathway_heatmap(pooled_gsea_results$scores$reactome, "REACTOME Pathways", annotation_col)
dev.off()

###################################
### 4. Specific-questions GSEA ###
###################################

### Select contrast desing
contrasts.m <- makeContrasts(
    # RANKL effect on infection response at 4h
    RANKL_effect_4h = (t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0),
    # RANKL effect on infection response at 24h
    RANKL_effect_24h = (t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0),
    # RANKL time-dependent effect
    Time_RANKL_effect = ((t24h_STm_100 - t24h_mock_100) - (t24h_STm_0 - t24h_mock_0)) - ((t4h_STm_100 - t4h_mock_100) - (t4h_STm_0 - t4h_mock_0)),
    # Design
    levels = design
)

# Statistical estimation
## Fit the model
fit <- edgeR::voomLmFit(
    counts = DGErankl,
    design = design,
    sample.weights = TRUE
)
## Compute contrasts from the linear model fit
fit <- contrasts.fit(fit, contrasts.m)
## Compute Bayes moderated t-statistics and p-values
fit <- eBayes(fit, robust = TRUE)
## Extract contrasts coefficients
### Early infection RANKL modulation
DE_rankl_4h <- topTable(fit,
    coef = "RANKL_effect_4h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Late infection RANKL modulation
DE_rankl_24h <- topTable(fit,
    coef = "RANKL_effect_24h",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)
### Time-dependent RANKL modulation of the infection effect
DE_rankl_time <- topTable(fit,
    coef = "Time_RANKL_effect",
    sort.by = "t",
    adjust.method = "fdr",
    n = Inf
)

### Save results
write.csv(DE_rankl_4h, "3_Results/DE_tables/DE_rankl_4h.csv")
write.csv(DE_rankl_24h, "3_Results/DE_tables/DE_rankl_24h.csv")
write.csv(DE_rankl_time, "3_Results/DE_tables/DE_rankl_time.csv")

# Volcano plots 
## Early infection
p1 <- create_volcano_plot(
  DE_rankl_4h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection response change at 4 hours"
)

## Late infection
p2 <- create_volcano_plot(
  DE_rankl_24h,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 8,
  label_method = "sig",
  title = "RANKL+infection response change at 24 hours"
)

## Time-dependent RANKL effect
p3 <- create_volcano_plot(
  DE_rankl_time,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  x_breaks = 2,
  max.overlaps = 10,
  label_method = "sig",
  title = "RANKL+infection DEGs changing over time: 24h vs 4h "
)

## Combined volcano 
p4 <- create_combined_volcano_plots(
    DE_rankl_4h,
    DE_rankl_24h,
    DE_rankl_time,
    p_cutoff = 0.05,
    fc_cutoff = 2.0,
    max.overlaps = 15
)


# GSEA 
## Time-dependent RANKL effect 
### 1. REACTOME
gsea_reactome_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:REACTOME",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)


### 2. HALLMARK
gsea_hallmark_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "H",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

### 3. GO:BP
gsea_gobp_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

### 4. KEGG
gsea_kegg_time <- runGSEA(
  DE_rankl_time,
  rank_metric = "t",
  species = "Mus musculus",
  category = "C2",
  subcategory = "CP:KEGG",
  padj_method = "fdr",
  nperm = 100000,
  pvalue_cutoff = 0.05
)

# GSEA visualisation
## 1. REACTOME
GSEA_dotplot(
  gsea_obj = gsea_reactome_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched REACTOME time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2,
  output_dir = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/"
)

runSumGSEAplot(
  gsea_reactome_time,
  gene_set_ids = 1:5,
  title = "REACTOME running sum for top 5 pathways"
)

## 2. HALLMARK
GSEA_dotplot(
  gsea_obj = gsea_hallmark_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched HALLMARK time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2,
  output_dir = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/"
)

runSumGSEAplot(
  gsea_hallmark_time,
  gene_set_ids = 1:5,
  title = "HALLMARK running sum for top 5 pathways"
)

## 3. GO:BP
GSEA_dotplot(
  gsea_obj = gsea_gobp_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched GOBP time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2,
  output_dir = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/"
)

runSumGSEAplot(
  gsea_gobp_time,
  gene_set_ids = 1:5,
  title = "GOBP running sum for top 5 pathways"
)

## 4. KEGG 
GSEA_dotplot(
  gsea_obj = gsea_kegg_time,
  showCategory = 15,
  font.size = 10,
  title = "Top 15 KO enriched KEGG time-dependent gene sets",
  sortBy = "GeneRatio",
  filterBy = "NES",
  q_cut = 0.05,
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 2,
  output_dir = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/"
)

runSumGSEAplot(
  gsea_kegg_time,
  gene_set_ids = 1:5,
  title = "KEGG running sum for top 5 pathways"
)

# Get pathway genes for each database
hallmark_genes <- get_pathway_genes(gsea_hallmark_time)
kegg_genes <- get_pathway_genes(gsea_kegg_time)
gobp_genes <- get_pathway_genes(gsea_gobp_time)
reactome_genes <- get_pathway_genes(gsea_reactome_time)

# Calculate pathway scores
hallmark_scores <- calculate_pathway_scores(norm_counts, hallmark_genes)
kegg_scores <- calculate_pathway_scores(norm_counts, kegg_genes)
gobp_scores <- calculate_pathway_scores(norm_counts, gobp_genes)
reactome_scores <- calculate_pathway_scores(norm_counts, reactome_genes)

# Create and save heatmaps
pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/hallmark_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(hallmark_scores, "Hallmark Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/kegg_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(kegg_scores, "KEGG Pathways", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/gobp_heatmap.pdf", width = 12, height = 8)
plot_pathway_heatmap(gobp_scores, "GO Biological Processes", annotation_col)
dev.off()

pdf("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/B.Granular_questions/reactome_heatmap.pdf", width = 20, height = 8)
plot_pathway_heatmap(reactome_scores, "REACTOME Pathways", annotation_col)
dev.off()


#---------------------#---------------------#---------------------#-------------
# Follow-up analysis №1
### 5. Which pathways RANKL is (up or) downregulating in STm infected cells? 
###  i.e. RANKL effect on infection, 4h or 24h individually  
#---------------------#---------------------#---------------------#-------------

source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/run_GSEA_analysis.R")


# Define databases to analyze
db_configs <- list(
    HALLMARK = list(category = "H", subcategory = NULL),
    REACTOME = list(category = "C2", subcategory = "CP:REACTOME"),
    KEGG = list(category = "C2", subcategory = "CP:KEGG"),
    GOBP = list(category = "C5", subcategory = "GO:BP")
)

# Run analysis for both timepoints with custom output directory
custom_output_dir <- "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/"

# Run analysis for both timepoints
res_4h <- run_gsea_analysis(DE_rankl_4h, "4h", output_dir = custom_output_dir)
res_24h <- run_gsea_analysis(DE_rankl_24h, "24h", output_dir = custom_output_dir)

str(res_4h)

#---------------------#---------------------#---------------------#-------------
# Follow-up analysis №2
# - separate NES neg and pos plots 
# - pathway gene heatmaps
#---------------------#---------------------#---------------------#-------------

source("/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/1_Scripts/plot_single_pathway_heatmap.R", encoding = "UTF-8")

# Assuming you already have:
# gsea_result_4h_kegg (the GSEA object for 4h KEGG)
# norm_counts (expression matrix)
# sample_order (vector of sample names in correct order)
# annotation_col, ann_colors

# 1) 4h KEGG downregulated: TOLL LIKE RECEPTOR SIGNALING PATHWAY
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["KEGG"]],
  pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  expression_data = norm_counts,
  sample_order = sample_order,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_KEGG_downregulated"
)

# 2) 4h GO BP downregulated pathways
# We might have gsea_result_4h_gobp = runGSEA(...) for 4h GO:BP
for(path_name in c(
  "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  "GOBP_CELLULAR_RESPONSE_TO_BIOTIC_STIMULUS",
  "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN"
)) {
  plot_single_pathway_heatmap(
    gsea_obj = res_4h[["GOBP"]],
    pathway_name = path_name,
    expression_data = norm_counts,
    sample_order = sample_order,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_GOBP_downregulated"
  )
}

# 3) 4h GO BP upregulated: GOBP NEGATIVE REGULATION OF UBIQUITIN PROTEIN LIGASE ACTIVITY
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["GOBP"]],
  pathway_name = "GOBP_NEGATIVE_REGULATION_OF_UBIQUITIN_PROTEIN_LIGASE_ACTIVITY",
  expression_data = norm_counts,
  sample_order = sample_order,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_GOBP_upregulated"
)

## Follow-up 3

# ## TOLL LIKE RECEPTOR SIGNALING PATHWAY
# # 1) 4h KEGG downregulated: TOLL LIKE RECEPTOR SIGNALING PATHWAY
# plot_single_pathway_heatmap(
#   gsea_obj = res_4h[["KEGG"]],
#   pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
#   expression_data = norm_counts,
#   sample_order = sample_order,
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_TLR_downregulated"
# )

# # 2) 24h KEGG downregulated: TOLL LIKE RECEPTOR SIGNALING PATHWAY
# plot_single_pathway_heatmap(
#   gsea_obj = res_24h[["KEGG"]],
#   pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
#   expression_data = norm_counts,
#   sample_order = sample_order,
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/24h_TLR_downregulated"
# )

# ## CELLULAR RESPONSE TO MOLECULE OF BACTERIAL ORIGIN
# # 1) 4h GOBP downregulated: CELLULAR RESPONSE TO MOLECULE OF BACTERIAL ORIGIN
# plot_single_pathway_heatmap(
#   gsea_obj = res_4h[["GOBP"]],
#   pathway_name = "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
#   expression_data = norm_counts,
#   sample_order = sample_order,
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_Bacterial_response_downregulated"
# )

# # 2) 24h GOBP downregulated: CELLULAR RESPONSE TO MOLECULE OF BACTERIAL ORIGIN
# plot_single_pathway_heatmap(
#   gsea_obj = res_24h[["GOBP"]],
#   pathway_name = "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
#   expression_data = norm_counts,
#   sample_order = sample_order,
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/24h_Bacterial_response_downregulated"
# )


#### Separate the 4 and 24h column samples 

# First, create a time-based sample order
time_based_sample_order <- with(DGErankl$samples, {
  # First get all 4h samples
  time_4h_samples <- character()
  for(r in c(0, 100)) {
    for(tr in c("mock", "STm")) {
      idx <- (rankl == r & treatment == tr & timepoint == 4)
      these_cells <- rownames(DGErankl$samples)[idx]
      time_4h_samples <- c(time_4h_samples, these_cells)
    }
  }
  
  # Then get all 24h samples
  time_24h_samples <- character()
  for(r in c(0, 100)) {
    for(tr in c("mock", "STm")) {
      idx <- (rankl == r & treatment == tr & timepoint == 24)
      these_cells <- rownames(DGErankl$samples)[idx]
      time_24h_samples <- c(time_24h_samples, these_cells)
    }
  }
  
  # Combine them
  c(time_4h_samples, time_24h_samples)
})


# For TLR signaling pathway
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["KEGG"]],
  pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  expression_data = norm_counts,
  sample_order = time_based_sample_order,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/Combined_TLR",
  gaps_col = 12
)

# For cellular response to bacterial origin
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["GOBP"]],
  pathway_name = "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  expression_data = norm_counts,
  sample_order = time_based_sample_order,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/Combined_Bacterial_response",
  gaps_col = 12
)



### 4h and 24h separate 
# Get 4h samples only
samples_4h <- rownames(DGErankl$samples)[DGErankl$samples$timepoint == 4]

# Get 24h samples only
samples_24h <- rownames(DGErankl$samples)[DGErankl$samples$timepoint == 24]

# Create 4h-only heatmap for TLR signaling
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["KEGG"]],
  pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  expression_data = norm_counts,
  sample_order = samples_4h,  # Only 4h samples
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_only_TLR"
)

# Create 24h-only heatmap for TLR signaling
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["KEGG"]],
  pathway_name = "KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
  expression_data = norm_counts,
  sample_order = samples_24h,  # Only 24h samples
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/24h_only_TLR"
)

# Similarly for bacterial response pathway
# 4h only
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["GOBP"]],
  pathway_name = "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  expression_data = norm_counts,
  sample_order = samples_4h,  # Only 4h samples
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/4h_only_Bacterial_response"
)

# 24h only
plot_single_pathway_heatmap(
  gsea_obj = res_4h[["GOBP"]],
  pathway_name = "GOBP_CELLULAR_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
  expression_data = norm_counts,
  sample_order = samples_24h,  # Only 24h samples
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  output_prefix = "/Users/tony/My Drive (anton.bioinf.md@gmail.com)/Data_Analysis/ClaraRANKL_STAR/3_Results/imgs/GSEA/C.4_24h_RANKL_separate/Pathways/24h_only_Bacterial_response"
)