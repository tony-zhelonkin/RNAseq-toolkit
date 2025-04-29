# R GSEA Visualizations

This repository contains a collection of R functions for comprehensive RNA-seq data analysis, with a focus on Gene Set Enrichment Analysis (GSEA) and visualization. These tools are designed to work with differential expression results from edgeR/limma and GSEA results from the ClusterProfiler R package.

## Overview

The scripts are organized into modular components for:
1. Differential Expression (DE) analysis visualization
2. Gene Set Enrichment Analysis (GSEA) processing and visualization
3. Custom visualization themes

## Script Organization

```
scripts/
├── custom_minimal_theme.R       # Custom ggplot2 theme for consistent visualization
├── DE/                          # Differential Expression visualization
│   ├── analyzePathVolcanoViz.R  # Multi-style pathway volcano plots
│   ├── plot_standard_volcano.R  # Standard volcano plots for DE results
│   ├── plotPCA.R                # 2D PCA plots
│   └── plotPCA3d.R              # 3D PCA plots using plotly
└── GSEA/                        # GSEA analysis and visualization
    ├── GSEA_processing/         # Core GSEA analysis functions
    │   ├── calculate_pathway_scores.R    # Calculate pathway scores from expression data
    │   ├── get_pathway_genes_all.R       # Extract pathway genes across contrasts
    │   ├── get_significant_pathways.R    # Identify significant pathways
    │   ├── runGSEA.R                     # Run GSEA on DE results
    │   └── runGSEA_pool.R                # Run pooled GSEA across multiple contrasts
    └── GSEA_plotting/           # GSEA visualization functions
        ├── dotplot.R                     # Custom GSEA dotplots
        ├── dotplot_compare.R             # Compare GSEA results between two conditions
        ├── nesOutlier.R                  # NES comparison scatter plots
        ├── plot_RunSumGSEA.R             # Running sum enrichment plots
        └── plot_single_pathway_heatmap.R # Heatmaps for specific pathways
```

## Key Functions



## GSEA NESscatter Comparison Plot Function

This repository contains an R function `gsea_comparison_plot` for comparing and visualizing Gene Set Enrichment Analysis (GSEA) results between two datasets.

### Function Description

The `gsea_comparison_plot` function generates a comprehensive comparison of GSEA results from two datasets, including a scatter plot visualization, data preparation, and extraction of common pathways.

#### Parameters

- `x`, `y`: GSEA result dataframes for comparison
- `x_label`, `y_label`: Labels for the two datasets
- `use_normalized`: Boolean to use Normalized Enrichment Score (NES) instead of Enrichment Score (default: TRUE, i.e. use NES)
- `percentile_threshold`: Threshold for determining significant pathways (default: 0.95)
- `color`: Color scheme for the plot ("all", "common", or "distinct")
- `max_overlaps`: Maximum number of label overlaps in the plot (default: 20)

#### Dependencies

The function requires the following R packages:

- dplyr
- ggplot2
- ggrepel
- tidyr

### Usage

```R
source("gsea_comparison_plot.R")

# Assuming you have two GSEA result dataframes named 'gsea_result1' and 'gsea_result2'
comparison_results <- gsea_comparison_plot(gsea_result1, gsea_result2,
                                           x_label = "Dataset1", y_label = "Dataset2",
                                           use_normalized = TRUE,
                                           percentile_threshold = 0.95,
                                           color = "all")

# Access the plot
print(comparison_results$plot)

# Access the results dataframe
head(comparison_results$data)

# Access common pathways
common_pathways <- comparison_results$common_pathways

# Access GMT-formatted dataframes
gmt_data <- comparison_results$gmt
```


### Example usage

```R
res <- gsea_comparison_plot(
  gsea_sum_GOMF_th1,
  gsea_sum_GOMF_th17,
  "Th1",
  "Th17",
  color = "all",
  max_overlaps = Inf,
  percentile_threshold = 0.90,
  use_normalized = TRUE
  )

# Display the plot
print(res$plot)
```

![alt_text](https://github.com/tony-zhelonkin/R_GSEA_visualisations/blob/main/examples/NESoutlier/Enrichment%20Score%20Outlier%20Analysis.png?raw=true)

### Features

1. Generates a scatter plot comparing enrichment scores or NES between two datasets
2. Identifies and categorizes common and distinct pathways
3. Provides flexible color schemes for visualizing results
4. Implements a colorblind-friendly palette
5. Labels top pathways using ggrepel for improved readability
6. Generates GMT-formatted dataframes for further analysis
7. Returns a comprehensive list object with plot, data, common pathways, and GMT-formatted results

### Customization

You can customize the plot and analysis by adjusting various parameters:

- `use_normalized` to toggle between Enrichment Score and Normalized Enrichment Score
- `percentile_threshold` to adjust the significance threshold
- `color` to change the color scheme of the plot
- `max_overlaps` to control the density of labels in the plot

### Output

The function returns a list containing:

1. `data`: A dataframe with combined results from both datasets
2. `plot`: The ggplot object for the comparison plot
3. `common_pathways`: A list of common pathways (positive, negative, mixed, and all). This may be used downstream for further visualisation methods
4. `gmt`: A list of GMT-formatted dataframes for different pathway categories.

**Note!**: Various packages might expect the gmt matrix to be in specific formats to work with.

For example, to feed the GMT matrix to the library `TcGSA` I had to convert the output gmt into a specific format. I used this script of mine to convert it:

```R
# First, convert your matrix to a data.frame
gmt_df <- as.data.frame(res$gmt$gmt_mix) # Taking the GMT dataframe out of the NESscatter results function

# Assuming your GMT-like dataframe is `gmt_df`
gmt_list <- gmt_df %>%
  # Gather all the gene columns into a long format, removing NAs
  pivot_longer(cols = starts_with("Gene"), names_to = NULL, values_to = "Gene", values_drop_na = TRUE) %>%
  # Group by ID and Description to form the list for genesets
  group_by(ID, Description) %>%
  summarise(genesets = list(Gene), .groups = 'drop') %>%
  ungroup() %>%
  # Convert the result into a list containing the three required components
  summarise(
    genesets = list(genesets),
    geneset.names = list(ID),
    geneset.descriptions = list(Description)
  ) %>%
  # Turn it into a named list that TcGSA can work with
  {
    list(
      genesets = .$genesets[[1]],
      geneset.names = .$geneset.names[[1]],
      geneset.descriptions = .$geneset.descriptions[[1]]
    )
  }
```

For the `GSVA` package I had to convert it once again using this little data tyding script below:

```R
gmt_df <- as.data.frame(res$gmt$gmt_mix)

# Ensure the necessary libraries are loaded
require(dplyr)

# Assuming your GMT dataframe is called `gmt_df`
gmt_list <- gmt_df %>%
  # Pivot the gene columns into a long format
  pivot_longer(cols = starts_with("Gene"), names_to = NULL, values_to = "Gene", values_drop_na = TRUE) %>%
  # Group by the ID (gene set names)
  group_by(ID) %>%
  # Summarize into a list of gene identifiers for each gene set
  summarise(geneset = list(Gene)) %>%
  # Turn into a named list where the element names are the gene set names
  tibble::deframe()

# Now `gmt_list` is a list where each element is a vector of gene identifiers, and the names are the gene set names.
```


## Custom GSEA Dotplot Function

The R function `custom_dotplot` is a custom `ggplot2` GSEA dotplot. I was unsatisied with the standard `ClusterProfiler` functon `dotplot`, and so this function was born.
### Function Description

The `custom_dotplot` function generates a visually appealing and informative dotplot from GSEA results, allowing for various customizations including text processing, filtering, and sorting options.

#### Parameters

- `gsea_obj`: A GSEA result object
- `showCategory`: Number of categories to display (default: 10)
- `font.size`: Base font size for the plot (default: 10)
- `title`: Plot title (default: "GSEA Dotplot")
- `replace_`: Boolean to replace underscores with spaces in descriptions (default: TRUE)
- `capitalize_1`: Boolean to capitalize the first word in descriptions (default: TRUE)
- `capitalize_all`: Boolean to capitalize all words in descriptions (default: FALSE)
- `filterBy`: Method to filter results ("qvalue", "NES", "NES_positive", "NES_negative")
- `sortBy`: Method to sort results ("GeneRatio" or "qvalue")
- `q_cut`: Q-value cutoff for significance (default: 0.05)
- `min.dotSize`: Minimum dot size in the plot (default: 2)

#### Dependencies

The function requires the following R packages:

- ggplot2
- dplyr
- stringr
- scales

### Usage

```R
source("custom_dotplot.R")

# Assuming you have a GSEA result object named 'gsea_result'
plot <- custom_dotplot(gsea_result, showCategory = 15, filterBy = "NES", sortBy = "GeneRatio")
print(plot)
```

### Example
![alt_text](https://github.com/tony-zhelonkin/R_GSEA_visualisations/blob/main/examples/dotplot/1way.png?raw=true)


### Features

1. Flexible filtering options based on q-value, NES (Normalized Enrichment Score), or direction of enrichment
2. Customizable sorting of results by Gene Ratio or q-value
3. Text processing options for pathway descriptions (capitalization, underscore replacement)
4. Visually appealing dotplot with:
   - Dot size representing significance (-log10(q-value))
   - Color indicating direction of enrichment (positive or negative NES)
   - Customizable appearance (font size, title, etc.)

### Customization

You can customize the plot by adjusting various parameters:

- `showCategory` to control the number of pathways displayed
- `filterBy` and `sortBy` for different result filtering and sorting strategies
- `q_cut` to set the significance threshold
- `font.size` and `title` for visual customization

### Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes.




## Custom GSEA Dotplot Comparison Function

This repository contains an R function `custom_dotplot_comparison` for creating a comparative dotplot visualization of Gene Set Enrichment Analysis (GSEA) results from two samples.

### Function Description

The `custom_dotplot_comparison` function generates a side-by-side dotplot comparison of GSEA results from two samples, allowing for easy visualization of enriched pathways across both datasets.

#### Parameters

- `gsea_obj_x`, `gsea_obj_y`: GSEA result objects for comparison
- `pathway_ids`: Vector of pathway IDs to include in the plot
- `font.size`: Base font size for the plot (default: 10)
- `title`: Plot title (default: "Combined GSEA Dotplot")
- `replace_`: Boolean to replace underscores with spaces in descriptions (default: TRUE)
- `capitalize_1`: Boolean to capitalize the first word in descriptions (default: FALSE)
- `capitalize_all`: Boolean to capitalize all words in descriptions (default: FALSE)
- `min.dotSize`: Minimum dot size in the plot (default: 2)
- `sample_x_name`, `sample_y_name`: Names for the two samples (default: "Sample X" and "Sample Y")
- `sortBy`: Method to sort pathways ("importance_score", "qvalue", "x", or "y") (default: "importance_score")
**Note!** Choosing the sorting method might be tricky, so just take an look to see what works best for your specific use case. The x and y methods sort either for x or y samples. The `qvalue` sorting usually looks the most awkward. The `importance_score` sorting is a multiply of `qvalue` and `GeneRatio` and so results in a visual hybrid sorting.

#### Dependencies

The function requires the following R packages:

- ggplot2
- dplyr
- stringr
- tidyr

### Usage

```R
source("custom_dotplot_comparison.R")

# Assuming you have two GSEA result objects named 'gsea_result1' and 'gsea_result2'
# and a vector of pathway IDs you want to compare
plot <- custom_dotplot_comparison(
  gsea_obj_x = gsea_result1,
  gsea_obj_y = gsea_result2,
  pathway_ids = c("PATHWAY1", "PATHWAY2", "PATHWAY3"),
  sample_x_name = "Treatment",
  sample_y_name = "Control",
  sortBy = "importance_score"
)
print(plot)
```

### Example Usage

**Note!** That I\`m using the `common_mix` pathways from the NESscatter function as an input. Thus I define the pathways, which I want to visually compare between two GSEA runs. You may choose, whether you want to compare the `common_mix` (top 10 negative & top 10 positive pathways), or just negative or positive pathways only.

```R
# Example usage with your GSEA objects
plot <- custom_dotplot_comparison(
  gsea_obj_x = gsea_GOMF_th1,
  gsea_obj_y = gsea_GOMF_th17,
  pathway_ids = res$common_pathways$common_mix,
  font.size = 10,
  title = "Combined common GSEA pathways dotplot
  Top 10 common negative & Top 10 common positive pathways",
  sortBy = "y",
  replace_ = TRUE,
  capitalize_1 = FALSE,
  capitalize_all = FALSE,
  min.dotSize = 3,
  sample_x_name = "Th1 Sample",
  sample_y_name = "Th17 Sample"
)

# Display the plot
print(plot)
```

![alt_text](https://github.com/tony-zhelonkin/R_GSEA_visualisations/blob/main/examples/dotplot/2way.png?raw=true)

### Features

1. Generates a side-by-side dotplot comparison of GSEA results from two samples
2. Allows customization of pathway descriptions (replacing underscores, capitalization)
3. Provides flexible sorting options for pathways
4. Uses color to distinguish between positive and negative Normalized Enrichment Scores (NES)
5. Represents significance with dot size (-log10(q-value))
6. Handles missing pathways gracefully with warnings

### Customization

You can customize the plot by adjusting various parameters:

- `font.size` for changing the base font size
- `replace_`, `capitalize_1`, and `capitalize_all` for pathway description formatting
- `min.dotSize` to set the minimum dot size
- `sortBy` to change the ordering of pathways in the plot

### Output

The function returns a ggplot object representing the comparative dotplot. The plot includes:

- Two facets, one for each sample
- Dots representing pathways, with size indicating significance and color indicating NES direction
- Customizable pathway descriptions on the y-axis
- Gene Ratio on the x-axis




## Multi-Style Pathway Volcano Plot Function

R function `analyze_pathway_volcano` creates customized volcano plots that highlight genes from a specific GSEA analysis pathway within the context of differential expression results. The function offers multiple styling options to suit different visualization preferences. It`s usefeul if you need to show the \`behaviour\` of pathway genes against the universe of all the background genes.

### Function Description

The `analyze_pathway_volcano` function generates visually informative volcano plots that emphasize genes from a specified pathway, allowing for various customizations including significance thresholds, labeling methods, visual parameters, and plotting styles.

### Example

![alt text](https://github.com/tony-zhelonkin/R_GSEA_visualisations/blob/main/examples/analyzePathVolcano/ssDNA_DE_projection_clean.png?raw=true)



### Parameters

- `pathway_name`: Name of the pathway to highlight. Must be equal to the DESCRIPTION of the pathway in the GSEA results
- `gsea_results`: GSEA results object containing pathway information
- `de_results`: Differential expression results data frame
- `p_cutoff`: P-value cutoff for significance (default: 0.05)
- `fc_cutoff`: Fold change cutoff for significance (default: 2.0)
- `label_method`: Method to determine which genes to label ("default", "fc", "p", or "all")
- `max_overlaps`: Maximum number of label overlaps allowed (default: 100)
- `style`: Plotting style to use ("clean", "claude", or "gpt")

### Dependencies

The function requires the following R packages:
- ggplot2
- dplyr
- stringr
- ggrepel

### Usage

```R
source("analyze_pathway_volcano.R")
# Assuming you have GSEA results and DE results objects
plot <- analyze_pathway_volcano("Pathway Name", gsea_results, de_results,
                                p_cutoff = 0.05, fc_cutoff = 2.0,
                                label_method = "all", style = "clean")
print(plot)
```

### Example Usage

```R
p <- analyze_pathway_volcano(
  pathway_name = "single-stranded DNA binding",
  gsea_results = gsea_sum_GOMF_th1,
  de_results = limma_r_Th17.39_vs_Th1.39,
  p_cutoff = 0.05,
  fc_cutoff = 2.0,
  label_method = "p",  # or "fc"
  max_overlaps = 10,
  style = "claude"
)

print(p)
```


![alt text](https://github.com/tony-zhelonkin/R_GSEA_visualisations/blob/main/examples/analyzePathVolcano/ssDNA_DE_projection_claude.png?raw=true)


### Features

1. Highlights genes from a specific pathway within the context of background DE genes. Be carefull choosing your background genes, and DE results, so that the visuals would make biological sense.
2. Flexible labeling options for genes of interest
3. Customizable significance thresholds for p-value and fold change
4. Multiple styling options: "clean", "claude", "gpt"
5. Visually appealing volcano plots with:
   - Color-coded points for pathway genes and significant genes
   - Dashed lines indicating significance thresholds
   - Labels for selected genes using `ggrepel` for optimal placement
   - Customizable appearance (title, subtitle, axis labels)

### Customization

You can customize the plot by adjusting various parameters:
- `p_cutoff` and `fc_cutoff` to set different significance thresholds
- `label_method` to control which genes are labeled
- `max_overlaps` to adjust label density and prevent overcrowding
- `style` to choose between different visual presentations

### Plot Styles

1. Clean Style ("clean"):
   - Minimalist design with no legend
   - Pathway genes highlighted in distinct colors based on significance
   - Non-pathway genes in grey

2. Claude Style ("claude"):
   - Enhanced volcano-style legend
   - Colorblind-friendly palette
   - Legend positioned in the top-right corner
   - Generated by Claude 3.5 Sonnet

3. GPT Style ("gpt"):
   - Distinct colors for different gene statuses
   - Legend positioned near the top-right corner
   - Generated by ChatGPT-4o




### Contributing

Feel free to fork this repository and submit pull requests with improvements or bug fixes. We welcome contributions for new styling options or additional customization features.


### License

MIT License

Copyright (c) [2024] [Anton Zhelonkin]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
