#' Enhanced GSEA Dotplot with Continuous NES Gradient
#'
#' Creates a dotplot showing gene ratio vs pathway descriptions with improved filtering
#' and highlighting of significant results. Uses continuous NES gradient coloring
#' (Blue-White-Orange colorblind-safe palette).
#'
#' @param gsea_obj GSEA result object
#' @param filterBy Method to filter results ("NES_positive", "NES_negative", "p.adjust", "NES")
#' @param sortBy How to sort results ("GeneRatio" or "p.adjust")
#' @param showCategory Number of pathways to show
#' @param padj_cutoff Adjusted p-value cutoff
#' @param title Plot title
#' @param wrap_width Width for text wrapping
#' @param neg_color Color for negative NES (default: colorblind-safe blue #2166AC)
#' @param mid_color Color for zero NES (default: white #F7F7F7)
#' @param pos_color Color for positive NES (default: colorblind-safe orange #B35806)
#' @param nes_limits NES color scale limits (default: c(-3.5, 3.5))
#' @param min.dotSize Minimum dot size
#' @param max.dotSize Maximum dot size
#' @param highlight_sig Whether to highlight highly significant points with black outline
#' @param strip_prefix Logical, whether to strip common prefixes like "HALLMARK_"
#'
#' @return A ggplot2 object
#' @export
#'
#' @note Requires format_pathway_name() function to be available in environment.
#'       This is typically sourced by run_gsea_analysis() before calling this function.
#'
#' @note Color scheme updated 2025-12-02 to use continuous NES gradient
#'       (colorblind-safe Blue-White-Orange) matching Python publication figures.

gsea_dotplot <- function(
    gsea_obj,
    filterBy = "p.adjust",
    sortBy = "GeneRatio",
    showCategory = 10,
    padj_cutoff = 0.05,
    title = "GSEA Dotplot",
    wrap_width = 50,
    neg_color = "#2166AC",
    mid_color = "#F7F7F7",
    pos_color = "#B35806",
    nes_limits = c(-3.5, 3.5),
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    strip_prefix = TRUE) {
    # Extract and filter data
    gsea_data <- as.data.frame(gsea_obj@result)

    # Use qvalue if present, otherwise p.adjust
    sig_col <- if ("qvalue" %in% colnames(gsea_data)) "qvalue" else "p.adjust"

    # Calculate Gene Ratio
    gsea_data$count <- stringr::str_count(gsea_data$core_enrichment, "/") +
        ifelse(nchar(gsea_data$core_enrichment) > 0, 1, 0)
    gsea_data$GeneRatio <- gsea_data$count / gsea_data$setSize
    gsea_data$negLogPval <- -log10(gsea_data[[sig_col]])

    # Format pathway names using smart capitalization with biological exceptions
    gsea_data$Description <- format_pathway_name(
        gsea_data$Description,
        use_formatting = TRUE,
        strip_prefix = strip_prefix
    )

    # Apply custom text wrapping
    gsea_data$Description <- sapply(gsea_data$Description, function(txt) {
        words <- unlist(strsplit(txt, " "))
        if (length(words) <= 1) {
            return(txt)
        }

        total_chars <- nchar(txt)

        if (total_chars > wrap_width * 1.5) {
            # Very long text: split into three parts
            third_point <- ceiling(length(words) / 3)
            two_thirds <- third_point * 2

            part1 <- paste(words[1:third_point], collapse = " ")
            part2 <- paste(words[(third_point + 1):two_thirds], collapse = " ")
            part3 <- paste(words[(two_thirds + 1):length(words)], collapse = " ")

            return(paste(part1, part2, part3, sep = "\n"))
        } else if (total_chars > wrap_width) {
            # Moderately long text: split in half
            mid_point <- ceiling(length(words) / 2)

            part1 <- paste(words[1:mid_point], collapse = " ")
            part2 <- paste(words[(mid_point + 1):length(words)], collapse = " ")

            return(paste(part1, part2, sep = "\n"))
        }

        return(txt)
    }, USE.NAMES = FALSE)

    # Initial filtering by significance
    gsea_data_filtered <- gsea_data[gsea_data[[sig_col]] < padj_cutoff, ]

    if (nrow(gsea_data_filtered) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::labs(title = paste(title, "(No significant pathways)")))
    }

    # Apply direction filtering if needed
    if (filterBy == "NES_positive") {
        # Filter for positive NES values
        pos_data <- gsea_data_filtered[gsea_data_filtered$NES > 0, ]
        if (nrow(pos_data) > 0) {
            gsea_data_filtered <- pos_data[order(pos_data$NES, decreasing = TRUE), ]
        } else {
            # If no positive NES values, return empty plot with message
            return(ggplot2::ggplot() +
                ggplot2::labs(title = paste(title, "(No positive NES pathways)")))
        }
    } else if (filterBy == "NES_negative") {
        # Filter for negative NES values
        neg_data <- gsea_data_filtered[gsea_data_filtered$NES < 0, ]
        if (nrow(neg_data) > 0) {
            gsea_data_filtered <- neg_data[order(neg_data$NES), ]
        } else {
            # If no negative NES values, return empty plot with message
            return(ggplot2::ggplot() +
                ggplot2::labs(title = paste(title, "(No negative NES pathways)")))
        }
    } else if (filterBy == "NES") {
        gsea_data_filtered <- gsea_data_filtered[order(abs(gsea_data_filtered$NES), decreasing = TRUE), ]
    } else {
        # Default to p.adjust sorting
        gsea_data_filtered <- gsea_data_filtered[order(gsea_data_filtered[[sig_col]]), ]
    }

    if (nrow(gsea_data_filtered) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::labs(title = paste(title, "(No matching pathways)")))
    }

    # Take top categories
    gsea_data_filtered <- utils::head(gsea_data_filtered, showCategory)

    # Sort for display
    if (sortBy == "GeneRatio") {
        plot_data <- gsea_data_filtered[order(gsea_data_filtered$GeneRatio, decreasing = TRUE), ]
    } else {
        plot_data <- gsea_data_filtered[order(gsea_data_filtered[[sig_col]]), ]
    }

    # Reorder factor levels for proper y-axis display
    plot_data$Description <- factor(plot_data$Description,
        levels = rev(plot_data$Description)
    )

    # Create base plot with continuous NES gradient coloring
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = GeneRatio, y = Description)) +
        ggplot2::geom_point(
            ggplot2::aes(
                size = negLogPval,
                color = NES
            )
        )

    # Add outline for highly significant points if requested
    if (highlight_sig) {
        # Ensure padj_cutoff is numeric before division
        padj_cutoff_num <- as.numeric(padj_cutoff)
        if (is.na(padj_cutoff_num)) {
            warning("padj_cutoff is not numeric, using default value of 0.05")
            padj_cutoff_num <- 0.05
        }
        high_sig_threshold <- padj_cutoff_num / 10 # More stringent threshold for highlighting
        highlight_data <- plot_data[plot_data[[sig_col]] < high_sig_threshold, ]

        if (nrow(highlight_data) > 0) {
            p <- p +
                ggplot2::geom_point(
                    data = highlight_data,
                    ggplot2::aes(size = negLogPval),
                    shape = 21, color = "black", fill = NA, stroke = 1
                )
        }
    }

    # Adjust font size based on number of categories
    y_font_size <- ifelse(nrow(plot_data) > 20, 8, 9)

    # Complete the plot with continuous NES gradient scale
    p <- p +
        ggplot2::scale_color_gradient2(
            low = neg_color,
            mid = mid_color,
            high = pos_color,
            midpoint = 0,
            limits = nes_limits,
            oob = scales::squish,
            name = "NES"
        ) +
        ggplot2::scale_size_continuous(
            name = if ("qvalue" %in% colnames(gsea_data)) {
                bquote(-log[10](q - value))
            } else {
                bquote(-log[10](p - value))
            },
            range = c(min.dotSize, max.dotSize)
        ) +
        ggplot2::labs(
            title = title,
            x = "Gene Ratio",
            y = NULL
        ) +
        custom_minimal_theme_with_grid() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            legend.position = "right",
            plot.margin = ggplot2::margin(10, 10, 10, 10),
            axis.text.y = ggplot2::element_text(size = y_font_size)
        )

    return(p)
}
