#' Enhanced GSEA Dotplot with Continuous NES Gradient
#'
#' Creates a dotplot showing gene ratio vs pathway descriptions.
#' Separates "what to plot" (top N pathways by metric) from "what to highlight"
#' (significance). Top pathways are always visible; significant ones are
#' highlighted with a black border. Uses continuous NES gradient coloring
#' (Blue-White-Orange colorblind-safe palette).
#'
#' @param gsea_obj GSEA result object
#' @param filterBy Method to sort/filter results:
#'        - "p.adjust": Sort by significance (default)
#'        - "NES": Sort by absolute NES magnitude
#'        - "NES_positive": Keep only positive NES, sort by NES
#'        - "NES_negative": Keep only negative NES, sort by NES
#' @param sortBy Secondary sort for display ("GeneRatio" or "p.adjust")
#' @param showCategory Number of top pathways to show
#' @param padj_cutoff Significance threshold for highlighting (drawing black outline)
#' @param title Plot title
#' @param wrap_width Width for text wrapping
#' @param neg_color Color for negative NES (default: colorblind-safe blue #2166AC)
#' @param mid_color Color for zero NES (default: white #F7F7F7)
#' @param pos_color Color for positive NES (default: colorblind-safe orange #B35806)
#' @param nes_limits Numeric vector of length 2 for symmetric NES limits (auto if NULL)
#' @param min.dotSize Minimum dot size
#' @param max.dotSize Maximum dot size
#' @param highlight_sig Whether to highlight significant points with black outline.
#'        Base points have no outline; legend shows solid black circles.
#' @param highlight_threshold FDR threshold for highlighting significant points.
#'        If NULL (default), uses padj_cutoff. Set explicitly to override.
#' @param strip_prefix Logical, whether to strip common prefixes like "HALLMARK_"
#' @param use_gradient Logical, use continuous gradient for NES (TRUE) or binary colors (FALSE)
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
    min.dotSize = 2,
    max.dotSize = 10,
    highlight_sig = TRUE,
    highlight_threshold = NULL,
    strip_prefix = TRUE,
    use_gradient = TRUE,
    nes_limits = NULL) {
    # ========================================================================
    # STAGE 1: DATA PREPARATION
    # Extract all pathways and compute display metrics (no filtering)
    # ========================================================================
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

    # ========================================================================
    # STAGE 2: PATHWAY SELECTION
    # Select which pathways to SHOW (not based on significance!)
    # - filterBy controls sorting method (NES, p.adjust, direction)
    # - showCategory controls how many pathways to display
    # Significance only controls HIGHLIGHTING (black outline), not visibility
    # ========================================================================
    gsea_data_filtered <- gsea_data
    if (filterBy == "NES_positive") {
        # Filter for positive NES values
        pos_data <- gsea_data_filtered[gsea_data_filtered$NES > 0, ]
        if (nrow(pos_data) > 0) {
            gsea_data_filtered <- pos_data[order(pos_data$NES, decreasing = TRUE), ]
        } else {
            return(ggplot2::ggplot() +
                ggplot2::labs(title = paste(title, "(No positive NES pathways)")))
        }
    } else if (filterBy == "NES_negative") {
        # Filter for negative NES values
        neg_data <- gsea_data_filtered[gsea_data_filtered$NES < 0, ]
        if (nrow(neg_data) > 0) {
            gsea_data_filtered <- neg_data[order(neg_data$NES), ]
        } else {
            return(ggplot2::ggplot() +
                ggplot2::labs(title = paste(title, "(No negative NES pathways)")))
        }
    } else if (filterBy == "NES") {
        # Sort by absolute NES magnitude
        gsea_data_filtered <- gsea_data_filtered[order(abs(gsea_data_filtered$NES), decreasing = TRUE), ]
    } else {
        # Default to p.adjust sorting (most significant first)
        gsea_data_filtered <- gsea_data_filtered[order(gsea_data_filtered[[sig_col]]), ]
    }

    if (nrow(gsea_data_filtered) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::labs(title = paste(title, "(No matching pathways)")))
    }

    # Take top categories
    gsea_data_filtered <- utils::head(gsea_data_filtered, showCategory)

    # ========================================================================
    # STAGE 3: DISPLAY ORDERING
    # Reorder selected pathways for y-axis arrangement
    # ========================================================================
    if (sortBy == "GeneRatio") {
        plot_data <- gsea_data_filtered[order(gsea_data_filtered$GeneRatio, decreasing = TRUE), ]
    } else {
        plot_data <- gsea_data_filtered[order(gsea_data_filtered[[sig_col]]), ]
    }

    # Reorder factor levels for proper y-axis display
    plot_data$Description <- factor(plot_data$Description,
        levels = rev(plot_data$Description)
    )

    # Calculate symmetric NES limits if not provided
    if (is.null(nes_limits)) {
        nes_max <- max(abs(plot_data$NES), na.rm = TRUE)
        nes_limits <- c(-nes_max, nes_max)
    }

    # ========================================================================
    # STAGE 4: VISUALIZATION
    # Build ggplot with two layers:
    # - Base layer: ALL selected pathways, colored by NES, NO outline
    # - Overlay layer: ONLY significant pathways get black outline
    # ========================================================================
    if (use_gradient) {
        # Use continuous fill for NES gradient
        # Base points have NO outline; black outline added only for significant points
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = GeneRatio, y = Description)) +
            ggplot2::geom_point(
                ggplot2::aes(
                    size = negLogPval,
                    fill = NES
                ),
                shape = 21,  # Filled circle with border
                stroke = 0,  # No outline on base points
                color = "transparent"  # Transparent border (NA causes removal in ggplot2 4.0+)
            )
    } else {
        # Use binary colors for NES direction
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = GeneRatio, y = Description)) +
            ggplot2::geom_point(
                ggplot2::aes(
                    size = negLogPval,
                    color = NES_sign
                )
            )
    }

    # Add outline for highly significant points if requested
    if (highlight_sig) {
        # Use explicit threshold if provided, otherwise default to padj_cutoff
        if (!is.null(highlight_threshold)) {
            high_sig_threshold <- as.numeric(highlight_threshold)
        } else {
            padj_cutoff_num <- as.numeric(padj_cutoff)
            if (is.na(padj_cutoff_num)) {
                warning("padj_cutoff is not numeric, using default value of 0.05")
                padj_cutoff_num <- 0.05
            }
            high_sig_threshold <- padj_cutoff_num
        }
        highlight_data <- plot_data[plot_data[[sig_col]] < high_sig_threshold, ]

        if (nrow(highlight_data) > 0) {
            p <- p +
                ggplot2::geom_point(
                    data = highlight_data,
                    ggplot2::aes(size = negLogPval),
                    shape = 21, color = "black", fill = NA, stroke = 2  # Thick black outline for significant
                )
        }
    }

    # Adjust font size based on number of categories
    y_font_size <- ifelse(nrow(plot_data) > 20, 8, 9)

    # Add appropriate color/fill scale based on mode
    if (use_gradient) {
        p <- p +
            ggplot2::scale_fill_gradient2(
                low = neg_color,
                mid = mid_color,
                high = pos_color,
                midpoint = 0,
                name = "NES",
                limits = nes_limits,
                oob = scales::squish
            )
    } else {
        p <- p +
            ggplot2::scale_color_manual(
                name = "Direction",
                values = c("Positive NES" = pos_color, "Negative NES" = neg_color)
            )
    }

    # Add remaining scales and theme
    p <- p +
        ggplot2::scale_size_continuous(
            name = if ("qvalue" %in% colnames(gsea_data)) {
                bquote(-log[10](q - value))
            } else {
                bquote(-log[10](p - value))
            },
            range = c(min.dotSize, max.dotSize)
        ) +
        # Legend bubbles: solid black filled circles with no outline
        ggplot2::guides(
            size = ggplot2::guide_legend(
                override.aes = list(
                    shape = 16,      # Solid circle (no border)
                    fill = "black",
                    color = "black"
                )
            )
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
