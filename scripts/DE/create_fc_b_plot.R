#' Create a B-statistic vs Log2FC Plot for Differential Expression Results
#'
#' Creates a customizable B-statistic vs fold change plot highlighting significant genes.
#' The B-statistic is the log-odds of differential expression, with values > 0
#' indicating genes are more likely than not to be differentially expressed.
#'
#' Assumes `de_results` rownames are gene identifiers and contains at least
#' the columns `logFC` and `B`.
#'
#' @param de_results Data frame with DE results (rownames = genes, requires 'logFC', 'B').
#' @param fc_cutoff Numeric, absolute log2 fold change cutoff for significance (default: 2.0).
#' @param B_cutoff Numeric, B-statistic cutoff for significance (default: 0).
#' @param max.overlaps Integer, maximum number of label overlaps allowed (default: 10).
#' @param top_n Integer, number of top genes by B-statistic to label on each side (default: 5).
#' @param highlight_gene Character vector, specific gene IDs to highlight (default: NULL).
#' @param label_method Character, method to determine which genes to label:
#'        "top" (top_n genes by B-statistic in each direction),
#'        "sig" (genes significant by both B & fold change),
#'        "B" (genes significant by B-statistic only),
#'        "log2fc" (genes significant by fold change only),
#'        "none" (no gene labels).
#' @param x_breaks Numeric, interval between x-axis breaks (default: 1).
#' @param title Character, plot title (default: "B-statistic vs Log2FC").
#' @param color_palette Named character vector of four colors for point categories:
#'        NS, Log2FC only, B-statistic only, Both.
#' @param show_grid Logical, whether to display grid lines (default: FALSE).
#'
#' @return A ggplot2 object representing the B vs FC plot.
#' @export
#' @import ggplot2 dplyr ggrepel
#'
#' @examples
#' # Assuming de_results is a data frame with columns: logFC, B and gene names as rownames
#' plot_obj <- create_B_FC_plot(de_results, fc_cutoff = 1.5)
create_B_FC_plot <- function(
    de_results,
    fc_cutoff = 2,
    B_cutoff = 0,
    max.overlaps = 10,
    top_n = 5,
    highlight_gene = NULL,
    label_method = "top",
    x_breaks = 2,
    title = "B-statistic vs Log2FC",
    color_palette = c(
        "NS" = "#7F7F7F", # grey
        "Log2FC" = "#0173B2", # blue
        "B-statistic" = "#029E73", # green
        "B-statistic & Log2FC" = "#D55E00" # orange
    ),
    y_padding = 1, # Added padding above/below for y-axis
    show_grid = FALSE) {
    # ───────────────────────────────── helpers ──────────────────────────
    shade <- function(hex, factor = .6) {
        rgb <- grDevices::col2rgb(hex) / 255 * factor
        grDevices::rgb(
            pmax(pmin(rgb, 1), 0)[1],
            pmax(pmin(rgb, 1), 0)[2],
            pmax(pmin(rgb, 1), 0)[3]
        )
    }

    custom_minimal_theme_with_grid <- if (file.exists("scripts/custom_minimal_theme.R")) {
        source("scripts/custom_minimal_theme.R", local = TRUE)
        custom_minimal_theme_with_grid
    } else {
        function() ggplot2::theme_minimal()
    }

    # ──────────────────────── sanity checks ────────────────────────────
    stopifnot(all(c("logFC", "B") %in% colnames(de_results)))

    # ────────────────── 1. annotate significance ───────────────────────
    df <- dplyr::mutate(de_results,
        sig_fc = abs(logFC) >= fc_cutoff,
        sig_B = B > B_cutoff,
        cat = dplyr::case_when(
            sig_fc & sig_B ~ "B-statistic & Log2FC",
            sig_fc ~ "Log2FC",
            sig_B ~ "B-statistic",
            TRUE ~ "NS"
        )
    )

    # ────────────────── 2. label selection ─────────────────────────────
    # Helper to select top_n genes by B-statistic on each side
    get_top <- function(side) {
        if (side == "up") {
            df |>
                dplyr::filter(.data$logFC > 0) |>
                dplyr::arrange(dplyr::desc(.data$B)) |>
                dplyr::slice_head(n = top_n)
        } else {
            df |>
                dplyr::filter(.data$logFC < 0) |>
                dplyr::arrange(dplyr::desc(.data$B)) |>
                dplyr::slice_head(n = top_n)
        }
    }

    if (label_method == "top") {
        lab_df <- dplyr::bind_rows(get_top("up"), get_top("down"))
    } else if (label_method == "sig") {
        lab_df <- df[df$cat == "B-statistic & Log2FC", ]
    } else if (label_method == "B") {
        lab_df <- df[df$sig_B, ]
    } else if (label_method == "log2fc") {
        lab_df <- df[df$sig_fc, ]
    } else {
        lab_df <- df[0, ]
    }

    if (!is.null(highlight_gene)) {
        lab_df <- dplyr::bind_rows(
            lab_df,
            df[rownames(df) %in% highlight_gene, ]
        ) |>
            dplyr::distinct()
    }

    # ────────────────── 3. axis limits & colours ───────────────────────
    xmax <- ceiling(max(abs(df$logFC)) / x_breaks) * x_breaks
    # Calculate y-axis limits ensuring the threshold is visible
    ymin <- min(floor(min(df$B)), B_cutoff - y_padding)
    ymax <- max(ceiling(max(df$B)), B_cutoff + y_padding)
    # Ensure there's always space above the threshold line
    if (ymax <= B_cutoff) {
        ymax <- B_cutoff + y_padding
    }

    # If no genes above threshold, ensure proper spacing
    if (max(df$B) < B_cutoff) {
        # Add more space above the threshold when no genes cross it
        ymax <- B_cutoff + y_padding * 2
    }

    dark_pal <- vapply(color_palette, shade, character(1))

    # ────────────────── 4. build ggplot ────────────────────────────────
    # Build ggplot with modified axis limits
    g <- ggplot2::ggplot(df, ggplot2::aes(logFC, B, colour = cat)) +
        ggplot2::geom_point(size = 2, alpha = .65) +
        ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
        ggplot2::geom_hline(yintercept = B_cutoff, linetype = "dashed") +
        ggplot2::scale_colour_manual(
            name = NULL,
            values = color_palette,
            breaks = names(color_palette),
            labels = c(
                "B-statistic & Log2FC" = sprintf("B > %.1f & |log2FC| ≥ %.1f", B_cutoff, fc_cutoff),
                "Log2FC"               = sprintf("|log2FC| ≥ %.1f", fc_cutoff),
                "B-statistic"          = sprintf("B > %.1f", B_cutoff),
                "NS"                   = "NS"
            )
        ) +
        ggplot2::scale_x_continuous(
            breaks = seq(-xmax, xmax, by = x_breaks),
            limits = c(-xmax, xmax)
        ) +
        ggplot2::coord_cartesian(ylim = c(ymin, ymax)) +
        ggplot2::labs(
            x = "log2(FC)",
            y = "B statistic (log-odds of DE)",
            title = title,
            caption = sprintf(
                "Dashed lines: horiz. – B > %.1f; vert. – |log2FC| ≥ %.1f",
                B_cutoff, fc_cutoff
            )
        )

    # Add threshold annotation for clarity when it's above most genes
    if (max(df$B) < B_cutoff) {
        # Add a text label to clearly indicate the B threshold line
        g <- g + ggplot2::annotate(
            "text",
            x = xmax * 0.85,
            y = B_cutoff + y_padding * 0.5,
            label = sprintf("B = %.1f threshold", B_cutoff),
            color = "black",
            size = 3.5,
            fontface = "italic"
        )
    }

    # Apply theming
    g <- g + custom_minimal_theme_with_grid()

    if (!show_grid) {
        g <- g + ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()
        )
    }

    # ────────────────── 5. labels ───────────────────────────────────────
    if (nrow(lab_df)) {
        g <- g + ggrepel::geom_text_repel(
            data = lab_df,
            ggplot2::aes(label = rownames(lab_df)),
            colour = ifelse(rownames(lab_df) %in% highlight_gene, "black",
                dark_pal[lab_df$cat]
            ),
            fontface = ifelse(rownames(lab_df) %in% highlight_gene, "bold", "plain"),
            size = 3.5,
            box.padding = .4,
            point.padding = .3,
            max.overlaps = max.overlaps,
            min.segment.length = 0,
            show.legend = FALSE
        )
    }

    return(g)
}
