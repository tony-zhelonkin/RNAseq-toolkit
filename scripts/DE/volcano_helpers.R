# R_GSEA_visualisations/scripts/DE/volcano_helpers.R
# --------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)                 # for multi-panel layout
# The custom_minimal_theme.R is now sourced by the master script

# ─────────────────────────────────────────────────────────────────────────────
#  • generate a single volcano (rotated: LogFC on Y, –log10P on X) •
# ─────────────────────────────────────────────────────────────────────────────
create_vertical_volcano <- function(de_results,
                                    fc_cutoff   = 2,
                                    p_cutoff    = 0.05,
                                    title       = "",
                                    x_breaks    = 2,
                                    max.overlaps = 10,
                                    label_method = "sig",
                                    color_pal   = c(NS    = "grey80",
                                                    Log2FC = "#009E73",
                                                    pval   = "#56B4E9",
                                                    both   = "#E69F00"))
{
  stopifnot(all(c("logFC","P.Value") %in% colnames(de_results)))

  df <- de_results %>%
        mutate(sig_fc = abs(logFC) > fc_cutoff,
               sig_p  = P.Value   < p_cutoff,
               cat    = case_when(sig_fc & sig_p ~ "both",
                                  sig_fc          ~ "Log2FC",
                                  sig_p           ~ "pval",
                                  TRUE            ~ "NS"))

  ## ── which genes to label? ─────────────────────────────────────────────
  lab_df <- switch(label_method,
                   sig    = df[df$cat == "both", ],
                   p      = df[df$sig_p , ],
                   log2fc = df[df$sig_fc, ],
                   none   = NULL,
                   df[df$cat == "both", ])

  xmax <- ceiling(max(-log10(df$P.Value))/x_breaks)*x_breaks
  ymax <- ceiling(max(abs(df$logFC)))

  p <- ggplot(df, aes(x = -log10(P.Value), y = logFC, colour = cat)) +
       geom_point(size = 1.8, alpha = .7) +
       geom_vline(xintercept = -log10(p_cutoff), linetype = "dashed") +
       geom_hline(yintercept =  c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       scale_colour_manual(values = color_pal, name = NULL,
                           breaks = c("both","Log2FC","pval","NS"),
                           labels = c(
                             both   = sprintf("p < %.2g & |logFC| > %.1f",
                                              p_cutoff, fc_cutoff),
                             Log2FC = sprintf("|logFC| > %.1f", fc_cutoff),
                             pval   = sprintf("p < %.2g", p_cutoff),
                             NS     = "NS")) +
       scale_x_continuous(breaks = seq(0, xmax, by = x_breaks),
                          expand = expansion(mult = 0.02)) +       # <─ tiny gap
       coord_cartesian(ylim  = c(-ymax, ymax)) +
       labs(x = expression(-log[10]*"(p)"),
            y = "logFC",
            title = title) +
       custom_minimal_theme_with_grid() +
       theme(legend.position = "right") +
       guides(colour = guide_legend(override.aes = list(size = 4, shape = 16))) # ☆

  ##  text labels ---------------------------------------------------------
  if (!is.null(lab_df) && nrow(lab_df) > 0) {
      p <- p +
           ggrepel::geom_text_repel(
                data  = lab_df,
                aes(label = rownames(lab_df)),
                size  = 3.5,
                segment.color = "black",
                box.padding   = .4,
                point.padding = .3,
                max.overlaps  = max.overlaps,
                min.segment.length = 0,
                show.legend   = FALSE)    # ☆ prevents the “a” bullet
  }

  return(p)
}

# ─────────────────────────────────────────────────────────────────────────────
#  • arrange four volcanoes in one row •
# ─────────────────────────────────────────────────────────────────────────────
combine_volcano_row <- function(volcano_list,
                                labels       = names(volcano_list),
                                max_overlaps = 30)
{
  ## 1. common limits ------------------------------------------------------
  global_y <- max(vapply(volcano_list,
                         \(p) max(abs(ggplot_build(p)$data[[1]]$y)), 0))
  global_x <- max(vapply(volcano_list,
                         \(p) max(ggplot_build(p)$data[[1]]$x), 0))

  ## 2. apply limits, keep labels, add margin -----------------------------
  volcano_list <- lapply(seq_along(volcano_list), function(i) {
      volcano_list[[i]] +
        coord_cartesian(xlim = c(0, global_x),
                        ylim = c(-global_y, global_y),
                        clip = "off") +          # let labels spill over
        ggtitle(labels[i]) +
        theme(plot.margin = margin(5, 20, 5, 5)) # 2 mm extra on the right
  })

  wrap_plots(volcano_list, nrow = 1, guides = "collect") &
      theme(legend.position = "bottom")
}



# ──────────────────────────────────────────────────────────────────────────
#  • conventional volcano •
# ──────────────────────────────────────────────────────────────────────────
create_standard_volcano <- function(de_results,
                                    fc_cutoff   = 2,
                                    p_cutoff    = 0.05,
                                    title       = "",
                                    x_breaks    = 1,
                                    max.overlaps = 10,
                                    label_method = "sig",
                                    highlight_gene = NULL,
                                    color_pal   = c(NS           = "grey80",
                                                    Log2FC       = "#009E73",
                                                    pval         = "#56B4E9",
                                                    both         = "#E69F00")) {

  stopifnot(all(c("logFC", "P.Value") %in% colnames(de_results)))

  df <- de_results %>%
    mutate(sig_fc = abs(logFC) > fc_cutoff,
           sig_p  = P.Value     < p_cutoff,
           cat    = case_when(sig_fc & sig_p ~ "both",
                              sig_fc          ~ "Log2FC",
                              sig_p           ~ "pval",
                              TRUE            ~ "NS"))

  # ── choose which genes to label ─────────────────────────────────────────
  lab_df <- switch(label_method,
                   sig    = df[df$cat == "both", ],
                   p      = df[df$sig_p, ],
                   log2fc = df[df$sig_fc, ],
                   none   = NULL,
                   df[df$cat == "both", ])   # default to sig

  if (!is.null(highlight_gene))
      lab_df <- unique(rbind(lab_df,
                             df[rownames(df) %in% highlight_gene, ]))

  # ── figure limits ───────────────────────────────────────────────────────
  xmax <- ceiling(max(abs(df$logFC))/x_breaks)*x_breaks
  ymax <- ceiling(max(-log10(df$P.Value)))

  # ── plot ────────────────────────────────────────────────────────────────
  g <- ggplot(df, aes(x = logFC, y = -log10(P.Value), colour = cat)) +
       geom_point(size = 1.8, alpha = .7) +
       geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
       geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
       scale_colour_manual(values = color_pal, name = NULL,
           breaks = c("both","Log2FC","pval","NS"),
           labels = c(
             both   = sprintf("p < %.2g  &  |logFC| > %.1f", p_cutoff, fc_cutoff),
             Log2FC = sprintf("|logFC| > %.1f", fc_cutoff),
             pval   = sprintf("p < %.2g", p_cutoff),
             NS     = "NS")) +
       scale_x_continuous(breaks = seq(-xmax, xmax, by = x_breaks),
                          limits = c(-xmax, xmax)) +
       coord_cartesian(ylim = c(0, ymax)) +
       labs(x = "logFC", y = expression(-log[10]*"(p)"), title = title) +
       custom_minimal_theme_with_grid() +
       theme(legend.position = "right")

  # ── text-only labels with ggrepel ───────────────────────────────────────
  if (!is.null(lab_df) && nrow(lab_df) > 0) {
     g <- g +
          ggrepel::geom_text_repel(
               data  = lab_df,
               aes(label = rownames(lab_df)),
               size  = 3.5,
               fontface = ifelse(rownames(lab_df) %in% highlight_gene, "bold", "plain"),
               segment.color = "black",
               box.padding = 0.4,
               point.padding = 0.3,
               max.overlaps = max.overlaps,
               min.segment.length = 0)
  }

  return(g)
}
