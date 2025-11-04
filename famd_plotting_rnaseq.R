## ============================================================================
## FAMD Plotting Functions for Bulk RNA-seq Analysis
## ============================================================================
##
## General-purpose FAMD visualization functions adapted for RNA-seq context.
## Can be used with GSVA pathway scores or gene expression data.
##
## Main functions:
##   - famd_biplot_rnaseq(): Samples + centroids + variable arrows/points
##   - famd_var_biplot_rnaseq(): Variables-only plot
##
## Dependencies: ggplot2, ggrepel, ggnewscale, FactoMineR, factoextra, dplyr
## ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(FactoMineR)
  library(factoextra)
  library(dplyr)
})

## ============================================================================
## Helper functions
## ============================================================================

## Null coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

## Save plot wrapper
save_plot <- function(p, path, w = 12, h = 9, bg = "white", dpi = 300) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(path, plot = p, width = w, height = h, bg = bg, dpi = dpi)
}

## Extract label from metadata
## Looks for short_label, then label, then returns original ID
label_from_meta <- function(var_id_raw, meta) {
  if (is.null(meta)) return(as.character(var_id_raw))
  base <- as.character(var_id_raw)
  row <- meta[match(base, meta$param_id), , drop = FALSE]

  for (col in c("short_label", "label", "name")) {
    if (col %in% names(meta)) {
      val <- row[[col]][1]
      if (!is.null(val) && !is.na(val) && nzchar(trimws(as.character(val)))) {
        return(as.character(val))
      }
    }
  }
  base
}

## Extract process category from metadata
## Returns "other" if not found
process_from_meta <- function(var_id_raw, meta) {
  if (is.null(meta)) return("other")
  base <- as.character(var_id_raw)
  p <- tolower(trimws(as.character(meta$process[match(base, meta$param_id)])))
  if (length(p) == 0 || is.na(p)) return("other")

  ## Categorize based on common biological processes
  dplyr::case_when(
    grepl("glutamine|glutamate", p) ~ "glutamine_metabolism",
    grepl("glycolysis|glucose", p) ~ "glycolysis",
    grepl("tca|krebs|citric", p) ~ "oxidative_phosphorylation",
    grepl("amino.*acid|aa_met", p) ~ "amino_acid_metabolism",
    grepl("immune|inflammatory|ifn|lps|cytokine", p) ~ "immune_response",
    grepl("dendritic|dc_|activation", p) ~ "dc_activation",
    grepl("ifng|interferon.*gamma", p) ~ "ifng_response",
    grepl("lps|lipopolysaccharide", p) ~ "lps_response",
    TRUE ~ "other"
  )
}

## Convert to numeric matrix
.as_num_mat <- function(M) {
  if (is.null(M)) return(NULL)
  M <- as.matrix(M)
  suppressWarnings(storage.mode(M) <- "numeric")
  M
}

## Compute scaling factor for variable arrows/points
## Ensures variables fit within the individual space with margin
.compute_var_scale <- function(ind_xy, var_xy, margin = 0.95) {
  ind_xy <- .as_num_mat(ind_xy)
  var_xy <- .as_num_mat(var_xy)

  if (is.null(ind_xy) || is.null(var_xy) || nrow(var_xy) == 0) return(1)

  ## Range of individuals
  rx <- range(ind_xy[, 1], na.rm = TRUE)
  ry <- range(ind_xy[, 2], na.rm = TRUE)

  ## Max extent of variables
  vx <- max(abs(var_xy[, 1]), na.rm = TRUE)
  vy <- max(abs(var_xy[, 2]), na.rm = TRUE)

  if (!is.finite(vx) || !is.finite(vy) || vx == 0 || vy == 0) return(1)

  ## Scale so variables fit within individual space
  margin * min(diff(rx) / (2 * vx), diff(ry) / (2 * vy))
}

## Build map of categorical variable base names and levels
## For FAMD categorical variables
.build_quali_map <- function(res_famd) {
  fac_df <- res_famd$call$quali.sup$quali.sup
  if (is.null(fac_df)) return(NULL)

  sup_idx <- res_famd$call$sup.var %||% integer(0)
  sup_names <- if (length(sup_idx)) colnames(res_famd$call$X)[sup_idx] else character(0)

  ## Active categorical variables (not supplementary)
  act_fac <- setdiff(colnames(fac_df), sup_names)
  if (!length(act_fac)) return(NULL)

  ## Create base/level mapping
  base_rep <- rep(act_fac, times = sapply(fac_df[act_fac], nlevels))
  level_rep <- unlist(lapply(fac_df[act_fac], levels), use.names = FALSE)

  tibble::tibble(base = base_rep, level = level_rep)
}

## ============================================================================
## Main plotting functions
## ============================================================================

#' FAMD Biplot for RNA-seq Data
#'
#' Creates a biplot showing samples, group centroids, and variable vectors.
#' Works with GSVA pathway scores or gene expression data from FAMD.
#'
#' @param res_famd FactoMineR FAMD result object
#' @param X Original data frame used for FAMD (with factor columns)
#' @param meta Optional metadata data frame with columns: param_id, short_label, process
#' @param axes Numeric vector of length 2, which dimensions to plot (default c(1,2))
#' @param group_var Name of grouping variable in X (default "group")
#' @param show_ellipses Logical, draw confidence ellipses around groups? (default TRUE)
#' @param group_palette Named vector of colors for groups
#' @param process_palette Named vector of colors for biological processes
#' @param label_size Size of text labels (default 3)
#' @param top_vars List with contrib or cos2 threshold for filtering variables
#' @param quali_show For categorical variables: "both", "only1", "only0"
#' @param quali_label For categorical variables: "with_level", "base", "base_if1"
#' @param debug Logical, print debug info? (default FALSE)
#'
#' @return ggplot object
#' @export
famd_biplot_rnaseq <- function(
    res_famd,
    X,
    meta = NULL,
    axes = c(1, 2),
    group_var = "group",
    show_ellipses = TRUE,
    group_palette = NULL,
    process_palette = NULL,
    label_size = 3,
    top_vars = NULL,
    quali_show = c("both", "only1", "only0"),
    quali_label = c("with_level", "base", "base_if1"),
    debug = FALSE
) {

  quali_show <- match.arg(quali_show)
  quali_label <- match.arg(quali_label)

  ## Ensure grouping variable is a factor
  if (group_var %in% names(X)) {
    X[[group_var]] <- as.factor(X[[group_var]])
  } else {
    stop("group_var '", group_var, "' not found in X")
  }

  ## Extract individual coordinates
  ind_xy <- factoextra::get_famd_ind(res_famd)$coord[, axes, drop = FALSE]
  colnames(ind_xy) <- c("x", "y")

  ## Combine with grouping variable
  ind <- cbind(X[group_var], as.data.frame(ind_xy))
  names(ind)[1] <- "group"  # Standardize name

  ## Calculate group centroids
  centroids <- ind %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(cx = mean(x), cy = mean(y), .groups = "drop")

  ## Extract variable coordinates
  vq <- factoextra::get_famd_var(res_famd, "quanti.var")$coord
  vl <- factoextra::get_famd_var(res_famd, "quali.var")$coord

  ## Compute scaling factor
  S <- .compute_var_scale(
    ind_xy,
    rbind(
      vq[, axes, drop = FALSE] %||% NULL,
      vl[, axes, drop = FALSE] %||% NULL
    ),
    margin = 0.95
  )

  ## Process quantitative variables (arrows)
  q_df <- NULL
  if (!is.null(vq)) {
    q_df <- as.data.frame(vq[, axes, drop = FALSE] * S)
    names(q_df) <- c("x", "y")
    q_df$orig_id <- rownames(vq)
    q_df$name <- rownames(vq)
    q_df$label <- vapply(q_df$name, label_from_meta, "", meta = meta)
    q_df$process <- vapply(q_df$name, process_from_meta, "", meta = meta)
    q_df$geom <- "arrow"
    q_df$alpha <- 1
  }

  ## Process qualitative variables (points)
  l_df <- NULL
  if (!is.null(vl)) {
    qualimap <- .build_quali_map(res_famd)

    if (!is.null(qualimap) && nrow(qualimap) == nrow(vl)) {
      l_df <- as.data.frame(vl[, axes, drop = FALSE] * S)
      names(l_df) <- c("x", "y")
      l_df$orig_id <- rownames(vl)
      l_df$base <- qualimap$base
      l_df$level <- qualimap$level
      l_df$name <- paste0(l_df$base, "=", l_df$level)

      ## Filter by level
      is_pos <- tolower(as.character(l_df$level)) %in% c("1", "true", "yes", "present")
      if (quali_show == "only1") l_df <- l_df[is_pos, , drop = FALSE]
      if (quali_show == "only0") l_df <- l_df[!is_pos, , drop = FALSE]

      ## Set alpha for visualization
      if (nrow(l_df) > 0) {
        is_pos <- tolower(as.character(l_df$level)) %in% c("1", "true", "yes", "present")
        l_df$alpha <- if (quali_show == "both") ifelse(is_pos, 1, 0.3) else 1
      }

      ## Create labels
      base_lbl <- vapply(l_df$base, label_from_meta, "", meta = meta)
      if (quali_label == "with_level") {
        l_df$label <- paste0(base_lbl, " = ", l_df$level)
      } else if (quali_label == "base") {
        l_df$label <- base_lbl
      } else {  # base_if1
        l_df$label <- ifelse(is_pos, base_lbl, paste0(base_lbl, " = ", l_df$level))
      }

      l_df$process <- vapply(l_df$base, process_from_meta, "", meta = meta)
      l_df$geom <- "level"
    }
  }

  ## Filter variables by contribution or cos2
  if (!is.null(top_vars)) {
    if (!is.null(top_vars$cos2)) {
      thr <- top_vars$cos2
      if (!is.null(q_df)) {
        keep_q <- names(which(rowSums(factoextra::get_famd_var(res_famd, "quanti.var")$cos2[, axes, drop = FALSE]) >= thr))
        q_df <- q_df[q_df$orig_id %in% keep_q, , drop = FALSE]
      }
      if (!is.null(l_df)) {
        keep_l <- names(which(rowSums(factoextra::get_famd_var(res_famd, "quali.var")$cos2[, axes, drop = FALSE]) >= thr))
        l_df <- l_df[l_df$orig_id %in% keep_l, , drop = FALSE]
      }
    }

    if (!is.null(top_vars$contrib)) {
      k <- top_vars$contrib
      if (!is.null(q_df)) {
        qq <- factoextra::get_famd_var(res_famd, "quanti.var")$contrib
        keep <- rownames(qq)[head(order(rowSums(qq[, axes, drop = FALSE]), decreasing = TRUE), k)]
        q_df <- q_df[q_df$orig_id %in% keep, , drop = FALSE]
      }
      if (!is.null(l_df)) {
        ll <- factoextra::get_famd_var(res_famd, "quali.var")$contrib
        keep <- rownames(ll)[head(order(rowSums(ll[, axes, drop = FALSE]), decreasing = TRUE), k)]
        l_df <- l_df[l_df$orig_id %in% keep, , drop = FALSE]
      }
    }
  }

  ## Combine all variables
  vars_all <- dplyr::bind_rows(q_df %||% NULL, l_df %||% NULL)

  ## Filter groups for ellipses (need at least 3 samples)
  groups_ok <- ind %>%
    dplyr::count(group) %>%
    dplyr::filter(n >= 3) %>%
    dplyr::pull(group)

  ## Default palettes
  if (is.null(group_palette)) {
    group_palette <- scales::hue_pal()(length(unique(ind$group)))
    names(group_palette) <- unique(ind$group)
  }

  if (is.null(process_palette)) {
    process_palette <- c(
      glutamine_metabolism = "#E63946",
      glycolysis = "#457B9D",
      oxidative_phosphorylation = "#2A9D8F",
      amino_acid_metabolism = "#F4A261",
      immune_response = "#E76F51",
      dc_activation = "#264653",
      ifng_response = "#E9C46A",
      lps_response = "#F4A261",
      other = "grey50"
    )
  }

  ## Build plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3) +
    ggplot2::geom_point(data = ind, ggplot2::aes(x, y, color = group),
                       size = 2.4, alpha = 0.9) +
    {
      if (show_ellipses && length(groups_ok) > 0) {
        ggplot2::stat_ellipse(
          data = ind %>% dplyr::filter(group %in% groups_ok),
          ggplot2::aes(x, y, color = group),
          level = 0.9
        )
      }
    } +
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(cx, cy, fill = group),
      shape = 21, size = 10, color = "white", stroke = 1.1
    ) +
    ggplot2::geom_text(
      data = centroids,
      ggplot2::aes(cx, cy, label = group),
      color = "white", fontface = "bold", size = 4, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = group_palette, name = "Group") +
    ggplot2::scale_fill_manual(values = group_palette, guide = "none") +
    ggnewscale::new_scale_color() +
    {
      if (!is.null(q_df) && nrow(q_df) > 0) {
        ggplot2::geom_segment(
          data = q_df,
          ggplot2::aes(x = 0, y = 0, xend = x, yend = y, color = process),
          arrow = grid::arrow(length = grid::unit(0.18, "cm")),
          linewidth = 0.6, alpha = 0.95
        )
      }
    } +
    {
      if (!is.null(l_df) && nrow(l_df) > 0) {
        ggplot2::geom_point(
          data = l_df,
          ggplot2::aes(x, y, color = process, alpha = alpha),
          size = 2.3
        )
      }
    } +
    {
      if (!is.null(vars_all) && nrow(vars_all) > 0) {
        ggrepel::geom_text_repel(
          data = vars_all,
          ggplot2::aes(x, y, label = label, color = process, alpha = alpha),
          size = label_size, segment.size = 0.2,
          box.padding = 0.3, max.overlaps = 300,
          show.legend = FALSE
        )
      }
    } +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_color_manual(
      values = {
        levs <- unique((vars_all$process %||% factor(levels = "other")) %>% as.character())
        base <- process_palette
        miss <- setdiff(levs, names(base))
        if (length(miss)) base <- c(base, setNames(rep("grey50", length(miss)), miss))
        base
      },
      name = "Process"
    ) +
    ggplot2::labs(
      title = "FAMD: Samples, Centroids, and Variables",
      subtitle = paste0("Quantitative variables as arrows; categorical levels as points; ",
                       "scaled to sample coordinates"),
      x = paste0("Dim ", axes[1], " (", round(res_famd$eig[axes[1], 2], 1), "%)"),
      y = paste0("Dim ", axes[2], " (", round(res_famd$eig[axes[2], 2], 1), "%)")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )

  p
}

#' FAMD Variables-Only Biplot
#'
#' Shows only variables in FAMD space (no samples).
#' Uses res_famd$var$coord (unsigned coordinates).
#'
#' @param res_famd FactoMineR FAMD result object
#' @param meta Optional metadata with param_id, short_label, process columns
#' @param axes Which dimensions to plot (default c(1,2))
#' @param process_palette Named vector of colors for processes
#' @param label_size Text label size (default 3)
#' @param top_vars List with contrib or cos2 threshold for filtering
#' @param debug Print debug info? (default FALSE)
#'
#' @return ggplot object
#' @export
famd_var_biplot_rnaseq <- function(
    res_famd,
    meta = NULL,
    axes = c(1, 2),
    process_palette = NULL,
    label_size = 3,
    top_vars = NULL,
    debug = FALSE
) {

  ## Extract all variable coordinates
  V <- res_famd$var$coord[, axes, drop = FALSE]

  ## Remove supplementary variables if present
  exclude_vars <- c("id", "name", "Sample.ID", "group", "genotype", "treatment", "batch")
  V <- V[!(rownames(V) %in% exclude_vars), , drop = FALSE]

  ## Build data frame
  df <- as.data.frame(V)
  names(df) <- c("x", "y")
  df$name <- rownames(V)
  df$label <- vapply(df$name, label_from_meta, "", meta = meta)
  df$process <- vapply(df$name, process_from_meta, "", meta = meta)

  ## Filter by top variables
  if (!is.null(top_vars)) {
    if (!is.null(top_vars$cos2)) {
      cs <- rowSums(res_famd$var$cos2[, axes, drop = FALSE])
      keep <- names(cs)[cs >= top_vars$cos2]
      df <- df[df$name %in% keep, , drop = FALSE]
    }

    if (!is.null(top_vars$contrib)) {
      k <- top_vars$contrib
      ct <- rowSums(res_famd$var$contrib[, axes, drop = FALSE])
      keep <- names(sort(ct, decreasing = TRUE))[seq_len(min(k, length(ct)))]
      df <- df[df$name %in% keep, , drop = FALSE]
    }
  }

  ## Default palette
  if (is.null(process_palette)) {
    process_palette <- c(
      glutamine_metabolism = "#E63946",
      glycolysis = "#457B9D",
      oxidative_phosphorylation = "#2A9D8F",
      amino_acid_metabolism = "#F4A261",
      immune_response = "#E76F51",
      dc_activation = "#264653",
      ifng_response = "#E9C46A",
      lps_response = "#F4A261",
      other = "grey50"
    )
  }

  ## Build plot
  ggplot2::ggplot(df, ggplot2::aes(x, y, color = process)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, linewidth = 0.3) +
    ggplot2::geom_point(size = 2.3, alpha = 0.95) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = label_size,
      box.padding = 0.3,
      segment.size = 0.2,
      max.overlaps = 300,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = {
        levs <- unique(df$process %||% "other")
        base <- process_palette
        miss <- setdiff(levs, names(base))
        if (length(miss)) base <- c(base, setNames(rep("grey50", length(miss)), miss))
        base
      },
      name = "Process"
    ) +
    ggplot2::labs(
      title = "FAMD: Variable Space (Unsigned Coordinates)",
      subtitle = "Coordinates from res_famd$var$coord",
      x = paste0("Dim ", axes[1], " (", round(res_famd$eig[axes[1], 2], 1), "%)"),
      y = paste0("Dim ", axes[2], " (", round(res_famd$eig[axes[2], 2], 1), "%)")
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}
