#!/usr/bin/env Rscript
## ============================================================================
## test_gsea_running_sum_plot.R
##
## Regression + smoke tests for gsea_running_sum_plot()
##
## Covers:
##   1. Multi-pathway MSigDB-style (default palette)
##   2. Palette keying liveness: unnamed and named palettes both run
##   3. Single pathway
##   4. Custom DB (pre-populated @geneSets, simulating SynGO/MitoPathways)
##   5. Edge case: pathway ID not in results → warning, not crash
##   6. COLOUR MAPPING REGRESSION — asserts the BUILT label→colour pairing
##      (tests 1-5 are liveness only and cannot see a permuted mapping)
##
## Usage:
##   cd /path/to/RNAseq-toolkit
##   Rscript tests/test_gsea_running_sum_plot.R
## ============================================================================

suppressPackageStartupMessages({
  library(methods)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(patchwork)
})

source("scripts/custom_minimal_theme.R")
source("scripts/GSEA/GSEA_plotting/format_pathway_names.R")
source("scripts/GSEA/GSEA_plotting/gsea_running_sum_plot.R")

## ============================================================================
## Test helpers
## ============================================================================

test_passed <- function(name) cat(sprintf("  [PASS] %s\n", name))
test_failed <- function(name, reason) cat(sprintf("  [FAIL] %s: %s\n", name, reason))

output_dir <- "tests/output/gsea_running_sum"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## Mock gseaResult builder
##
## Creates a self-consistent S4 gseaResult with:
##   - @geneList : named numeric vector (ranked genes, descending)
##   - @geneSets : named list of character vectors (pathway -> gene names)
##   - @result   : data.frame with required GSEA columns
## ============================================================================

make_mock_gsea <- function(n_genes = 300, n_pathways = 20, seed = 42) {
  set.seed(seed)

  genes  <- paste0("Gene", seq_len(n_genes))
  ranks  <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  pathway_ids <- paste0("MOCK_PATH_", seq_len(n_pathways))
  gene_sets   <- lapply(pathway_ids, function(pid) {
    sample(genes, size = sample(20:60, 1))
  })
  names(gene_sets) <- pathway_ids

  nes_vals  <- rnorm(n_pathways, sd = 1.5)
  pvals     <- runif(n_pathways, 0.001, 0.5)
  padj_vals <- p.adjust(pvals, method = "BH")

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = paste("Mock Pathway", seq_len(n_pathways)),
    setSize        = vapply(gene_sets, length, integer(1)),
    NES            = nes_vals,
    pvalue         = pvals,
    p.adjust       = padj_vals,
    qvalue         = padj_vals * 1.05,
    rank           = sample(50:250, n_pathways),
    leading_edge   = rep("tags=30%, list=50%", n_pathways),
    core_enrichment = vapply(gene_sets, function(gs) {
      paste(gs[seq_len(min(5L, length(gs)))], collapse = "/")
    }, character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids   # enrichplot uses rowname lookup
  )

  new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Mus musculus",
    params   = list(exponent = 1)   # required by gseaScores()
  )
}

gsea_mock <- make_mock_gsea()

## ============================================================================
## Test 1: Multi-pathway MSigDB-style, default palette (palette = NULL)
## ============================================================================

cat("\n=== Test 1: Multi-pathway, default palette ===\n")

p1 <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1:5,
    palette      = NULL,
    title        = "Test 1: Default palette, 5 pathways"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p1) && inherits(p1, "patchwork")) {
  test_passed("Returns patchwork object with default palette")
  ggsave(file.path(output_dir, "test1_multi_pathway_default_pal.pdf"),
         p1, width = 10, height = 8)
  cat("  Saved: test1_multi_pathway_default_pal.pdf\n")
} else if (!is.null(p1) && inherits(p1, "gg")) {
  test_passed("Returns gg object with default palette")
} else {
  test_failed("Multi-pathway default palette", "did not return patchwork/gg object")
}

## ============================================================================
## Test 2: Palette keying liveness (unnamed AND named must both run)
##
## History: the Feb 2026 bug was naming the palette by PATHWAY ID, which never
## matches gseaplot2's mapped aesthetic (`Description`) and so silently produced
## grey50 keys for custom DBs. The remedy at the time — "never name the palette"
## — banned the only keying that works, leaving unnamed palettes matched
## POSITIONALLY against alphabetically sorted levels, i.e. silently permuted.
## Both keyings are now accepted and re-keyed onto the final plotted label.
##
## NOTE: this test is LIVENESS ONLY — it cannot detect a wrong colour mapping.
## The mapping itself is asserted in Test 6.
## ============================================================================

cat("\n=== Test 2: Palette keying liveness (unnamed + named) ===\n")

unnamed_pal <- c("#E41A1C", "#377EB8", "#4DAF4A")

p2a <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1:3,
    palette      = unnamed_pal,
    title        = "Test 2a: Unnamed palette (3 pathways)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR (unnamed): %s\n", e$message)); NULL
})

if (!is.null(p2a)) {
  test_passed("Unnamed palette runs without error")
  ggsave(file.path(output_dir, "test2a_unnamed_palette.pdf"),
         p2a, width = 10, height = 8)
} else {
  test_failed("Unnamed palette regression", "function errored with unnamed palette")
}

# Also verify the internal palette is not named after pathway IDs
# (we can only check indirectly by verifying the function runs with custom DB mock)
named_pal <- c(MOCK_PATH_1 = "#E41A1C", MOCK_PATH_2 = "#377EB8", MOCK_PATH_3 = "#4DAF4A")

p2b <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = c("MOCK_PATH_1", "MOCK_PATH_2", "MOCK_PATH_3"),
    palette      = named_pal,
    title        = "Test 2b: Named palette (should also work)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR (named): %s\n", e$message)); NULL
})

if (!is.null(p2b)) {
  test_passed("Named palette also runs without error")
} else {
  test_failed("Named palette", "function errored with named palette (known limitation)")
}

## ============================================================================
## Test 3: Single pathway
## ============================================================================

cat("\n=== Test 3: Single pathway ===\n")

p3 <- tryCatch({
  gsea_running_sum_plot(
    gsea_mock,
    gene_set_ids = 1L,
    title        = "Test 3: Single pathway"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p3)) {
  test_passed("Single pathway runs without error")
  ggsave(file.path(output_dir, "test3_single_pathway.pdf"),
         p3, width = 10, height = 6)
} else {
  test_failed("Single pathway", "function errored on single pathway")
}

## ============================================================================
## Test 4: Custom DB (pre-populated @geneSets, SynGO/MitoPathways pattern)
##
## For custom databases, @geneSets must be set manually after GSEA run.
## This tests the fix that makes running sum work for non-MSigDB databases.
## ============================================================================

cat("\n=== Test 4: Custom DB with pre-populated @geneSets ===\n")

# Simulate a custom DB result: IDs are short codes, Description = same as ID
make_custom_db_gsea <- function(seed = 99) {
  set.seed(seed)
  n_genes <- 200
  genes   <- paste0("Gene", seq_len(n_genes))
  ranks   <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  # Custom DB pathway IDs (not HALLMARK_ style)
  pathway_ids <- c("SYNGO:123", "SYNGO:456", "MITO:789", "MITO:012", "TRANSPORT:345")
  gene_sets   <- lapply(pathway_ids, function(pid) sample(genes, size = 25))
  names(gene_sets) <- pathway_ids

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = c("Synaptic vesicle exocytosis", "Postsynaptic density",
                       "Mitochondrial complex I", "ATP synthase",
                       "Solute carrier transport"),
    setSize        = 25L,
    NES            = c(2.1, -1.8, 1.5, -1.2, 1.9),
    pvalue         = c(0.001, 0.01, 0.03, 0.08, 0.002),
    p.adjust       = c(0.005, 0.02, 0.05, 0.12, 0.008),
    qvalue         = c(0.005, 0.02, 0.05, 0.12, 0.008),
    rank           = c(80L, 120L, 95L, 140L, 75L),
    leading_edge   = rep("tags=30%, list=50%", 5),
    core_enrichment = vapply(gene_sets, function(gs) paste(gs[1:5], collapse="/"), character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids
  )

  gsea_obj <- new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Mus musculus",
    params   = list(exponent = 1)
  )
  gsea_obj
}

gsea_custom <- make_custom_db_gsea()

p4 <- tryCatch({
  gsea_running_sum_plot(
    gsea_custom,
    gene_set_ids = c("SYNGO:123", "MITO:789", "TRANSPORT:345"),
    title        = "Test 4: Custom DB (SynGO/Mito pattern)"
  )
}, error = function(e) {
  cat(sprintf("  ERROR: %s\n", e$message)); NULL
})

if (!is.null(p4)) {
  test_passed("Custom DB with pre-populated @geneSets runs without error")
  ggsave(file.path(output_dir, "test4_custom_db.pdf"),
         p4, width = 10, height = 8)
} else {
  test_failed("Custom DB", "function errored with custom DB gene sets")
}

## ============================================================================
## Test 5: Edge case — pathway ID not in results → warning, not crash
## ============================================================================

cat("\n=== Test 5: Invalid pathway ID → warning, not crash ===\n")

# Mix of valid and invalid IDs
valid_id   <- gsea_mock@result$ID[1]
invalid_id <- "NONEXISTENT_PATH_XYZ"

p5 <- withCallingHandlers(
  tryCatch({
    gsea_running_sum_plot(
      gsea_mock,
      gene_set_ids = c(valid_id, invalid_id),
      title        = "Test 5: Invalid pathway ID (should warn)"
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message)); NULL
  }),
  warning = function(w) {
    cat(sprintf("  [OK] Warning issued: %s\n", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
)

if (!is.null(p5)) {
  test_passed("Invalid pathway ID: warning issued, function recovers with valid ID")
  ggsave(file.path(output_dir, "test5_invalid_id_recovery.pdf"),
         p5, width = 10, height = 6)
} else {
  test_failed("Invalid pathway ID", "function crashed instead of warning + recovering")
}

## ============================================================================
## All-invalid IDs → expect informative stop
## ============================================================================

cat("\n=== Test 5b: All-invalid IDs → informative stop ===\n")

all_invalid <- tryCatch({
  gsea_running_sum_plot(gsea_mock, gene_set_ids = c("BAD_1", "BAD_2"))
  "no_error"
}, error = function(e) {
  cat(sprintf("  [OK] Error caught as expected: %s\n", e$message))
  "error_caught"
})

if (all_invalid == "error_caught") {
  test_passed("All-invalid IDs: informative stop (not silent failure)")
} else {
  test_failed("All-invalid IDs", "function did not stop on all-invalid IDs")
}

## ============================================================================
## Test 6: COLOUR MAPPING REGRESSION — assert the BUILT label -> colour pairing
##
## The defect this guards: `palette` was forwarded UNNAMED to
## enrichplot::gseaplot2(color = ), which does scale_color_manual(values = color).
## With unnamed values ggplot matches BY POSITION against discrete levels it has
## sorted ALPHABETICALLY. So an "up"/"down" pair declared as c(up, down) rendered
## SWAPPED, because "... down" sorts before "... up". The name-keyed `labels` kept
## the legend TEXT correct, so it looked like a deliberate colour choice, not a bug.
##
## Every assertion below reads the mapping back out of the BUILT plot: the ES
## panel's rendered `colour` per `group`, with `group` resolved through the
## alphabetically sorted `Description` levels. That is the only way to see a
## permutation — "runs without error" cannot.
##
## The fixture is chosen so alphabetical order DIFFERS from declared order (an
## up/down pair), and all three accepted keyings are covered:
##   (a) unnamed, zipped to gene_set_ids in the caller's declared order
##   (b) named by pathway ID
##   (c) named by the final plotted label
## plus legend break order, the rug panel, and the restyle closure's overrides.
## ============================================================================

cat("\n=== Test 6: Built colour mapping (up/down permutation regression) ===\n")

## Read the ACTUAL label -> colour pairing out of a built plot.
built_color_map <- function(p) {
  b  <- ggplot2::ggplot_build(p[[1]])                     # ES panel
  d  <- unique(b$data[[1]][, c("colour", "group")])
  lv <- levels(factor(b$plot$data$Description))           # ggplot's sorted levels
  stats::setNames(toupper(d$colour), lv[d$group])
}

## Legend key order as the colour scale will draw it.
built_legend_order <- function(p) {
  b <- ggplot2::ggplot_build(p[[1]])
  as.character(b$plot$scales$get_scales("colour")$get_breaks())
}

## Fixture: an up/down pair whose alphabetical label order is the REVERSE of the
## declared order — exactly the case a positional mapping gets wrong.
make_updown_gsea <- function(seed = 7) {
  set.seed(seed)
  n_genes <- 300
  genes   <- paste0("Gene", seq_len(n_genes))
  ranks   <- sort(rnorm(n_genes, sd = 2), decreasing = TRUE)
  names(ranks) <- genes

  pathway_ids <- c("WT_heat_up", "WT_heat_down")          # DECLARED order
  gene_sets   <- lapply(pathway_ids, function(pid) sample(genes, size = 40))
  names(gene_sets) <- pathway_ids

  result_df <- data.frame(
    ID             = pathway_ids,
    Description    = pathway_ids,
    setSize        = 40L,
    NES            = c(2.1, -1.9),
    pvalue         = c(0.001, 0.002),
    p.adjust       = c(0.01, 0.02),
    qvalue         = c(0.01, 0.02),
    rank           = c(80L, 120L),
    leading_edge   = rep("tags=30%, list=50%", 2),
    core_enrichment = vapply(gene_sets, function(gs) paste(gs[1:5], collapse = "/"),
                             character(1)),
    stringsAsFactors = FALSE,
    row.names        = pathway_ids
  )

  new("gseaResult",
    result   = result_df,
    geneSets = gene_sets,
    geneList = ranks,
    readable = FALSE,
    keytype  = "SYMBOL",
    organism = "Homo sapiens",
    params   = list(exponent = 1)
  )
}

gsea_updown <- make_updown_gsea()
ud_ids      <- c("WT_heat_up", "WT_heat_down")            # DECLARED order
ud_labels   <- c(WT_heat_up = "WT_heat up", WT_heat_down = "WT_heat down")
BROWN       <- "#A6611A"
BLUE        <- "#2166AC"
## The one true mapping, whatever the keying: first declared id -> first colour.
ud_expected <- c("WT_heat up" = toupper(BROWN), "WT_heat down" = toupper(BLUE))

## Guard the fixture itself: if alphabetical == declared this test proves nothing.
if (!identical(sort(names(ud_expected)), names(ud_expected))) {
  test_passed("Fixture is discriminating (alphabetical label order != declared order)")
} else {
  test_failed("Fixture is discriminating",
              "alphabetical order equals declared order - test cannot detect a swap")
}

fmt_map <- function(m) paste(sprintf("%s=%s", names(m), m), collapse = ", ")

check_mapping <- function(name, p) {
  if (is.null(p)) {
    test_failed(name, "plot build returned NULL"); return(invisible(FALSE))
  }
  got <- built_color_map(p)
  if (!all(names(ud_expected) %in% names(got))) {
    test_failed(name, sprintf("expected labels %s, got %s",
                              paste(names(ud_expected), collapse = "/"),
                              paste(names(got), collapse = "/")))
    return(invisible(FALSE))
  }
  got <- got[names(ud_expected)]
  if (identical(unname(got), unname(ud_expected))) {
    test_passed(sprintf("%s: %s", name, fmt_map(got)))
    invisible(TRUE)
  } else {
    test_failed(name, sprintf("mapping is %s but must be %s (colours permuted)",
                              fmt_map(got), fmt_map(ud_expected)))
    invisible(FALSE)
  }
}

build_ud <- function(pal) {
  withCallingHandlers(
    tryCatch(
      gsea_running_sum_plot(gsea_updown, gene_set_ids = ud_ids,
                            labels = ud_labels, palette = pal),
      error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }),
    warning = function(w) {
      cat(sprintf("  [warn] %s\n", conditionMessage(w))); invokeRestart("muffleWarning")
    })
}

## (a) unnamed palette, zipped to gene_set_ids in the caller's declared order.
##     This is the case that used to render SWAPPED.
p6a <- build_ud(c(BROWN, BLUE))
check_mapping("Unnamed palette follows gene_set_ids declared order", p6a)

## (b) palette named by PATHWAY ID (the Feb-2026 trap: re-keyed, not banned).
p6b <- build_ud(c(WT_heat_up = BROWN, WT_heat_down = BLUE))
check_mapping("ID-named palette re-keyed to the plotted label", p6b)

## (c) palette named by the FINAL PLOTTED LABEL (used as-is).
p6c <- build_ud(c(`WT_heat up` = BROWN, `WT_heat down` = BLUE))
check_mapping("Label-named palette honoured as-is", p6c)

## All three keyings must be indistinguishable in the built plot.
if (!is.null(p6a) && !is.null(p6b) && !is.null(p6c) &&
    identical(built_color_map(p6a)[names(ud_expected)],
              built_color_map(p6b)[names(ud_expected)]) &&
    identical(built_color_map(p6a)[names(ud_expected)],
              built_color_map(p6c)[names(ud_expected)])) {
  test_passed("All three palette keyings yield the same built mapping")
} else {
  test_failed("Palette keying equivalence", "keyings disagree on the built mapping")
}

## Legend key order follows the DECLARED order, not the alphabetical one.
if (!is.null(p6a) && identical(built_legend_order(p6a), unname(ud_labels[ud_ids]))) {
  test_passed("Legend breaks follow declared order (breaks = names(palette))")
} else {
  test_failed("Legend break order",
              sprintf("got %s, expected %s",
                      paste(built_legend_order(p6a), collapse = " | "),
                      paste(unname(ud_labels[ud_ids]), collapse = " | ")))
}

## The rug panel must carry the SAME mapping as the ES curve it annotates.
if (!is.null(p6a)) {
  rb   <- ggplot2::ggplot_build(p6a[[2]])
  rd   <- unique(rb$data[[1]][, c("colour", "group")])
  rlv  <- levels(factor(rb$plot$data$Description))
  rmap <- stats::setNames(toupper(rd$colour), rlv[rd$group])[names(ud_expected)]
  if (identical(unname(rmap), unname(ud_expected))) {
    test_passed("Rug panel shares the ES panel's label -> colour mapping")
  } else {
    test_failed("Rug panel mapping", sprintf("got %s", fmt_map(rmap)))
  }
}

## The restyle closure must be able to override the palette (it previously could
## not: the palette was baked into the raw panels before the closure existed).
if (!is.null(p6c)) {
  p6d <- attr(p6c, "grs_restyle")(palette = c("#111111", "#EEEEEE"))
  got <- built_color_map(p6d)[names(ud_expected)]
  if (identical(unname(got), c("#111111", "#EEEEEE"))) {
    test_passed("restyle(palette = ) overrides colours without a panel rebuild")
  } else {
    test_failed("restyle palette override", sprintf("got %s", fmt_map(got)))
  }

  ## ...and overriding `labels` (a build-time knob) must rebuild AND re-key.
  p6e <- attr(p6c, "grs_restyle")(labels = c(WT_heat_up = "UP set",
                                             WT_heat_down = "DOWN set"),
                                  palette = c(BROWN, BLUE))
  got <- built_color_map(p6e)[c("UP set", "DOWN set")]
  if (identical(unname(got), c(toupper(BROWN), toupper(BLUE)))) {
    test_passed("restyle(labels = ) rebuilds panels and re-keys the palette")
  } else {
    test_failed("restyle labels override", sprintf("got %s", fmt_map(got)))
  }
}

## Single set keeps gseaplot2's plain black rug (guard against a visual regression).
p6f <- tryCatch(
  gsea_running_sum_plot(gsea_updown, gene_set_ids = "WT_heat_up", palette = BROWN),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
if (!is.null(p6f)) {
  es_col  <- toupper(unique(ggplot2::ggplot_build(p6f[[1]])$data[[1]]$colour))
  rug_col <- unique(ggplot2::ggplot_build(p6f[[2]])$data[[1]]$colour)
  if (identical(es_col, toupper(BROWN)) && identical(rug_col, "black")) {
    test_passed("Single set: ES uses the palette colour, rug stays black")
  } else {
    test_failed("Single set colours",
                sprintf("ES=%s rug=%s", paste(es_col, collapse = "/"),
                        paste(rug_col, collapse = "/")))
  }
}

## Default palette (NULL) must still be label-keyed, not positional.
p6g <- tryCatch(
  gsea_running_sum_plot(gsea_updown, gene_set_ids = ud_ids, labels = ud_labels,
                        palette = NULL),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
if (!is.null(p6g)) {
  got <- built_color_map(p6g)[names(ud_expected)]
  if (identical(unname(got), c("#E41A1C", "#377EB8"))) {
    test_passed("Default palette is label-keyed in declared order")
  } else {
    test_failed("Default palette keying", sprintf("got %s", fmt_map(got)))
  }
}

## A palette LONGER than the number of plotted sets must still be honoured.
## (gseaplot2 only applies `color=` when length(color) == length(geneSetID), so
## the old code silently fell back to ggplot's default hue scale here.)
p6h <- tryCatch(
  gsea_running_sum_plot(gsea_updown, gene_set_ids = ud_ids, labels = ud_labels,
                        palette = c(BROWN, BLUE, "#4DAF4A", "#984EA3")),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
check_mapping("Over-long palette truncated to declared order (not ignored)", p6h)

## Label-keyed palette when the label is long enough to be SOFT-WRAPPED: the
## plotted Description carries "\n", so keying must be whitespace-insensitive.
long_labels <- c(WT_heat_up   = "WT heat shock response upregulated core module",
                 WT_heat_down = "WT heat shock response downregulated core module")
p6i <- tryCatch(
  gsea_running_sum_plot(gsea_updown, gene_set_ids = ud_ids, labels = long_labels,
                        max_name_length = 20,
                        palette = stats::setNames(c(BROWN, BLUE), long_labels)),
  error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
if (!is.null(p6i)) {
  m6i <- built_color_map(p6i)
  wrapped <- vapply(long_labels, function(x) paste(strwrap(x, width = 20), collapse = "\n"),
                    character(1))
  got <- m6i[unname(wrapped)]
  if (identical(unname(got), c(toupper(BROWN), toupper(BLUE))) && any(grepl("\n", names(m6i)))) {
    test_passed("Label-keyed palette matches even after the label is soft-wrapped")
  } else {
    test_failed("Wrapped label keying",
                sprintf("got %s", fmt_map(m6i)))
  }
}

ggsave(file.path(output_dir, "test6_updown_mapping.pdf"), p6a, width = 10, height = 8)

## ============================================================================
## Summary
## ============================================================================

cat("\n=== All tests complete! ===\n")
cat("Output PDFs in:", output_dir, "\n\n")
cat("VISUAL INSPECTION CHECKLIST:\n")
cat("  test1: 5-panel running sum curves, legend shows 5 pathway names\n")
cat("  test2a: 3-panel running sum, clean unnamed colors\n")
cat("  test3: Single running sum curve\n")
cat("  test4: Custom DB paths (SynGO/Mito IDs), readable legend\n")
cat("  test5: 1-panel running sum (invalid ID removed with warning)\n")
cat("  test6: 'WT_heat up' curve is BROWN (#A6611A), 'WT_heat down' is BLUE (#2166AC),\n")
cat("         legend lists up before down (declared, not alphabetical, order)\n")
