# Package-local CoReSh chunk index. Keys are normalized species-directory
# paths, so a refcache symlink that moves to a new snapshot gets a new entry.
.coresh_chunk_cache <- new.env(parent = emptyenv())

#' Normalise a CoReSh species label
#'
#' @param species One of `"human"`, `"hsa"`, `"mouse"`, or `"mmu"`.
#' @return `"hsa"` or `"mmu"`.
#' @keywords internal
.coresh_species_code <- function(species) {
  if (!is.character(species) || length(species) != 1L || is.na(species) ||
      !nzchar(species)) {
    stop("`species` must be one of \"human\", \"hsa\", \"mouse\", or ",
         "\"mmu\".", call. = FALSE)
  }
  aliases <- c(human = "hsa", hsa = "hsa", mouse = "mmu", mmu = "mmu")
  if (!species %in% names(aliases)) {
    stop("`species` must be one of \"human\", \"hsa\", \"mouse\", or ",
         "\"mmu\"; got ", sQuote(species), ".", call. = FALSE)
  }
  unname(aliases[[species]])
}

#' Validate a logical CoReSh argument
#'
#' @param x Value to validate.
#' @param name Argument name, without backticks.
#' @return `x`, invisibly.
#' @keywords internal
.coresh_logical_scalar <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  invisible(x)
}

#' Validate a positive whole-number CoReSh argument
#'
#' @param x Value to validate.
#' @param name Argument name, without backticks.
#' @return The value as an integer.
#' @keywords internal
.coresh_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x)) {
    stop("`", name, "` must be a positive whole number.", call. = FALSE)
  }
  as.integer(x)
}

# `geseca()` takes no seed and draws its internal seeds from R's RNG, while
# `bplapply()` changes the generator, so pinning the seed alone is not enough.
.coresh_with_seed <- function(seed, expr) {
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv,
                         inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  on.exit({
    # A warning here describes the caller's legacy sampler, which we are
    # restoring, rather than a sampler selected by bulkiRNA.
    suppressWarnings(RNGkind(old_kind[1L], old_kind[2L], old_kind[3L]))
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  force(expr)
}

#' Resolve a CoReSh species directory
#'
#' @param chunk_dir Optional explicit species-directory path.
#' @param species Species alias accepted by [.coresh_species_code()].
#' @return A path carrying the attributes supplied by [.ref_path()].
#' @keywords internal
.coresh_chunk_dir <- function(chunk_dir, species) {
  code <- .coresh_species_code(species)
  if (!is.null(chunk_dir) &&
      (!is.character(chunk_dir) || length(chunk_dir) != 1L ||
       is.na(chunk_dir) || !nzchar(chunk_dir))) {
    stop("`chunk_dir` must be a non-empty character scalar or `NULL`.",
         call. = FALSE)
  }
  path <- if (is.null(chunk_dir)) {
    .ref_path("coresh", "preprocessed_chunks", code)
  } else {
    .ref_path("coresh", path = chunk_dir)
  }
  if (!dir.exists(path)) {
    stop("CoReSh `chunk_dir` is not a directory: ", path, ".",
         call. = FALSE)
  }
  path
}

#' List CoReSh chunk files
#'
#' @param chunk_dir A resolved species directory.
#' @return Sorted full paths to chunk files.
#' @keywords internal
.coresh_chunk_files <- function(chunk_dir) {
  paths <- sort(list.files(
    chunk_dir,
    pattern = "_full_objects\\.qs2$",
    full.names = TRUE
  ))
  if (!length(paths)) {
    stop(
      "No `*_full_objects.qs2` files found in ", chunk_dir,
      ". See `references/synapse-data-setup.md` for setup instructions.",
      call. = FALSE
    )
  }
  paths
}

#' Read one CoReSh chunk
#'
#' @param path Path to a qs2 chunk.
#' @return The deserialized list of dataset objects.
#' @keywords internal
.coresh_read_chunk <- function(path) {
  .require_pkg("qs2", "Reading CoReSh chunks",
               'install.packages("qs2")')
  qs2::qs_read(path)
}

#' Validate a CoReSh dataset object
#'
#' @param obj A dataset object from a CoReSh chunk.
#' @param context Description used in an error message.
#' @return `TRUE`, invisibly.
#' @keywords internal
.coresh_validate_object <- function(obj, context = "`obj`") {
  expected <- c("gseId", "gplId", "E1024", "rownames", "totalVar")
  if (!is.list(obj)) {
    stop(context, " must be a list with fields ",
         paste0("`", expected, "`", collapse = ", "), ".",
         call. = FALSE)
  }
  missing <- setdiff(expected, names(obj))
  if (length(missing)) {
    stop(context, " is missing field(s): ",
         paste0("`", missing, "`", collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.matrix(obj$E1024) || !is.numeric(obj$E1024)) {
    stop(context, "$E1024 must be a numeric matrix.", call. = FALSE)
  }
  if (!is.integer(obj$rownames) || length(obj$rownames) != nrow(obj$E1024)) {
    stop(context, "$rownames must be an integer Entrez vector with one value ",
         "per E1024 row.", call. = FALSE)
  }
  if (!is.numeric(obj$totalVar) || length(obj$totalVar) != 1L ||
      is.na(obj$totalVar) || !is.finite(obj$totalVar) || obj$totalVar <= 0) {
    stop(context, "$totalVar must be one finite positive number.",
         call. = FALSE)
  }
  for (field in c("gseId", "gplId")) {
    value <- obj[[field]]
    if (length(value) != 1L || !is.atomic(value)) {
      stop(context, "$", field, " must be a scalar value.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Build CoReSh reference provenance
#'
#' @param path A path returned by [.coresh_chunk_dir()].
#' @param species Normalized species code.
#' @param n_chunks Number of chunk files.
#' @return A named list.
#' @keywords internal
.coresh_provenance <- function(path, species, n_chunks) {
  list(
    source = attr(path, "source"),
    snapshot = attr(path, "snapshot"),
    path = unname(as.character(path)),
    species = species,
    n_chunks = as.integer(n_chunks)
  )
}

#' Empty CoReSh search result
#'
#' @return A zero-row tibble with the search schema.
#' @keywords internal
.coresh_empty_search <- function() {
  tibble::tibble(
    query_name = character(),
    gse = character(),
    gpl = character(),
    pct_var = numeric(),
    p_value = numeric(),
    log2err = numeric(),
    size = integer(),
    rank = integer()
  )
}

#' Index the CoReSh chunk tree
#'
#' Reads every chunk for one species and returns one row per GEO dataset,
#' including the chunk file containing it. The index is cached in the package
#' namespace by normalized directory path; set `cache = FALSE` to rebuild it.
#'
#' @param chunk_dir Optional explicit path to an `hsa` or `mmu` chunk directory.
#'   When `NULL`, resolve the corresponding directory in the shared refcache.
#' @param species One of `"human"`, `"hsa"`, `"mouse"`, or `"mmu"`.
#' @param cache Logical. Reuse an index built earlier in this R session.
#' @return A tibble with columns `gse`, `gpl`, and `chunk`. Its `provenance`
#'   attribute records the resolved reference snapshot.
#' @examples
#' \dontrun{
#' index <- coresh_chunks(species = "human")
#' }
#' @export
coresh_chunks <- function(chunk_dir = NULL, species = "human", cache = TRUE) {
  .coresh_logical_scalar(cache, "cache")
  code <- .coresh_species_code(species)
  resolved <- .coresh_chunk_dir(chunk_dir, species)
  normalized <- normalizePath(resolved, winslash = "/", mustWork = TRUE)

  if (cache && exists(normalized, envir = .coresh_chunk_cache,
                      inherits = FALSE)) {
    return(get(normalized, envir = .coresh_chunk_cache, inherits = FALSE))
  }

  paths <- .coresh_chunk_files(resolved)
  pieces <- lapply(paths, function(path) {
    chunk <- .coresh_read_chunk(path)
    if (!is.list(chunk) || !length(chunk)) {
      stop("CoReSh chunk did not contain a non-empty object list: ", path,
           ".", call. = FALSE)
    }
    lapply(seq_along(chunk), function(i) {
      obj <- chunk[[i]]
      .coresh_validate_object(
        obj,
        sprintf("Dataset %d in chunk %s", i, basename(path))
      )
      tibble::tibble(
        gse = as.character(obj$gseId),
        gpl = as.character(obj$gplId),
        chunk = path
      )
    }) |>
      dplyr::bind_rows()
  })
  out <- dplyr::bind_rows(pieces)
  attr(out, "provenance") <- .coresh_provenance(resolved, code, length(paths))
  assign(normalized, out, envir = .coresh_chunk_cache)
  out
}

#' Score one dataset for a CoReSh query
#'
#' Computes the percentage of stored total variance explained by the query's
#' shared profile and, optionally, its GESECA p-value.
#' Query IDs are treated as a set: duplicates are removed once during argument
#' validation, with a message. This deliberately differs from the upstream
#' vignette because counting one gene twice inflates both the score and set size
#' while GESECA tests unique matrix rows, making `pct_var` and `p_value` describe
#' different sets.
#'
#' @param obj One CoReSh dataset object with fields `gseId`, `gplId`, `E1024`,
#'   `rownames`, and `totalVar`.
#' @param query A non-empty integer vector of Entrez IDs. Duplicate IDs are
#'   removed with a message before scoring.
#' @param pvalues Logical. Calculate a GESECA p-value.
#' @param sample_size Positive whole number passed to the GESECA multilevel
#'   estimator as `sampleSize`.
#' @param seed Whole-number RNG seed reset immediately before the GESECA call.
#'   The generator is also pinned, so results do not depend on the number of
#'   cores or the caller's `RNGkind()`.
#' @param eps Positive numeric tolerance passed to the GESECA multilevel
#'   estimator.
#' @return A one-row tibble with `gse`, `gpl`, `pct_var`, `p_value`, `log2err`,
#'   and `size`. `log2err` is `Inf` for sets beyond the estimator's
#'   reliable resolution and is reported rather than hidden. Because the same
#'   seed is reset for every dataset, p-values are not independent across the
#'   compendium; multiple-testing correction across them is invalid and
#'   `p_value` is for ranking.
#' @examples
#' obj <- list(
#'   gseId = "GSE1", gplId = "GPL1",
#'   E1024 = matrix(c(1024L, 0L, 0L, 1024L), nrow = 2L),
#'   rownames = c(1L, 2L), totalVar = 2
#' )
#' coresh_match(obj, c(1L, 2L))
#' @export
coresh_match <- function(obj, query, pvalues = FALSE,
                         sample_size = 21L, seed = 1L, eps = 1e-300) {
  .coresh_logical_scalar(pvalues, "pvalues")

  .coresh_validate_object(obj)
  if (!is.integer(query) || !length(query) || anyNA(query)) {
    stop("`query` must be a non-empty integer vector of Entrez IDs.",
         call. = FALSE)
  }
  deduplicated <- unique(query)
  n_duplicates <- length(query) - length(deduplicated)
  if (n_duplicates) {
    query <- deduplicated
    message("Dropped ", n_duplicates, " duplicate Entrez ID",
            if (n_duplicates == 1L) "" else "s", " from `query`.")
  }
  .coresh_positive_integer(sample_size, "sample_size")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("`seed` must be one finite whole number.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) ||
      !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one finite positive number.", call. = FALSE)
  }

  query_idx <- match(query, obj$rownames)
  query_idx <- query_idx[!is.na(query_idx)]
  k <- length(query_idx)
  if (!k) {
    return(tibble::tibble(
      gse = as.character(obj$gseId),
      gpl = as.character(obj$gplId),
      pct_var = 0,
      p_value = NA_real_,
      log2err = NA_real_,
      size = 0L
    ))
  }

  E <- obj$E1024 / 1024
  profile <- colSums(E[query_idx, , drop = FALSE])
  query_var <- sum(profile^2)
  p_value <- NA_real_
  log2err <- NA_real_
  if (pvalues) {
    rn <- as.character(obj$rownames)
    missing <- is.na(rn)
    if (any(missing)) {
      rn[missing] <- paste0("unmapped_", seq_len(sum(missing)))
    }
    rownames(E) <- make.unique(rn)
    query_rows <- rownames(E)[query_idx]

    # Match upstream's literal seed for every dataset. Reusing its permutation
    # stream gives common-random-number precision for comparisons and, because
    # it is reset per call, makes chunk scheduling irrelevant. The resulting
    # cross-dataset p-values are correlated, as documented in the return value.
    geseca_result <- .coresh_with_seed(
      seed,
      fgsea::geseca(
        pathways = list(query = query_rows),
        E = E,
        minSize = 1L,
        maxSize = min(k, nrow(E) - 1L),
        center = FALSE,
        scale = FALSE,
        sampleSize = sample_size,
        eps = eps,
        nproc = 1L
      )
    )
    if (!is.data.frame(geseca_result) || nrow(geseca_result) != 1L) {
      result_description <- if (is.data.frame(geseca_result)) {
        paste(nrow(geseca_result), "rows")
      } else {
        "a non-tabular result"
      }
      stop(
        "GESECA returned ", result_description,
        " for dataset gse=", sQuote(as.character(obj$gseId)),
        ", gpl=", sQuote(as.character(obj$gplId)),
        " (k = ", k, ", nrow(E) = ", nrow(E), ").",
        call. = FALSE
      )
    }
    p_value <- as.numeric(geseca_result$pval[[1L]])
    log2err <- as.numeric(geseca_result$log2err[[1L]])
  }

  tibble::tibble(
    gse = as.character(obj$gseId),
    gpl = as.character(obj$gplId),
    pct_var = query_var / k / obj$totalVar * 100,
    p_value = p_value,
    log2err = log2err,
    size = as.integer(k)
  )
}

#' Score one CoReSh chunk file
#'
#' @param path Chunk file path.
#' @param queries Validated named list of Entrez vectors.
#' @param pvalues Logical. Calculate GESECA p-values.
#' @param sample_size Positive whole number passed to [coresh_match()].
#' @param seed Whole-number RNG seed passed to [coresh_match()].
#' @param eps Positive numeric tolerance passed to [coresh_match()].
#' @return A search-shaped tibble without ranks.
#' @keywords internal
.coresh_score_file <- function(path, queries, pvalues, sample_size, seed, eps) {
  chunk <- .coresh_read_chunk(path)
  if (!is.list(chunk) || !length(chunk)) {
    stop("CoReSh chunk did not contain a non-empty object list: ", path, ".",
         call. = FALSE)
  }
  rows <- lapply(names(queries), function(query_name) {
    scored <- lapply(
      chunk,
      coresh_match,
      query = queries[[query_name]],
      pvalues = pvalues,
      sample_size = sample_size,
      seed = seed,
      eps = eps
    ) |>
      dplyr::bind_rows()
    scored$query_name <- query_name
    scored[c(
      "query_name", "gse", "gpl", "pct_var", "p_value", "log2err", "size"
    )]
  })
  dplyr::bind_rows(rows)
}

#' Search the CoReSh compendium
#'
#' Scores each named Entrez query against every dataset for one species. With
#' more than one core, work is parallelized over chunk files so each worker
#' holds only one chunk at a time.
#' Query IDs are treated as sets: duplicates are removed once during argument
#' validation, with one message per search. This deliberately differs from the
#' upstream vignette because duplicate genes would inflate `pct_var` and `size`
#' while GESECA tests unique matrix rows, so the reported score and p-value
#' would describe different sets.
#'
#' @param queries A non-empty named list of integer Entrez vectors, each of
#'   length at least three after duplicate IDs are removed.
#' @param chunk_dir Optional explicit path to an `hsa` or `mmu` chunk directory.
#' @param species One of `"human"`, `"hsa"`, `"mouse"`, or `"mmu"`.
#' @param n_cores Positive whole number. `1` uses base R and does not require
#'   BiocParallel.
#' @param pvalues Logical. Calculate GESECA p-values and rank by ascending
#'   `p_value` instead of descending `pct_var`.
#' @param sample_size Positive whole number passed to the GESECA multilevel
#'   estimator as `sampleSize`.
#' @param seed Whole-number RNG seed reset immediately before every dataset's
#'   GESECA call. The generator is also pinned, so results do not depend on the
#'   number of cores or the caller's `RNGkind()`.
#' @param eps Positive numeric tolerance passed to the GESECA multilevel
#'   estimator.
#' @return A tibble with columns `query_name`, `gse`, `gpl`, `pct_var`,
#'   `p_value`, `log2err`, `size`, and `rank`, ordered within query by ascending
#'   `p_value` when requested and descending `pct_var` otherwise. `log2err` is
#'   `Inf` for sets beyond the estimator's reliable resolution and is reported
#'   rather than hidden. Because the same seed is reset for every dataset,
#'   p-values are not independent across the compendium; multiple-testing
#'   correction across them is invalid and `p_value` is for ranking. The
#'   `provenance` attribute records the reference snapshot.
#' @examples
#' \dontrun{
#' hits <- coresh_search(
#'   list(iron = c(7037L, 4891L, 55240L)),
#'   species = "human",
#'   n_cores = 1L
#' )
#' }
#' @export
coresh_search <- function(queries, chunk_dir = NULL, species = "human",
                          n_cores = 4L, pvalues = FALSE,
                          sample_size = 21L, seed = 1L, eps = 1e-300) {
  .coresh_logical_scalar(pvalues, "pvalues")

  if (!is.list(queries) || !length(queries) || is.null(names(queries))) {
    stop("`queries` must be a non-empty named list of integer Entrez vectors.",
         call. = FALSE)
  }
  query_names <- names(queries)
  if (anyNA(query_names) || any(!nzchar(query_names))) {
    stop("`queries` must have non-missing, non-empty names.", call. = FALSE)
  }
  if (anyDuplicated(query_names)) {
    duplicate <- unique(query_names[duplicated(query_names)])[[1L]]
    stop("`queries` names must be unique; duplicated query ",
         sQuote(duplicate), ".", call. = FALSE)
  }
  n_duplicates <- 0L
  for (query_name in query_names) {
    query <- queries[[query_name]]
    if (!is.integer(query) || length(query) < 3L || anyNA(query)) {
      stop("Query ", sQuote(query_name), " must be an integer Entrez vector ",
           "of length at least 3 with no missing values.", call. = FALSE)
    }
    deduplicated <- unique(query)
    n_duplicates <- n_duplicates + length(query) - length(deduplicated)
    queries[[query_name]] <- deduplicated
  }
  if (n_duplicates) {
    message("Dropped ", n_duplicates, " duplicate Entrez ID",
            if (n_duplicates == 1L) "" else "s", " from `queries`.")
  }
  for (query_name in query_names) {
    if (length(queries[[query_name]]) < 3L) {
      stop("Query ", sQuote(query_name), " must contain at least 3 unique ",
           "Entrez IDs after duplicates are removed.", call. = FALSE)
    }
  }
  n_cores <- .coresh_positive_integer(n_cores, "n_cores")
  sample_size <- .coresh_positive_integer(sample_size, "sample_size")
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed)) {
    stop("`seed` must be one finite whole number.", call. = FALSE)
  }
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) ||
      !is.finite(eps) || eps <= 0) {
    stop("`eps` must be one finite positive number.", call. = FALSE)
  }
  code <- .coresh_species_code(species)
  resolved <- .coresh_chunk_dir(chunk_dir, species)
  paths <- .coresh_chunk_files(resolved)

  if (n_cores == 1L) {
    pieces <- lapply(
      paths,
      .coresh_score_file,
      queries = queries,
      pvalues = pvalues,
      sample_size = sample_size,
      seed = seed,
      eps = eps
    )
  } else {
    .require_pkg(
      "BiocParallel",
      "Parallel CoReSh search",
      'BiocManager::install("BiocParallel")'
    )
    param <- BiocParallel::MulticoreParam(n_cores, progressbar = TRUE)
    pieces <- BiocParallel::bplapply(
      paths,
      .coresh_score_file,
      queries = queries,
      pvalues = pvalues,
      sample_size = sample_size,
      seed = seed,
      eps = eps,
      BPPARAM = param
    )
  }

  unranked <- dplyr::bind_rows(pieces)
  ranked <- lapply(query_names, function(query_name) {
    part <- unranked[unranked$query_name == query_name, , drop = FALSE]
    # Radix ordering makes character tie-breaks independent of LC_COLLATE.
    ordering <- if (pvalues) {
      order(part$p_value, -part$pct_var, part$gse, part$gpl,
            na.last = TRUE, method = "radix")
    } else {
      order(-part$pct_var, part$gse, part$gpl, na.last = TRUE,
            method = "radix")
    }
    part <- part[ordering, , drop = FALSE]
    part$rank <- seq_len(nrow(part))
    part
  }) |>
    dplyr::bind_rows()
  if (!nrow(ranked)) ranked <- .coresh_empty_search()
  ranked <- ranked[c(
    "query_name", "gse", "gpl", "pct_var", "p_value", "log2err", "size",
    "rank"
  )]
  attr(ranked, "provenance") <- .coresh_provenance(
    resolved, code, length(paths)
  )
  ranked
}

#' Find datasets supported by multiple CoReSh queries
#'
#' @param ranking A data frame returned by [coresh_search()].
#' @param top_n Positive whole number of top hits to consider per query.
#' @param min_queries Positive whole number of distinct agreeing queries.
#' @return A tibble with `gse`, `n_queries`, comma-joined `queries`,
#'   `best_rank`, and `mean_pct_var`, ordered by agreement and best rank.
#' @examples
#' ranking <- tibble::tibble(
#'   query_name = c("a", "b"), gse = c("GSE1", "GSE1"),
#'   pct_var = c(4, 6), rank = c(1L, 2L)
#' )
#' coresh_convergence(ranking)
#' @export
coresh_convergence <- function(ranking, top_n = 10L, min_queries = 2L) {
  required <- c("query_name", "gse", "pct_var", "rank")
  if (!is.data.frame(ranking)) {
    stop("`ranking` must be a data frame returned by `coresh_search()`.",
         call. = FALSE)
  }
  missing <- setdiff(required, names(ranking))
  if (length(missing)) {
    stop("`ranking` is missing column(s): ",
         paste0("`", missing, "`", collapse = ", "), ".",
         call. = FALSE)
  }
  top_n <- .coresh_positive_integer(top_n, "top_n")
  min_queries <- .coresh_positive_integer(min_queries, "min_queries")
  if (!is.character(ranking$query_name) || anyNA(ranking$query_name) ||
      any(!nzchar(ranking$query_name))) {
    stop("`ranking$query_name` must contain non-missing, non-empty strings.",
         call. = FALSE)
  }
  if (!is.character(ranking$gse) || anyNA(ranking$gse) ||
      any(!nzchar(ranking$gse))) {
    stop("`ranking$gse` must contain non-missing, non-empty strings.",
         call. = FALSE)
  }
  if (!is.numeric(ranking$pct_var) || anyNA(ranking$pct_var) ||
      any(!is.finite(ranking$pct_var))) {
    stop("`ranking$pct_var` must contain finite numeric values.",
         call. = FALSE)
  }
  if (!is.numeric(ranking$rank) || anyNA(ranking$rank) ||
      any(!is.finite(ranking$rank)) || any(ranking$rank < 1) ||
      any(ranking$rank != floor(ranking$rank))) {
    stop("`ranking$rank` must contain positive whole numbers.",
         call. = FALSE)
  }

  empty <- tibble::tibble(
    gse = character(),
    n_queries = integer(),
    queries = character(),
    best_rank = integer(),
    mean_pct_var = numeric()
  )
  top <- ranking[ranking$rank <= top_n, required, drop = FALSE]
  if (!nrow(top)) return(empty)

  # A GSE can have more than one platform row for one query. Keep that
  # query's best row so it contributes once to both the count and the mean.
  # Radix ordering makes platform tie-breaks independent of LC_COLLATE.
  top <- top[order(top$gse, top$query_name, top$rank, -top$pct_var,
                   method = "radix"),
             , drop = FALSE]
  top <- top[!duplicated(top[c("gse", "query_name")]), , drop = FALSE]
  groups <- split(seq_len(nrow(top)), top$gse)
  rows <- lapply(names(groups), function(gse) {
    part <- top[groups[[gse]], , drop = FALSE]
    query_names <- sort(unique(part$query_name), method = "radix")
    tibble::tibble(
      gse = gse,
      n_queries = as.integer(length(query_names)),
      queries = paste(query_names, collapse = ", "),
      best_rank = as.integer(min(part$rank)),
      mean_pct_var = mean(part$pct_var)
    )
  }) |>
    dplyr::bind_rows()
  rows <- rows[rows$n_queries >= min_queries, , drop = FALSE]
  if (!nrow(rows)) return(empty)
  rows[order(-rows$n_queries, rows$best_rank, rows$gse, method = "radix"),
       , drop = FALSE]
}

#' Check whether CoReSh search prerequisites are available
#'
#' This preflight reports every check instead of stopping at the first failed
#' one. In particular, it reports installation of the upstream `coresh`
#' package while making clear that version 0.1.0 contains no callable R
#' functions.
#'
#' @param chunk_dir Optional explicit path to an `hsa` or `mmu` chunk directory.
#' @param species One of `"human"`, `"hsa"`, `"mouse"`, or `"mmu"`.
#' @return Invisibly, a tibble with columns `check`, `ok`, and `detail`.
#' @examples
#' \dontrun{
#' coresh_validate(species = "human")
#' }
#' @export
coresh_validate <- function(chunk_dir = NULL, species = "human") {
  .coresh_species_code(species)
  if (!is.null(chunk_dir) &&
      (!is.character(chunk_dir) || length(chunk_dir) != 1L ||
       is.na(chunk_dir) || !nzchar(chunk_dir))) {
    stop("`chunk_dir` must be a non-empty character scalar or `NULL`.",
         call. = FALSE)
  }

  rows <- list()
  add_check <- function(check, ok, detail) {
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      check = check,
      ok = isTRUE(ok),
      detail = as.character(detail)
    )
  }

  packages <- c("qs2", "coresh", "BiocParallel", "org.Hs.eg.db",
                "org.Mm.eg.db")
  installed <- vapply(packages, requireNamespace, quietly = TRUE,
                      FUN.VALUE = logical(1L))
  fixes <- c(
    qs2 = 'Install with install.packages("qs2").',
    coresh = paste(
      "Not required: coresh 0.1.0 ships no R code or callable functions;",
      "CoReSh p-values use `fgsea::geseca()` directly."
    ),
    BiocParallel = 'Install with BiocManager::install("BiocParallel").',
    org.Hs.eg.db = 'Install with BiocManager::install("org.Hs.eg.db").',
    org.Mm.eg.db = 'Install with BiocManager::install("org.Mm.eg.db").'
  )
  for (pkg in packages) {
    if (pkg == "coresh" && installed[[pkg]]) {
      version <- tryCatch(
        as.character(utils::packageVersion(pkg)),
        error = function(e) "unknown"
      )
      detail <- paste0(
        "Installed (version ", version, "), but coresh 0.1.0 ships no R ",
        "code or callable functions and is not required; CoReSh p-values ",
        "use `fgsea::geseca()` directly."
      )
    } else if (installed[[pkg]]) {
      detail <- paste0("Installed (version ",
                       as.character(utils::packageVersion(pkg)), ").")
    } else {
      detail <- fixes[[pkg]]
    }
    add_check(paste0("package: ", pkg), installed[[pkg]], detail)
  }

  resolved <- tryCatch(
    .coresh_chunk_dir(chunk_dir, species),
    error = function(e) e
  )
  directory_ok <- !inherits(resolved, "error")
  add_check(
    "chunk directory",
    directory_ok,
    if (directory_ok) {
      paste0("Resolved to ", unname(as.character(resolved)), ".")
    } else {
      paste0(conditionMessage(resolved), " Fix the refcache or pass `chunk_dir=`.")
    }
  )

  paths <- if (directory_ok) {
    sort(list.files(resolved, pattern = "_full_objects\\.qs2$",
                    full.names = TRUE))
  } else {
    character()
  }
  files_ok <- length(paths) > 0L
  add_check(
    "chunk files",
    files_ok,
    if (files_ok) {
      paste(length(paths), "chunk file(s) found.")
    } else {
      paste(
        "No `*_full_objects.qs2` files found; see",
        "`references/synapse-data-setup.md` for setup instructions."
      )
    }
  )

  structure_ok <- FALSE
  structure_detail <- "Not checked because no chunk file is available."
  if (files_ok && !installed[["qs2"]]) {
    structure_detail <- 'Not checked; install qs2 with install.packages("qs2").'
  } else if (files_ok) {
    checked <- tryCatch({
      first <- .coresh_read_chunk(paths[[1L]])
      if (!is.list(first) || !length(first)) {
        stop("the first chunk is not a non-empty list", call. = FALSE)
      }
      .coresh_validate_object(first[[1L]], "The first dataset object")
      TRUE
    }, error = function(e) e)
    structure_ok <- isTRUE(checked)
    structure_detail <- if (structure_ok) {
      "The first dataset carries gseId, gplId, E1024, rownames, and totalVar."
    } else {
      paste0(conditionMessage(checked), " Rebuild or refresh the chunk snapshot.")
    }
  }
  add_check("first chunk structure", structure_ok, structure_detail)

  out <- dplyr::bind_rows(rows)
  print(out, n = Inf)
  invisible(out)
}
