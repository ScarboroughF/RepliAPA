.restore_random_seed <- function(had_seed, old_seed) {
  if (had_seed) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}

.generate_ntc_contrasts <- function(guides, guide_sizes, n_contrasts, seed) {
  if (!length(guide_sizes) || any(guide_sizes < 2L) || any(guide_sizes > length(guides) - 2L)) {
    stop("Every guide size must leave at least two pseudo-target and two control guides.",
         call. = FALSE)
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(.restore_random_seed(had_seed, old_seed), add = TRUE)
  set.seed(seed)
  contrast_seeds <- sample.int(.Machine$integer.max, n_contrasts)
  sizes <- rep(as.integer(guide_sizes), length.out = n_contrasts)
  selected <- lapply(seq_len(n_contrasts), function(i) {
    set.seed(contrast_seeds[[i]])
    sort(sample(guides, sizes[[i]], replace = FALSE))
  })
  data.frame(
    contrast_id = sprintf("NTC_%04d", seq_len(n_contrasts)),
    seed = contrast_seeds, target_matched_guide_n = sizes,
    pseudo_target_guides = vapply(selected, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

#' Calibrate the method with guide-unit NTC pseudo-target contrasts
#'
#' Non-targeting guide identities are relabeled as pseudo-targets; guide cells
#' are never split. The exact test in [repliapa_fit()] preserves all requested
#' biological replicate observations belonging to each guide.
#'
#' @inheritParams repliapa_fit
#' @param contrasts Optional predeclared data frame with `contrast_id` and a
#'   semicolon-delimited `pseudo_target_guides` column. Optional columns are
#'   `seed` and `target_matched_guide_n`.
#' @param guide_sizes Pseudo-target guide counts used when generating contrasts.
#' @param n_contrasts Number of generated contrasts.
#' @param seed Seed used only to generate and record contrast-specific seeds.
#' @return A `repliapa_null_calibration` object with contrasts, gene results,
#'   and calibration summaries.
#' @export
repliapa_ntc_calibration <- function(counts, sample_data, pas_data, contrasts = NULL,
                                     guide_sizes = c(2L, 3L, 4L), n_contrasts = 500L,
                                     seed = 1L, gene = "gene", replicate = "replicate",
                                     guide = "guide", target_col = "target", ntc = "ntc",
                                     replicates = NULL, pseudocount = 0.5,
                                     min_gene_total = 20, max_combinations = 1000000) {
  repliapa_validate_pseudobulk(counts, sample_data, replicate, guide, target_col, ntc)
  replicates <- .validate_replicates(sample_data, replicate, replicates)
  ntc_rows <- sample_data[[ntc]] & as.character(sample_data[[replicate]]) %in% replicates
  ntc_guides <- sort(.complete_guides(sample_data, replicate, guide, ntc_rows, replicates))
  if (length(ntc_guides) < 4L) stop("At least four complete NTC guides are required.", call. = FALSE)
  if (is.null(contrasts)) {
    n_contrasts <- .validate_scalar_integer(
      n_contrasts, "n_contrasts", minimum = 1L
    )
    seed <- .validate_scalar_integer(seed, "seed", minimum = 0L)
    if (!is.numeric(guide_sizes) || !length(guide_sizes) ||
        anyNA(guide_sizes) || any(!is.finite(guide_sizes)) ||
        any(abs(guide_sizes - round(guide_sizes)) > 1e-8) ||
        anyDuplicated(guide_sizes)) {
      stop("`guide_sizes` must contain unique finite integers.", call. = FALSE)
    }
    guide_sizes <- as.integer(guide_sizes)
    contrasts <- .generate_ntc_contrasts(ntc_guides, guide_sizes, n_contrasts, seed)
  } else {
    if (!is.data.frame(contrasts)) stop("`contrasts` must be a data.frame.", call. = FALSE)
    .require_columns(contrasts, c("contrast_id", "pseudo_target_guides"), "contrasts")
    contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
    contrast_ids <- as.character(contrasts$contrast_id)
    guide_strings <- as.character(contrasts$pseudo_target_guides)
    if (!nrow(contrasts) || anyNA(contrast_ids) || any(contrast_ids == "") ||
        anyDuplicated(contrast_ids) || anyNA(guide_strings) || any(guide_strings == "")) {
      stop("Contrast IDs must be unique and guide declarations complete.", call. = FALSE)
    }
    contrasts$contrast_id <- contrast_ids
    contrasts$pseudo_target_guides <- guide_strings
    if (!"seed" %in% names(contrasts)) contrasts$seed <- NA_integer_
    parsed_n <- lengths(strsplit(contrasts$pseudo_target_guides, ";", fixed = TRUE))
    if (!"target_matched_guide_n" %in% names(contrasts)) {
      contrasts$target_matched_guide_n <- parsed_n
    }
    declared_n <- contrasts$target_matched_guide_n
    if (!is.numeric(declared_n) || anyNA(declared_n) ||
        any(!is.finite(declared_n)) ||
        any(abs(declared_n - round(declared_n)) > 1e-8) ||
        any(parsed_n != declared_n)) {
      stop("Declared guide counts do not match pseudo-target guide strings.", call. = FALSE)
    }
  }
  parsed <- lapply(contrasts$pseudo_target_guides, function(x) sort(strsplit(x, ";", fixed = TRUE)[[1]]))
  valid <- vapply(parsed, function(x) length(x) >= 2L &&
                    length(setdiff(ntc_guides, x)) >= 2L && all(x %in% ntc_guides), logical(1))
  if (!all(valid)) stop("A contrast contains unavailable guides or inadequate groups.", call. = FALSE)
  contrasts$split_key <- vapply(parsed, paste, character(1), collapse = ";")

  ntc_data <- sample_data[ntc_rows & as.character(sample_data[[guide]]) %in% ntc_guides,
                          , drop = FALSE]
  ntc_counts <- counts[, rownames(ntc_data), drop = FALSE]
  unique_splits <- unique(contrasts[c("split_key", "target_matched_guide_n")])
  split_results <- lapply(seq_len(nrow(unique_splits)), function(i) {
    selected <- strsplit(unique_splits$split_key[[i]], ";", fixed = TRUE)[[1]]
    pseudo_data <- ntc_data
    is_selected <- as.character(pseudo_data[[guide]]) %in% selected
    pseudo_data[[target_col]] <- ifelse(is_selected, "PSEUDO_TARGET", "NT")
    pseudo_data[[ntc]] <- !is_selected
    fit <- repliapa_fit(ntc_counts, pseudo_data, pas_data, target = "PSEUDO_TARGET",
                        gene = gene, replicate = replicate, guide = guide,
                        target_col = target_col, ntc = ntc, replicates = replicates,
                        pseudocount = pseudocount, min_gene_total = min_gene_total,
                        max_combinations = max_combinations)
    result <- fit$gene_table
    result$split_key <- unique_splits$split_key[[i]]
    result
  })
  unique_results <- do.call(rbind, split_results)
  results <- merge(contrasts, unique_results, by = "split_key", all.x = TRUE,
                   sort = FALSE)
  result_order <- match(contrasts$contrast_id, results$contrast_id)
  if (anyNA(result_order)) stop("Failed to expand null contrasts.", call. = FALSE)
  results <- results[order(match(results$contrast_id, contrasts$contrast_id), results$gene), ]
  tested <- results$tested
  discoveries <- stats::aggregate(significant ~ contrast_id, results, sum)
  summary <- data.frame(
    contrasts = nrow(contrasts), tested_gene_results = sum(tested),
    empirical_FPR_0.01 = mean(results$p_value[tested] <= 0.01),
    empirical_FPR_0.05 = mean(results$p_value[tested] <= 0.05),
    empirical_FPR_0.10 = mean(results$p_value[tested] <= 0.10),
    mean_BH_false_discoveries_0.05 = mean(discoveries$significant),
    abstention_rate = mean(results$abstained), filtered_rate = mean(results$filtered),
    stringsAsFactors = FALSE
  )
  output <- list(contrasts = contrasts, results = results, summary = summary,
                 ntc_guides = ntc_guides, replicates = replicates)
  class(output) <- "repliapa_null_calibration"
  output
}

#' @export
print.repliapa_null_calibration <- function(x, ...) {
  cat("RepliAPA NTC null calibration\n")
  cat("  contrasts: ", x$summary$contrasts, "\n", sep = "")
  cat("  empirical FPR @ 0.05: ", format(x$summary$empirical_FPR_0.05, digits = 4), "\n", sep = "")
  cat("  mean BH false discoveries: ",
      format(x$summary$mean_BH_false_discoveries_0.05, digits = 4), "\n", sep = "")
  invisible(x)
}
