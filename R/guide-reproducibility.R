.standardize_matching_columns <- function(guide_data, replicate, match_columns) {
  out <- guide_data
  for (column in match_columns) {
    if (!is.numeric(out[[column]])) {
      stop("Every matching column must be numeric.", call. = FALSE)
    }
    standardized <- numeric(nrow(out))
    for (rep_value in unique(as.character(out[[replicate]]))) {
      idx <- which(as.character(out[[replicate]]) == rep_value)
      values <- out[[column]][idx]
      finite <- is.finite(values)
      if (!any(finite)) {
        stop("A matching column has no finite value within a replicate.", call. = FALSE)
      }
      values[!finite] <- stats::median(values[finite])
      scale_value <- stats::sd(values)
      standardized[idx] <- if (!is.finite(scale_value) || scale_value == 0) 0 else
        (values - mean(values)) / scale_value
    }
    out[[paste0(".match_", column)]] <- standardized
  }
  out
}

.effect_cosines <- function(first, second, value_names) {
  joined <- merge(first[c("gene", "pas_id", "effect")],
                  second[c("gene", "pas_id", "effect")],
                  by = c("gene", "pas_id"), suffixes = c("_first", "_second"))
  if (!nrow(joined)) return(data.frame())
  by_gene <- split(joined, joined$gene)
  result <- do.call(rbind, lapply(by_gene, function(x) {
    data.frame(gene = x$gene[[1]],
               value = .cosine(x$effect_first, x$effect_second),
               stringsAsFactors = FALSE)
  }))
  names(result)[2] <- value_names
  result
}

.restore_rng_state <- function(old_kind, had_seed, old_seed) {
  do.call(RNGkind, as.list(old_kind))
  .restore_random_seed(had_seed, old_seed)
}

#' Collect lossless long guide effects from a target batch
#'
#' @param fits A `repliapa_fit_list` from [repliapa_fit_targets()] or a named
#'   list of `repliapa_fit` objects.
#' @return A long data frame with target, replicate, guide, gene, PAS, and effect.
#' @export
repliapa_collect_guide_effects <- function(fits) {
  if (inherits(fits, "repliapa_fit_list")) fits <- fits$fits
  if (!is.list(fits) || !length(fits) ||
      !all(vapply(fits, inherits, logical(1), what = "repliapa_fit"))) {
    stop("`fits` must contain at least one repliapa_fit object.", call. = FALSE)
  }
  result <- do.call(rbind, lapply(fits, function(x) {
    value <- x$guide_effects
    if (is.null(value) || !nrow(value)) return(NULL)
    value$target <- x$target
    value[c("target", "replicate", "guide", "gene", "pas_id", "effect")]
  }))
  if (is.null(result) || !nrow(result)) {
    stop("No tested guide effects are available.", call. = FALSE)
  }
  rownames(result) <- NULL
  result
}

#' Evaluate same-target guide reproducibility against a matched null
#'
#' Same-target guide pairs are formed within each biological replicate. Each
#' guide is matched deterministically to a guide from another target using
#' standardized sample-level covariates; the two matched guides must themselves
#' represent different targets. Cells are never resampled or treated as units.
#'
#' @param guide_effects Long table with columns `target`, `replicate`, `guide`,
#'   `gene`, `pas_id`, and numeric `effect`.
#' @param guide_data One row per guide and biological replicate, including the
#'   target label and numeric matching covariates.
#' @param match_columns Numeric columns used for nearest-neighbour matching.
#' @param replicate,target,guide Column names in `guide_data`.
#' @param n_boot Number of guide-pair and target bootstrap resamples. Set to zero
#'   to omit confidence intervals.
#' @param seed Bootstrap seed. Matching itself is deterministic.
#' @return A `repliapa_guide_reproducibility` object with gene-level pairs,
#'   pair and target summaries, matching records, and bootstrap intervals.
#' @export
repliapa_guide_reproducibility <- function(guide_effects, guide_data, match_columns,
                                           replicate = "replicate", target = "target",
                                           guide = "guide", n_boot = 2000L,
                                           seed = 20260812L) {
  if (!is.data.frame(guide_effects) || !is.data.frame(guide_data)) {
    stop("`guide_effects` and `guide_data` must be data frames.", call. = FALSE)
  }
  .require_columns(guide_effects,
                   c("target", "replicate", "guide", "gene", "pas_id", "effect"),
                   "guide_effects")
  .require_columns(guide_data, c(replicate, target, guide, match_columns), "guide_data")
  if (!length(match_columns) || anyDuplicated(match_columns)) {
    stop("`match_columns` must name at least one unique numeric covariate.", call. = FALSE)
  }
  if (anyDuplicated(guide_data[c(replicate, guide)])) {
    stop("`guide_data` must have one row per guide and replicate.", call. = FALSE)
  }
  if (!is.numeric(guide_effects$effect) || any(!is.finite(guide_effects$effect))) {
    stop("Guide effects must be finite numbers.", call. = FALSE)
  }
  if (anyDuplicated(guide_effects[c("target", "replicate", "guide", "gene", "pas_id")])) {
    stop("Guide effect rows must be unique by target/replicate/guide/gene/PAS.", call. = FALSE)
  }
  n_boot <- .validate_scalar_integer(n_boot, "n_boot", minimum = 0L)
  seed <- .validate_scalar_integer(seed, "seed", minimum = 0L)

  matching_data <- .standardize_matching_columns(guide_data, replicate, match_columns)
  z_columns <- paste0(".match_", match_columns)
  matching_data$.replicate <- as.character(matching_data[[replicate]])
  matching_data$.target <- as.character(matching_data[[target]])
  matching_data$.guide <- as.character(matching_data[[guide]])
  matching_data <- matching_data[order(matching_data$.replicate, matching_data$.target,
                                       matching_data$.guide), , drop = FALSE]

  pair_groups <- split(matching_data, interaction(matching_data$.replicate,
                                                   matching_data$.target,
                                                   drop = TRUE, lex.order = TRUE))
  same_pairs <- do.call(rbind, lapply(pair_groups, function(x) {
    guides <- sort(unique(x$.guide))
    if (length(guides) < 2L) return(NULL)
    combinations <- utils::combn(guides, 2L)
    data.frame(replicate = x$.replicate[[1]], target = x$.target[[1]],
               guide_1 = combinations[1, ], guide_2 = combinations[2, ],
               stringsAsFactors = FALSE)
  }))
  if (is.null(same_pairs) || !nrow(same_pairs)) {
    stop("No within-replicate same-target guide pairs are available.", call. = FALSE)
  }
  same_pairs <- same_pairs[order(same_pairs$replicate, same_pairs$target,
                                 same_pairs$guide_1, same_pairs$guide_2), , drop = FALSE]

  nearest <- function(query_guide, query_replicate, original_target,
                      excluded_target = character()) {
    query <- matching_data[matching_data$.guide == query_guide &
                             matching_data$.replicate == query_replicate, , drop = FALSE]
    candidates <- matching_data[
      matching_data$.replicate == query_replicate &
        matching_data$.target != original_target &
        !matching_data$.target %in% excluded_target, , drop = FALSE
    ]
    if (nrow(query) != 1L || !nrow(candidates)) {
      stop("A guide pair has no valid matched different-target candidate.", call. = FALSE)
    }
    query_vector <- unlist(query[1, z_columns, drop = FALSE], use.names = FALSE)
    candidate_matrix <- as.matrix(candidates[z_columns])
    candidates$.distance <- sqrt(rowSums((candidate_matrix - matrix(
      query_vector, nrow(candidate_matrix), length(query_vector), byrow = TRUE
    ))^2))
    candidates <- candidates[order(candidates$.distance, candidates$.target,
                                   candidates$.guide), , drop = FALSE]
    candidates[1, c(".guide", ".target", ".distance"), drop = FALSE]
  }

  matching <- do.call(rbind, lapply(seq_len(nrow(same_pairs)), function(i) {
    row <- same_pairs[i, ]
    first <- nearest(row$guide_1, row$replicate, row$target)
    second <- nearest(row$guide_2, row$replicate, row$target, first$.target)
    data.frame(
      replicate = row$replicate, target = row$target,
      guide_1 = row$guide_1, guide_2 = row$guide_2,
      null_guide_1 = first$.guide, null_target_1 = first$.target,
      null_guide_2 = second$.guide, null_target_2 = second$.target,
      match_distance_1 = first$.distance, match_distance_2 = second$.distance,
      stringsAsFactors = FALSE
    )
  }))
  matching$pair_id <- sprintf("P%04d", seq_len(nrow(matching)))
  matching$pair_index <- seq_len(nrow(matching))

  effects_by_guide <- split(guide_effects, paste(guide_effects$replicate,
                                                 guide_effects$guide, sep = "\r"))
  pair_results <- lapply(seq_len(nrow(matching)), function(i) {
    row <- matching[i, ]
    get_effect <- function(guide_value) effects_by_guide[[paste(row$replicate,
                                                                 guide_value, sep = "\r")]]
    first <- get_effect(row$guide_1)
    second <- get_effect(row$guide_2)
    null_first <- get_effect(row$null_guide_1)
    null_second <- get_effect(row$null_guide_2)
    if (any(vapply(list(first, second, null_first, null_second), is.null, logical(1)))) {
      return(NULL)
    }
    same <- .effect_cosines(first, second, "same_target_similarity")
    null <- .effect_cosines(null_first, null_second,
                            "matched_different_target_similarity")
    result <- merge(same, null, by = "gene")
    if (!nrow(result)) return(NULL)
    result$similarity_difference <- result$same_target_similarity -
      result$matched_different_target_similarity
    cbind(row[rep(1L, nrow(result)), , drop = FALSE], result, row.names = NULL)
  })
  pair_results <- do.call(rbind, pair_results)
  if (is.null(pair_results) || !nrow(pair_results)) {
    stop("No matched guide effect vectors were available.", call. = FALSE)
  }
  pair_results <- pair_results[c(
    "pair_id", "replicate", "target", "guide_1", "guide_2",
    "null_guide_1", "null_target_1", "null_guide_2", "null_target_2",
    "match_distance_1", "match_distance_2", "gene",
    "same_target_similarity", "matched_different_target_similarity",
    "similarity_difference", "pair_index"
  )]

  pair_groups <- split(pair_results, pair_results$pair_id)
  pair_summary <- do.call(rbind, lapply(pair_groups, function(x) data.frame(
    pair_id = x$pair_id[[1]], replicate = x$replicate[[1]], target = x$target[[1]],
    guide_1 = x$guide_1[[1]], guide_2 = x$guide_2[[1]],
    gene_n = sum(is.finite(x$similarity_difference)),
    mean_same_similarity = mean(x$same_target_similarity, na.rm = TRUE),
    mean_matched_null_similarity = mean(x$matched_different_target_similarity, na.rm = TRUE),
    mean_difference = mean(x$similarity_difference, na.rm = TRUE),
    median_difference = stats::median(x$similarity_difference, na.rm = TRUE),
    positive_fraction = mean(x$similarity_difference > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )))
  pair_summary <- pair_summary[order(match(pair_summary$pair_id, matching$pair_id)), , drop = FALSE]
  target_groups <- split(pair_summary, interaction(pair_summary$target,
                                                    pair_summary$replicate,
                                                    drop = TRUE, lex.order = TRUE))
  target_summary <- do.call(rbind, lapply(target_groups, function(x) data.frame(
    target = x$target[[1]], replicate = x$replicate[[1]], pair_n = nrow(x),
    mean_same_similarity = mean(x$mean_same_similarity),
    mean_matched_null_similarity = mean(x$mean_matched_null_similarity),
    mean_difference = mean(x$mean_difference),
    positive_pair_fraction = mean(x$mean_difference > 0),
    stringsAsFactors = FALSE
  )))

  bootstrap <- data.frame()
  if (n_boot > 0L) {
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(.restore_rng_state(old_kind, had_seed, old_seed), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    pair_boot <- replicate(n_boot, {
      selected <- sample.int(nrow(pair_summary), nrow(pair_summary), replace = TRUE)
      mean(pair_summary$mean_difference[selected])
    })
    targets <- unique(pair_summary$target)
    target_boot <- replicate(n_boot, {
      selected <- sample(targets, length(targets), replace = TRUE)
      rows <- unlist(lapply(selected, function(x) which(pair_summary$target == x)),
                     use.names = FALSE)
      mean(pair_summary$mean_difference[rows])
    })
    bootstrap <- data.frame(
      resampling_unit = c("guide_pair", "target"),
      estimate = mean(pair_summary$mean_difference),
      CI_low = c(stats::quantile(pair_boot, 0.025), stats::quantile(target_boot, 0.025)),
      CI_high = c(stats::quantile(pair_boot, 0.975), stats::quantile(target_boot, 0.975)),
      resamples = n_boot, seed = seed, stringsAsFactors = FALSE
    )
  }
  output <- list(
    results = pair_results, pair_summary = pair_summary,
    target_summary = target_summary, matching = matching,
    bootstrap = bootstrap,
    settings = list(match_columns = match_columns, n_boot = n_boot, seed = seed)
  )
  class(output) <- "repliapa_guide_reproducibility"
  output
}

#' @export
print.repliapa_guide_reproducibility <- function(x, ...) {
  cat("RepliAPA matched guide reproducibility\n")
  cat("  guide pairs / targets: ", nrow(x$pair_summary), " / ",
      length(unique(x$pair_summary$target)), "\n", sep = "")
  cat("  mean same-target minus matched-null cosine: ",
      format(mean(x$pair_summary$mean_difference), digits = 4), "\n", sep = "")
  invisible(x)
}
