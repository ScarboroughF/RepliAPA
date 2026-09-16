.complete_guides <- function(sample_data, replicate_col, guide_col, rows, replicates) {
  x <- sample_data[rows, c(replicate_col, guide_col), drop = FALSE]
  table_value <- table(x[[guide_col]], factor(x[[replicate_col]], levels = replicates))
  required <- table_value[, replicates, drop = FALSE]
  rownames(required)[rowSums(required == 1L) == length(replicates)]
}

.fit_gene <- function(gene, idx, counts, pas_data, context_data, target_guides,
                      ntc_guides, replicate_col, guide_col, gene_col,
                      pseudocount, max_combinations) {
  replicates <- sort(unique(as.character(context_data[[replicate_col]])))
  guide_pool <- c(target_guides, ntc_guides)
  target_n <- length(target_guides)
  permutation_n <- choose(length(guide_pool), target_n)
  total_count <- sum(counts[idx, , drop = FALSE])
  if (!is.finite(permutation_n) || permutation_n > max_combinations) {
    return(list(gene = data.frame(
      gene = gene, pas_n = length(idx), total_count = total_count,
      target_guide_n = target_n, ntc_guide_n = length(ntc_guides),
      permutation_n = permutation_n, statistic = NA_real_, p_value = NA_real_,
      effect_norm = NA_real_, coordinate_shift = NA_real_, tested = FALSE,
      filtered = FALSE, abstained = TRUE, reason = "TOO_MANY_ADMISSIBLE_PERMUTATIONS",
      stringsAsFactors = FALSE), pas = NULL, guide = NULL, heterogeneity = NULL))
  }

  mat <- as.matrix(counts[idx, rownames(context_data), drop = FALSE])
  k <- nrow(mat)
  proportions <- .composition(mat, pseudocount)
  transformed <- sqrt(proportions)
  combinations <- utils::combn(guide_pool, target_n, simplify = FALSE)
  keys <- vapply(combinations, function(x) paste(sort(x), collapse = "\r"), character(1))
  observed_key <- paste(sort(target_guides), collapse = "\r")
  observed_index <- match(observed_key, keys)
  selection <- vapply(combinations, function(x) as.numeric(guide_pool %in% x),
                      numeric(length(guide_pool)))

  flatten_by_guide <- function(values) {
    do.call(rbind, lapply(replicates, function(rep_value) {
      wanted <- paste(rep_value, guide_pool, sep = "\r")
      available <- paste(context_data[[replicate_col]], context_data[[guide_col]], sep = "\r")
      columns <- match(wanted, available)
      if (anyNA(columns)) stop("Incomplete guide-by-replicate context.", call. = FALSE)
      values[, columns, drop = FALSE]
    }))
  }

  h <- flatten_by_guide(transformed)
  selected_sum <- h %*% selection
  total_sum <- rowSums(h)
  deltas <- selected_sum / target_n -
    (total_sum - selected_sum) / length(ntc_guides)
  permutation_statistics <- colSums(deltas^2)
  statistic <- permutation_statistics[[observed_index]]
  p_value <- mean(permutation_statistics >= statistic - 1e-15)

  raw <- flatten_by_guide(proportions)
  raw_target <- raw %*% selection[, observed_index]
  raw_delta <- raw_target / target_n -
    (rowSums(raw) - raw_target) / length(ntc_guides)
  replicate_effect <- matrix(raw_delta, nrow = k, ncol = length(replicates))
  effect <- rowMeans(replicate_effect)

  coordinate_shift <- NA_real_
  if (all(c("position", "strand") %in% colnames(pas_data))) {
    strand <- unique(as.character(pas_data[["strand"]][idx]))
    if (length(strand) == 1L && strand %in% c("+", "-")) {
      coordinate <- pas_data[["position"]][idx] * if (strand == "-") -1 else 1
      coordinate_shift <- sum(effect * coordinate)
    }
  }

  pas_result <- data.frame(
    gene = gene, pas_id = rownames(pas_data)[idx], effect = effect,
    stringsAsFactors = FALSE
  )
  for (r in seq_along(replicates)) {
    pas_result[[paste0("effect_", replicates[[r]])]] <- replicate_effect[, r]
  }

  guide_rows <- lapply(replicates, function(rep_value) {
    rep_rows <- context_data[[replicate_col]] == rep_value
    rep_guides <- as.character(context_data[[guide_col]][rep_rows])
    rep_props <- proportions[, rep_rows, drop = FALSE]
    ntc_mean <- rowMeans(rep_props[, rep_guides %in% ntc_guides, drop = FALSE])
    do.call(rbind, lapply(target_guides, function(guide_value) {
      column <- match(guide_value, rep_guides)
      data.frame(gene = gene, pas_id = rownames(pas_data)[idx],
                 replicate = rep_value, guide = guide_value,
                 effect = rep_props[, column] - ntc_mean,
                 stringsAsFactors = FALSE)
    }))
  })
  guide_result <- do.call(rbind, guide_rows)

  heterogeneity <- do.call(rbind, lapply(replicates, function(rep_value) {
    x <- guide_result[guide_result$replicate == rep_value, , drop = FALSE]
    vectors <- split(x$effect, x$guide)
    pairs <- utils::combn(names(vectors), 2L, simplify = FALSE)
    similarities <- vapply(pairs, function(pair) .cosine(vectors[[pair[[1]]]],
                                                          vectors[[pair[[2]]]]), numeric(1))
    data.frame(gene = gene, replicate = rep_value, guide_n = length(vectors),
               mean_pairwise_cosine = mean(similarities, na.rm = TRUE),
               effect_dispersion = mean(vapply(seq_len(k), function(site) {
                 stats::var(vapply(vectors, `[[`, numeric(1), site), na.rm = TRUE)
               }, numeric(1)), na.rm = TRUE), stringsAsFactors = FALSE)
  }))

  list(
    gene = data.frame(
      gene = gene, pas_n = k, total_count = total_count,
      target_guide_n = target_n, ntc_guide_n = length(ntc_guides),
      permutation_n = permutation_n, statistic = statistic, p_value = p_value,
      effect_norm = sqrt(sum(effect^2)), coordinate_shift = coordinate_shift,
      tested = TRUE, filtered = FALSE, abstained = FALSE, reason = "NONE",
      stringsAsFactors = FALSE),
    pas = pas_result, guide = guide_result, heterogeneity = heterogeneity
  )
}

#' Fit replicate-aware all-PAS composition shifts
#'
#' Fits one perturbation target against non-targeting guides. Guide identities,
#' with all selected biological replicates attached, are the permutation units.
#' Cells are not accepted by this function; use [repliapa_pseudobulk()] first.
#'
#' @param counts PAS-by-guide-replicate pseudobulk matrix.
#' @param sample_data Sample metadata aligned to columns of `counts`.
#' @param pas_data PAS metadata aligned to rows of `counts`. It must contain a
#'   gene column and may contain `position` and `strand`.
#' @param target Target label to test.
#' @param gene,replicate,guide,target_col,ntc Column names.
#' @param replicates Optional biological replicates to include.
#' @param pseudocount PAS-count pseudocount used before composition and square root.
#' @param min_gene_total Minimum total PAS count in the tested context.
#' @param max_combinations Maximum number of exhaustively enumerated admissible
#'   guide-label permutations.
#'   Genes exceeding it abstain; labels are never sampled or split silently.
#' @return A `repliapa_fit` object with gene, PAS, guide, and heterogeneity tables.
#' @export
repliapa_fit <- function(counts, sample_data, pas_data, target,
                         gene = "gene", replicate = "replicate", guide = "guide",
                         target_col = "target", ntc = "ntc", replicates = NULL,
                         pseudocount = 0.5, min_gene_total = 20,
                         max_combinations = 1000000) {
  repliapa_validate_pseudobulk(counts, sample_data, replicate, guide, target_col, ntc)
  .aligned_data(pas_data, rownames(counts), "pas_data")
  .require_columns(pas_data, gene, "pas_data")
  gene_value <- as.character(pas_data[[gene]])
  if (anyNA(gene_value) || any(gene_value == "")) {
    stop("The PAS-to-gene map is incomplete.", call. = FALSE)
  }
  pseudocount <- .validate_scalar_number(
    pseudocount, "pseudocount", lower = 0, lower_inclusive = FALSE
  )
  min_gene_total <- .validate_scalar_number(
    min_gene_total, "min_gene_total", lower = 0
  )
  max_combinations <- .validate_scalar_integer(
    max_combinations, "max_combinations", minimum = 1L
  )
  target <- .validate_target_label(target)
  replicates <- .validate_replicates(sample_data, replicate, replicates)

  target_rows <- as.character(sample_data[[target_col]]) == target & !sample_data[[ntc]] &
    as.character(sample_data[[replicate]]) %in% replicates
  ntc_rows <- sample_data[[ntc]] & as.character(sample_data[[replicate]]) %in% replicates
  target_guides <- sort(.complete_guides(sample_data, replicate, guide, target_rows, replicates))
  ntc_guides <- sort(.complete_guides(sample_data, replicate, guide, ntc_rows, replicates))
  if (length(target_guides) < 2L) {
    stop("At least two target guides complete across selected replicates are required.", call. = FALSE)
  }
  if (length(ntc_guides) < 2L) {
    stop("At least two NTC guides complete across selected replicates are required.", call. = FALSE)
  }
  keep <- as.character(sample_data[[replicate]]) %in% replicates &
    as.character(sample_data[[guide]]) %in% c(target_guides, ntc_guides)
  context_data <- sample_data[keep, , drop = FALSE]
  context_counts <- counts[, rownames(context_data), drop = FALSE]
  context_data <- context_data[order(as.character(context_data[[replicate]]),
                                     match(as.character(context_data[[guide]]),
                                           c(target_guides, ntc_guides))), , drop = FALSE]
  context_counts <- context_counts[, rownames(context_data), drop = FALSE]

  gene_indices <- split(seq_along(gene_value), gene_value)
  fitted <- lapply(names(gene_indices), function(gene_name) {
    idx <- gene_indices[[gene_name]]
    if (length(idx) < 2L) {
      return(list(gene = data.frame(
        gene = gene_name, pas_n = length(idx), total_count = sum(context_counts[idx, , drop = FALSE]),
        target_guide_n = length(target_guides), ntc_guide_n = length(ntc_guides),
        permutation_n = NA_real_, statistic = NA_real_, p_value = NA_real_,
        effect_norm = NA_real_, coordinate_shift = NA_real_, tested = FALSE,
        filtered = TRUE, abstained = FALSE, reason = "SINGLE_PAS_GENE",
        stringsAsFactors = FALSE), pas = NULL, guide = NULL, heterogeneity = NULL))
    }
    if (sum(context_counts[idx, , drop = FALSE]) < min_gene_total) {
      return(list(gene = data.frame(
        gene = gene_name, pas_n = length(idx), total_count = sum(context_counts[idx, , drop = FALSE]),
        target_guide_n = length(target_guides), ntc_guide_n = length(ntc_guides),
        permutation_n = choose(length(target_guides) + length(ntc_guides), length(target_guides)),
        statistic = NA_real_, p_value = NA_real_, effect_norm = NA_real_,
        coordinate_shift = NA_real_, tested = FALSE, filtered = TRUE, abstained = FALSE,
        reason = "LOW_GENE_TOTAL", stringsAsFactors = FALSE),
        pas = NULL, guide = NULL, heterogeneity = NULL))
    }
    .fit_gene(gene_name, idx, context_counts, pas_data, context_data,
              target_guides, ntc_guides, replicate, guide, gene,
              pseudocount, max_combinations)
  })
  gene_table <- do.call(rbind, lapply(fitted, `[[`, "gene"))
  gene_table$q_value <- NA_real_
  gene_table$q_value[gene_table$tested] <- stats::p.adjust(gene_table$p_value[gene_table$tested], "BH")
  gene_table$significant <- gene_table$tested & gene_table$q_value <= 0.05
  output <- list(
    target = target, replicates = replicates,
    target_guides = target_guides, ntc_guides = ntc_guides,
    gene_table = gene_table,
    pas_effects = do.call(rbind, lapply(fitted, `[[`, "pas")),
    guide_effects = do.call(rbind, lapply(fitted, `[[`, "guide")),
    heterogeneity = do.call(rbind, lapply(fitted, `[[`, "heterogeneity")),
    settings = list(pseudocount = pseudocount, min_gene_total = min_gene_total,
                    max_combinations = max_combinations)
  )
  class(output) <- "repliapa_fit"
  output
}

#' @export
print.repliapa_fit <- function(x, ...) {
  cat("RepliAPA fit\n")
  cat("  target: ", x$target, "\n", sep = "")
  cat("  biological replicates: ", paste(x$replicates, collapse = ", "), "\n", sep = "")
  cat("  target / NTC guides: ", length(x$target_guides), " / ", length(x$ntc_guides), "\n", sep = "")
  cat("  genes tested / filtered / abstained: ", sum(x$gene_table$tested), " / ",
      sum(x$gene_table$filtered), " / ", sum(x$gene_table$abstained), "\n", sep = "")
  cat("  BH discoveries: ", sum(x$gene_table$significant, na.rm = TRUE), "\n", sep = "")
  invisible(x)
}
