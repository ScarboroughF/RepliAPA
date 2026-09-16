#' Compare effects across independently fitted biological replicates
#'
#' The two inputs should be fitted separately, for example one development
#' replicate and one frozen validation replicate. No cells are split.
#'
#' @param development,validation Two `repliapa_fit` objects for the same target.
#' @param fdr Development/validation BH threshold.
#' @return A data frame with gene-wise effect-vector and hit replication metrics.
#'   A target-level summary is stored in the `summary` attribute.
#' @export
repliapa_reproducibility <- function(development, validation, fdr = 0.05) {
  if (!inherits(development, "repliapa_fit") || !inherits(validation, "repliapa_fit")) {
    stop("Both inputs must be `repliapa_fit` objects.", call. = FALSE)
  }
  if (!identical(development$target, validation$target)) {
    stop("Fits must represent the same target.", call. = FALSE)
  }
  overlap <- intersect(as.character(development$replicates),
                       as.character(validation$replicates))
  if (length(overlap)) {
    stop("Development and validation fits must use disjoint biological replicates.",
         call. = FALSE)
  }
  fdr <- .validate_scalar_number(
    fdr, "fdr", lower = 0, lower_inclusive = FALSE,
    upper = 1, upper_inclusive = FALSE
  )
  d <- development$pas_effects[c("gene", "pas_id", "effect")]
  v <- validation$pas_effects[c("gene", "pas_id", "effect")]
  names(d)[3] <- "effect_development"
  names(v)[3] <- "effect_validation"
  effects <- merge(d, v, by = c("gene", "pas_id"))
  by_gene <- split(effects, effects$gene)
  out <- do.call(rbind, lapply(by_gene, function(x) {
    pearson <- if (nrow(x) >= 3L && stats::sd(x$effect_development) > 0 &&
                   stats::sd(x$effect_validation) > 0) {
      stats::cor(x$effect_development, x$effect_validation)
    } else NA_real_
    data.frame(gene = x$gene[[1]], pas_n = nrow(x),
               effect_vector_cosine = .cosine(x$effect_development, x$effect_validation),
               effect_vector_pearson = pearson, stringsAsFactors = FALSE)
  }))
  d_gene <- development$gene_table[c("gene", "p_value", "q_value", "tested")]
  v_gene <- validation$gene_table[c("gene", "p_value", "q_value", "tested")]
  names(d_gene)[-1] <- paste0(names(d_gene)[-1], "_development")
  names(v_gene)[-1] <- paste0(names(v_gene)[-1], "_validation")
  out <- merge(merge(d_gene, v_gene, by = "gene", all = TRUE), out, by = "gene", all.x = TRUE)
  out$development_hit <- out$tested_development & out$q_value_development <= fdr
  out$validation_hit <- out$tested_validation & out$q_value_validation <= fdr
  out$direction_replicated <- out$effect_vector_cosine > 0
  out$significant_direction_replication <- out$development_hit & out$validation_hit &
    out$direction_replicated
  tested_both <- out$tested_development & out$tested_validation
  rank_correlation <- suppressWarnings(stats::cor(
    -log10(out$p_value_development[tested_both]),
    -log10(out$p_value_validation[tested_both]), method = "spearman",
    use = "pairwise.complete.obs"))
  summary <- data.frame(
    target = development$target,
    tested_in_both = sum(tested_both, na.rm = TRUE),
    development_discoveries = sum(out$development_hit, na.rm = TRUE),
    validation_discoveries = sum(out$validation_hit, na.rm = TRUE),
    directional_replication_rate = if (sum(out$development_hit, na.rm = TRUE))
      mean(out$direction_replicated[out$development_hit], na.rm = TRUE) else NA_real_,
    significant_replication_rate = if (sum(out$development_hit, na.rm = TRUE))
      mean((out$validation_hit & out$direction_replicated)[out$development_hit], na.rm = TRUE) else NA_real_,
    rank_spearman = rank_correlation,
    stringsAsFactors = FALSE
  )
  attr(out, "summary") <- summary
  out
}
