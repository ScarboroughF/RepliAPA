#' Export a frozen pseudobulk contrast as DRIMSeq input
#'
#' Creates the two tables expected by `DRIMSeq::dmDSdata()` without fitting or
#' filtering the baseline. Target and non-targeting guide identities remain
#' attached to their biological-replicate rows. Only guides observed in every
#' requested replicate are exported.
#'
#' @inheritParams repliapa_fit
#' @param pas_ids Optional predeclared PAS identifiers to export. Every
#'   identifier must exist; the function never silently intersects site sets.
#' @return A `repliapa_drimseq_input` object containing `counts`, `samples`,
#'   `design`, and an explicit `accounting` table. The `counts` and `samples`
#'   components can be passed directly to `DRIMSeq::dmDSdata()`.
#' @export
repliapa_drimseq_input <- function(counts, sample_data, pas_data, target,
                                   gene = "gene", replicate = "replicate",
                                   guide = "guide", target_col = "target",
                                   ntc = "ntc", replicates = NULL,
                                   pas_ids = NULL) {
  repliapa_validate_pseudobulk(counts, sample_data, replicate, guide, target_col, ntc)
  .aligned_data(pas_data, rownames(counts), "pas_data")
  .require_columns(pas_data, gene, "pas_data")
  gene_value <- as.character(pas_data[[gene]])
  if (anyNA(gene_value) || any(gene_value == "")) {
    stop("The PAS-to-gene map is incomplete.", call. = FALSE)
  }
  target <- .validate_target_label(target)
  replicates <- .validate_replicates(sample_data, replicate, replicates)

  target_rows <- as.character(sample_data[[target_col]]) == target &
    !sample_data[[ntc]] & as.character(sample_data[[replicate]]) %in% replicates
  ntc_rows <- sample_data[[ntc]] &
    as.character(sample_data[[replicate]]) %in% replicates
  target_guides <- sort(.complete_guides(
    sample_data, replicate, guide, target_rows, replicates
  ))
  ntc_guides <- sort(.complete_guides(
    sample_data, replicate, guide, ntc_rows, replicates
  ))
  if (length(target_guides) < 2L) {
    stop("At least two target guides complete across selected replicates are required.",
         call. = FALSE)
  }
  if (length(ntc_guides) < 2L) {
    stop("At least two NTC guides complete across selected replicates are required.",
         call. = FALSE)
  }

  keep_sample <- as.character(sample_data[[replicate]]) %in% replicates &
    as.character(sample_data[[guide]]) %in% c(target_guides, ntc_guides)
  context_data <- sample_data[keep_sample, , drop = FALSE]
  context_data <- context_data[
    order(as.character(context_data[[replicate]]),
          match(as.character(context_data[[guide]]), c(target_guides, ntc_guides))),
    , drop = FALSE
  ]

  if (is.null(pas_ids)) {
    pas_ids <- rownames(counts)
  } else {
    pas_ids <- as.character(pas_ids)
    if (!length(pas_ids) || anyNA(pas_ids) || any(pas_ids == "") || anyDuplicated(pas_ids)) {
      stop("`pas_ids` must contain unique non-missing identifiers.", call. = FALSE)
    }
    missing_pas <- setdiff(pas_ids, rownames(counts))
    if (length(missing_pas)) {
      stop("Requested PAS identifiers are unavailable: ",
           paste(utils::head(missing_pas, 5L), collapse = ", "), call. = FALSE)
    }
  }
  pas_index <- match(pas_ids, rownames(counts))
  context_counts <- counts[pas_index, rownames(context_data), drop = FALSE]
  context_pas <- pas_data[pas_index, , drop = FALSE]

  samples <- data.frame(
    sample_id = rownames(context_data),
    replicate = factor(as.character(context_data[[replicate]]), levels = replicates),
    guide = as.character(context_data[[guide]]),
    group = as.integer(as.character(context_data[[guide]]) %in% target_guides),
    stringsAsFactors = FALSE,
    row.names = rownames(context_data)
  )
  design <- if (length(replicates) > 1L) {
    stats::model.matrix(~ replicate + group, data = samples)
  } else {
    stats::model.matrix(~ group, data = samples)
  }
  drim_counts <- data.frame(
    gene_id = as.character(context_pas[[gene]]),
    feature_id = rownames(context_pas),
    as.matrix(context_counts),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  exported_total <- sum(as.matrix(drim_counts[, -(1:2), drop = FALSE]))
  context_total <- sum(context_counts)
  if (!isTRUE(all.equal(exported_total, context_total, tolerance = 0))) {
    stop("DRIMSeq export did not conserve selected pseudobulk counts.", call. = FALSE)
  }

  accounting <- data.frame(
    target = target,
    replicate_n = length(replicates),
    target_guide_n = length(target_guides),
    ntc_guide_n = length(ntc_guides),
    input_sample_n = ncol(counts),
    exported_sample_n = nrow(samples),
    excluded_context_sample_n = ncol(counts) - nrow(samples),
    input_pas_n = nrow(counts),
    exported_pas_n = nrow(drim_counts),
    excluded_pas_n = nrow(counts) - nrow(drim_counts),
    input_gene_n = length(unique(gene_value)),
    exported_gene_n = length(unique(drim_counts$gene_id)),
    exported_total_count = as.numeric(exported_total),
    export_filtering = "NONE",
    count_conserved = TRUE,
    stringsAsFactors = FALSE
  )
  output <- list(
    counts = drim_counts,
    samples = samples,
    design = design,
    accounting = accounting,
    target = target,
    replicates = replicates,
    target_guides = target_guides,
    ntc_guides = ntc_guides
  )
  class(output) <- "repliapa_drimseq_input"
  output
}

#' @export
print.repliapa_drimseq_input <- function(x, ...) {
  cat("RepliAPA DRIMSeq input\n")
  cat("  target: ", x$target, "\n", sep = "")
  cat("  target / NTC guides: ", length(x$target_guides), " / ",
      length(x$ntc_guides), "\n", sep = "")
  cat("  samples / PAS / genes: ", nrow(x$samples), " / ",
      nrow(x$counts), " / ", length(unique(x$counts$gene_id)), "\n", sep = "")
  cat("  export filtering: NONE; counts conserved: TRUE\n")
  invisible(x)
}

#' Compare a RepliAPA fit with a sample-aware baseline
#'
#' Aligns gene-level results without hiding method-native filtering. Native
#' summaries use each method's own tested set; matched summaries use only genes
#' tested by both methods. This function does not run or reinterpret the
#' baseline model.
#'
#' @param repliapa A fitted `repliapa_fit` object.
#' @param baseline A data frame containing one row per baseline gene.
#' @param baseline_gene,baseline_p Column names for baseline gene identifiers
#'   and nominal P values.
#' @param baseline_q Optional column containing baseline adjusted P values. If
#'   `NULL`, BH values are computed over the baseline-native tested set.
#' @param baseline_tested Optional logical tested-status column. If `NULL`, a
#'   finite nominal P value defines a tested gene.
#' @param baseline_name Non-empty label for the baseline method.
#' @param alpha Discovery threshold applied to both methods' adjusted P values.
#' @return A `repliapa_baseline_comparison` containing aligned gene results,
#'   method-native accounting, matched-set summaries, and agreement metrics.
#' @export
repliapa_compare_baseline <- function(repliapa, baseline,
                                      baseline_gene = "gene_id",
                                      baseline_p = "pvalue",
                                      baseline_q = NULL,
                                      baseline_tested = NULL,
                                      baseline_name = "DRIMSeq",
                                      alpha = 0.05) {
  if (!inherits(repliapa, "repliapa_fit")) {
    stop("`repliapa` must be a repliapa_fit object.", call. = FALSE)
  }
  if (!is.data.frame(baseline)) stop("`baseline` must be a data.frame.", call. = FALSE)
  .require_columns(baseline, c(baseline_gene, baseline_p), "baseline")
  if (!is.null(baseline_q)) .require_columns(baseline, baseline_q, "baseline")
  if (!is.null(baseline_tested)) .require_columns(baseline, baseline_tested, "baseline")
  if (length(baseline_name) != 1L || is.na(baseline_name) || baseline_name == "") {
    stop("`baseline_name` must be one non-missing label.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one number strictly between zero and one.", call. = FALSE)
  }

  baseline_gene_value <- as.character(baseline[[baseline_gene]])
  if (anyNA(baseline_gene_value) || any(baseline_gene_value == "") ||
      anyDuplicated(baseline_gene_value)) {
    stop("Baseline gene identifiers must be unique and non-missing.", call. = FALSE)
  }
  if (!is.numeric(baseline[[baseline_p]])) {
    stop("The baseline P-value column must be numeric.", call. = FALSE)
  }
  baseline_p_value <- as.numeric(baseline[[baseline_p]])
  baseline_tested_value <- if (is.null(baseline_tested)) {
    is.finite(baseline_p_value)
  } else {
    value <- baseline[[baseline_tested]]
    if (!is.logical(value) || anyNA(value)) {
      stop("The baseline tested-status column must be non-missing logical.", call. = FALSE)
    }
    value
  }
  if (any(!is.finite(baseline_p_value[baseline_tested_value])) ||
      any(baseline_p_value[baseline_tested_value] < 0 |
          baseline_p_value[baseline_tested_value] > 1)) {
    stop("Tested baseline P values must be finite and between zero and one.", call. = FALSE)
  }
  baseline_q_value <- rep(NA_real_, nrow(baseline))
  if (is.null(baseline_q)) {
    baseline_q_value[baseline_tested_value] <- stats::p.adjust(
      baseline_p_value[baseline_tested_value], method = "BH"
    )
  } else {
    if (!is.numeric(baseline[[baseline_q]])) {
      stop("The baseline adjusted P-value column must be numeric.", call. = FALSE)
    }
    baseline_q_value <- as.numeric(baseline[[baseline_q]])
    if (any(!is.finite(baseline_q_value[baseline_tested_value])) ||
        any(baseline_q_value[baseline_tested_value] < 0 |
            baseline_q_value[baseline_tested_value] > 1)) {
      stop("Tested baseline adjusted P values must be finite and between zero and one.",
           call. = FALSE)
    }
  }

  repliapa_table <- repliapa$gene_table
  .require_columns(repliapa_table,
                   c("gene", "p_value", "q_value", "tested", "filtered", "abstained", "reason"),
                   "repliapa$gene_table")
  if (anyDuplicated(as.character(repliapa_table$gene))) {
    stop("RepliAPA gene results are not unique.", call. = FALSE)
  }
  repliapa_aligned <- data.frame(
    gene = as.character(repliapa_table$gene),
    repliapa_present = TRUE,
    repliapa_tested = repliapa_table$tested,
    repliapa_p_value = repliapa_table$p_value,
    repliapa_q_value = repliapa_table$q_value,
    repliapa_filtered = repliapa_table$filtered,
    repliapa_abstained = repliapa_table$abstained,
    repliapa_reason = as.character(repliapa_table$reason),
    stringsAsFactors = FALSE
  )
  baseline_aligned <- data.frame(
    gene = baseline_gene_value,
    baseline_present = TRUE,
    baseline_tested = baseline_tested_value,
    baseline_p_value = baseline_p_value,
    baseline_q_value = baseline_q_value,
    stringsAsFactors = FALSE
  )
  aligned <- merge(repliapa_aligned, baseline_aligned, by = "gene", all = TRUE, sort = TRUE)
  aligned$repliapa_present[is.na(aligned$repliapa_present)] <- FALSE
  aligned$baseline_present[is.na(aligned$baseline_present)] <- FALSE
  aligned$repliapa_tested[is.na(aligned$repliapa_tested)] <- FALSE
  aligned$baseline_tested[is.na(aligned$baseline_tested)] <- FALSE
  aligned$repliapa_filtered[is.na(aligned$repliapa_filtered)] <- FALSE
  aligned$repliapa_abstained[is.na(aligned$repliapa_abstained)] <- FALSE
  aligned$repliapa_significant <- aligned$repliapa_tested &
    !is.na(aligned$repliapa_q_value) & aligned$repliapa_q_value <= alpha
  aligned$baseline_significant <- aligned$baseline_tested &
    !is.na(aligned$baseline_q_value) & aligned$baseline_q_value <= alpha
  aligned$matched_tested <- aligned$repliapa_tested & aligned$baseline_tested

  native_summary <- data.frame(
    method = c("RepliAPA", baseline_name),
    native_gene_n = c(nrow(repliapa_aligned), nrow(baseline_aligned)),
    tested_gene_n = c(sum(repliapa_aligned$repliapa_tested), sum(baseline_tested_value)),
    untested_gene_n = c(sum(!repliapa_aligned$repliapa_tested), sum(!baseline_tested_value)),
    filtered_gene_n = c(sum(repliapa_aligned$repliapa_filtered), NA_integer_),
    abstained_gene_n = c(sum(repliapa_aligned$repliapa_abstained), NA_integer_),
    discovery_n = c(
      sum(repliapa_aligned$repliapa_tested &
            repliapa_aligned$repliapa_q_value <= alpha, na.rm = TRUE),
      sum(baseline_tested_value & baseline_q_value <= alpha, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  matched <- aligned$matched_tested
  matched_n <- sum(matched)
  matched_summary <- data.frame(
    method = c("RepliAPA", baseline_name),
    matched_tested_gene_n = matched_n,
    discovery_n = c(sum(aligned$repliapa_significant[matched]),
                    sum(aligned$baseline_significant[matched])),
    discovery_rate = if (matched_n) {
      c(mean(aligned$repliapa_significant[matched]),
        mean(aligned$baseline_significant[matched]))
    } else c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  both_discovered <- aligned$repliapa_significant[matched] &
    aligned$baseline_significant[matched]
  either_discovered <- aligned$repliapa_significant[matched] |
    aligned$baseline_significant[matched]
  rank_correlation <- if (matched_n >= 2L &&
      stats::sd(aligned$repliapa_p_value[matched]) > 0 &&
      stats::sd(aligned$baseline_p_value[matched]) > 0) {
    stats::cor(aligned$repliapa_p_value[matched], aligned$baseline_p_value[matched],
               method = "spearman")
  } else NA_real_
  agreement <- data.frame(
    matched_tested_gene_n = matched_n,
    rank_spearman = rank_correlation,
    discovery_overlap_n = sum(both_discovered),
    discovery_union_n = sum(either_discovered),
    discovery_jaccard = if (any(either_discovered)) {
      sum(both_discovered) / sum(either_discovered)
    } else NA_real_,
    repliapa_only_discovery_n = sum(
      aligned$repliapa_significant[matched] & !aligned$baseline_significant[matched]
    ),
    baseline_only_discovery_n = sum(
      !aligned$repliapa_significant[matched] & aligned$baseline_significant[matched]
    ),
    stringsAsFactors = FALSE
  )
  output <- list(
    gene_table = aligned,
    native_summary = native_summary,
    matched_summary = matched_summary,
    agreement = agreement,
    baseline_name = baseline_name,
    alpha = alpha
  )
  class(output) <- "repliapa_baseline_comparison"
  output
}

#' @export
print.repliapa_baseline_comparison <- function(x, ...) {
  cat("RepliAPA baseline comparison\n")
  cat("  baseline: ", x$baseline_name, "\n", sep = "")
  cat("  matched tested genes: ", x$agreement$matched_tested_gene_n, "\n", sep = "")
  cat("  native discoveries (RepliAPA / baseline): ",
      paste(x$native_summary$discovery_n, collapse = " / "), "\n", sep = "")
  cat("  matched discoveries (RepliAPA / baseline): ",
      paste(x$matched_summary$discovery_n, collapse = " / "), "\n", sep = "")
  invisible(x)
}
