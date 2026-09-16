#' Estimate equal-guide gene-expression effects separately from APA
#'
#' @param rna_counts Gene-by-guide-replicate pseudobulk RNA count matrix.
#' @param sample_data Sample metadata aligned to columns of `rna_counts` and
#'   containing a positive `cell_n` column.
#' @inheritParams repliapa_fit
#' @param replicate,guide,target_col,ntc Column names in `sample_data`.
#' @param cell_n Column containing the number of cells in each sample unit.
#' @return A data frame of replicate-specific and mean log2 fold changes.
#' @export
repliapa_gene_expression_effect <- function(rna_counts, sample_data, target,
                                            replicate = "replicate", guide = "guide",
                                            target_col = "target", ntc = "ntc",
                                            cell_n = "cell_n", replicates = NULL,
                                            pseudocount = 0.5) {
  repliapa_validate_pseudobulk(rna_counts, sample_data, replicate, guide, target_col, ntc)
  .require_columns(sample_data, cell_n, "sample_data")
  cells <- sample_data[[cell_n]]
  if (!is.numeric(cells) || any(!is.finite(cells)) || any(cells <= 0)) {
    stop("`cell_n` must contain positive finite numbers.", call. = FALSE)
  }
  target <- .validate_target_label(target)
  replicates <- .validate_replicates(sample_data, replicate, replicates)
  pseudocount <- .validate_scalar_number(
    pseudocount, "pseudocount", lower = 0, lower_inclusive = FALSE
  )
  target_rows <- as.character(sample_data[[target_col]]) == target & !sample_data[[ntc]] &
    as.character(sample_data[[replicate]]) %in% replicates
  ntc_rows <- sample_data[[ntc]] & as.character(sample_data[[replicate]]) %in% replicates
  target_guides <- sort(.complete_guides(sample_data, replicate, guide, target_rows, replicates))
  ntc_guides <- sort(.complete_guides(sample_data, replicate, guide, ntc_rows, replicates))
  if (length(target_guides) < 2L || length(ntc_guides) < 2L) {
    stop("At least two complete target and NTC guides are required.", call. = FALSE)
  }
  effects <- vapply(replicates, function(rep_value) {
    in_rep <- as.character(sample_data[[replicate]]) == rep_value
    target_columns <- in_rep & as.character(sample_data[[guide]]) %in% target_guides
    ntc_columns <- in_rep & as.character(sample_data[[guide]]) %in% ntc_guides
    target_rate <- sweep(as.matrix(rna_counts[, target_columns, drop = FALSE]), 2L,
                         cells[target_columns], "/")
    ntc_rate <- sweep(as.matrix(rna_counts[, ntc_columns, drop = FALSE]), 2L,
                      cells[ntc_columns], "/")
    log2((rowMeans(target_rate) + pseudocount) / (rowMeans(ntc_rate) + pseudocount))
  }, numeric(nrow(rna_counts)))
  if (is.null(dim(effects))) effects <- matrix(effects, ncol = 1L)
  colnames(effects) <- paste0("log2FC_", replicates)
  data.frame(gene = rownames(rna_counts), log2FC = rowMeans(effects), effects,
             stringsAsFactors = FALSE, check.names = FALSE)
}

#' Classify APA and gene-expression effects into four quadrants
#'
#' @param fit A `repliapa_fit` object.
#' @param gene_expression Output of [repliapa_gene_expression_effect()] or a
#'   data frame containing `gene` and `log2FC`.
#' @param apa_fdr BH threshold for the APA global test.
#' @param ge_abs_log2fc Absolute gene-expression effect threshold.
#' @param ge_significant Optional logical vector overriding the effect threshold.
#' @return A data frame with `APA_ONLY`, `GE_ONLY`, `APA_AND_GE`, or `NEITHER`.
#' @export
repliapa_classify_effects <- function(fit, gene_expression, apa_fdr = 0.05,
                                      ge_abs_log2fc = 0.25, ge_significant = NULL) {
  if (!inherits(fit, "repliapa_fit")) stop("`fit` must be a repliapa_fit object.", call. = FALSE)
  if (!is.data.frame(gene_expression)) {
    stop("`gene_expression` must be a data frame.", call. = FALSE)
  }
  .require_columns(gene_expression, c("gene", "log2FC"), "gene_expression")
  apa_fdr <- .validate_scalar_number(
    apa_fdr, "apa_fdr", lower = 0, lower_inclusive = FALSE,
    upper = 1, upper_inclusive = FALSE
  )
  ge_abs_log2fc <- .validate_scalar_number(
    ge_abs_log2fc, "ge_abs_log2fc", lower = 0
  )
  gene_value <- as.character(gene_expression$gene)
  if (anyNA(gene_value) || any(gene_value == "") || anyDuplicated(gene_value)) {
    stop("Gene-expression identifiers must be unique and non-missing.", call. = FALSE)
  }
  if (!is.numeric(gene_expression$log2FC) ||
      any(!is.finite(gene_expression$log2FC))) {
    stop("Gene-expression log2 fold changes must be finite numbers.", call. = FALSE)
  }
  if (!is.null(ge_significant) &&
      (!is.logical(ge_significant) || length(ge_significant) != nrow(gene_expression) ||
       anyNA(ge_significant))) {
    stop("`ge_significant` must be a logical vector aligned to gene_expression.", call. = FALSE)
  }
  ge <- gene_expression[c("gene", "log2FC")]
  ge$GE_significant <- if (is.null(ge_significant)) abs(ge$log2FC) >= ge_abs_log2fc else ge_significant
  out <- merge(fit$gene_table, ge, by = "gene", all.x = TRUE)
  out$APA_significant <- out$tested & out$q_value <= apa_fdr
  out$GE_significant[is.na(out$GE_significant)] <- FALSE
  out$effect_class <- ifelse(out$APA_significant & out$GE_significant, "APA_AND_GE",
    ifelse(out$APA_significant, "APA_ONLY",
      ifelse(out$GE_significant, "GE_ONLY", "NEITHER")))
  out$GE_rule <- if (is.null(ge_significant))
    paste0("abs_log2FC_at_least_", ge_abs_log2fc) else "user_supplied_significance"
  out
}
