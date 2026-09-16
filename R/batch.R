#' Fit multiple perturbation targets with explicit support accounting
#'
#' @inheritParams repliapa_fit
#' @param targets Target labels to fit. By default, all non-NTC target labels.
#' @param strict If `TRUE`, stop at the first unsupported target. Otherwise,
#'   record it in `accounting` and continue.
#' @return A `repliapa_fit_list` containing named fits and target accounting.
#' @export
repliapa_fit_targets <- function(counts, sample_data, pas_data, targets = NULL,
                                 gene = "gene", replicate = "replicate", guide = "guide",
                                 target_col = "target", ntc = "ntc", replicates = NULL,
                                 pseudocount = 0.5, min_gene_total = 20,
                                 max_combinations = 1000000, strict = FALSE) {
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
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  replicates <- .validate_replicates(sample_data, replicate, replicates)
  if (is.null(targets)) {
    targets <- sort(unique(as.character(sample_data[[target_col]][!sample_data[[ntc]]])))
  }
  targets <- unique(as.character(targets))
  if (!length(targets) || anyNA(targets) || any(targets == "")) {
    stop("`targets` must contain non-missing target labels.", call. = FALSE)
  }
  target_results <- lapply(targets, function(target_value) {
    value <- tryCatch(
      repliapa_fit(counts, sample_data, pas_data, target = target_value,
                   gene = gene, replicate = replicate, guide = guide,
                   target_col = target_col, ntc = ntc, replicates = replicates,
                   pseudocount = pseudocount, min_gene_total = min_gene_total,
                   max_combinations = max_combinations),
      error = identity
    )
    if (inherits(value, "error")) {
      if (isTRUE(strict)) stop(value)
      return(list(
        fit = NULL,
        accounting = data.frame(
          target = target_value, status = "UNSUPPORTED",
          reason = conditionMessage(value), target_guide_n = NA_integer_,
          ntc_guide_n = NA_integer_, tested_gene_n = 0L,
          filtered_gene_n = 0L, abstained_gene_n = 0L,
          discovery_n = 0L, stringsAsFactors = FALSE
        )
      ))
    }
    list(
      fit = value,
      accounting = data.frame(
        target = target_value, status = "FITTED", reason = "NONE",
        target_guide_n = length(value$target_guides),
        ntc_guide_n = length(value$ntc_guides),
        tested_gene_n = sum(value$gene_table$tested),
        filtered_gene_n = sum(value$gene_table$filtered),
        abstained_gene_n = sum(value$gene_table$abstained),
        discovery_n = sum(value$gene_table$significant, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    )
  })
  fitted <- vapply(target_results, function(x) !is.null(x$fit), logical(1))
  fits <- lapply(target_results[fitted], `[[`, "fit")
  names(fits) <- targets[fitted]
  accounting <- do.call(rbind, lapply(target_results, `[[`, "accounting"))
  output <- list(fits = fits, accounting = accounting, targets = targets)
  class(output) <- "repliapa_fit_list"
  output
}

#' @export
print.repliapa_fit_list <- function(x, ...) {
  cat("RepliAPA target batch\n")
  cat("  fitted / unsupported: ", sum(x$accounting$status == "FITTED"), " / ",
      sum(x$accounting$status == "UNSUPPORTED"), "\n", sep = "")
  cat("  tested genes: ", sum(x$accounting$tested_gene_n), "\n", sep = "")
  cat("  BH discoveries: ", sum(x$accounting$discovery_n), "\n", sep = "")
  invisible(x)
}
