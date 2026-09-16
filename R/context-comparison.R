.normalize_fit_collection <- function(x, label) {
  if (inherits(x, "repliapa_fit_list")) {
    fits <- x$fits
    accounting <- x$accounting
  } else if (inherits(x, "repliapa_fit")) {
    fits <- stats::setNames(list(x), x$target)
    accounting <- data.frame(target = x$target, status = "FITTED",
                             reason = "NONE", stringsAsFactors = FALSE)
  } else if (is.list(x) && length(x) &&
             all(vapply(x, inherits, logical(1), what = "repliapa_fit"))) {
    fits <- x
    fit_targets <- vapply(fits, `[[`, character(1), "target")
    if (anyNA(fit_targets) || any(fit_targets == "") || anyDuplicated(fit_targets)) {
      stop("Fit target labels must be complete and unique.", call. = FALSE)
    }
    if (!is.null(names(fits)) && all(names(fits) != "") &&
        !identical(as.character(names(fits)), fit_targets)) {
      stop("Named fit collections must use names identical to fit target labels.",
           call. = FALSE)
    }
    names(fits) <- fit_targets
    accounting <- data.frame(target = fit_targets, status = "FITTED",
                             reason = "NONE", stringsAsFactors = FALSE)
  } else {
    stop("`", label, "` must contain RepliAPA fits.", call. = FALSE)
  }
  if (anyDuplicated(names(fits)) || anyDuplicated(accounting$target)) {
    stop("Target labels must be unique within each context.", call. = FALSE)
  }
  list(fits = fits, accounting = accounting)
}

.context_target_status <- function(collection, target) {
  row <- collection$accounting[as.character(collection$accounting$target) == target,
                               , drop = FALSE]
  if (!nrow(row)) return(list(status = "NOT_REQUESTED", reason = "TARGET_NOT_REQUESTED"))
  list(status = as.character(row$status[[1]]), reason = as.character(row$reason[[1]]))
}

.compare_context_target <- function(primary, external, primary_label, external_label, fdr) {
  primary_gene <- primary$gene_table[c("gene", "pas_n", "p_value", "q_value", "tested",
                                       "filtered", "abstained", "effect_norm",
                                       "coordinate_shift")]
  external_gene <- external$gene_table[c("gene", "pas_n", "p_value", "q_value", "tested",
                                         "filtered", "abstained", "effect_norm",
                                         "coordinate_shift")]
  names(primary_gene)[-1] <- paste0(primary_label, "_", names(primary_gene)[-1])
  names(external_gene)[-1] <- paste0(external_label, "_", names(external_gene)[-1])
  out <- merge(primary_gene, external_gene, by = "gene", all = TRUE)

  primary_effects <- split(primary$pas_effects, primary$pas_effects$gene)
  external_effects <- split(external$pas_effects, external$pas_effects$gene)
  common_pas_n <- integer(nrow(out))
  cosine <- rep(NA_real_, nrow(out))
  site_status <- character(nrow(out))
  for (i in seq_len(nrow(out))) {
    gene_value <- out$gene[[i]]
    first <- primary_effects[[gene_value]]
    second <- external_effects[[gene_value]]
    if (is.null(first) || is.null(second)) {
      site_status[[i]] <- "MISSING_TESTED_EFFECT"
      next
    }
    joined <- merge(first[c("pas_id", "effect")], second[c("pas_id", "effect")],
                    by = "pas_id", suffixes = c("_primary", "_external"))
    common_pas_n[[i]] <- nrow(joined)
    if (nrow(joined) < 2L) {
      site_status[[i]] <- "INSUFFICIENT_COMMON_PAS"
      next
    }
    site_status[[i]] <- "SUPPORTED_COMMON_SITE_SET"
    cosine[[i]] <- .cosine(joined$effect_primary, joined$effect_external)
  }
  out$common_pas_n <- common_pas_n
  out$site_status <- site_status
  out$effect_vector_cosine <- cosine
  primary_tested <- out[[paste0(primary_label, "_tested")]]
  external_tested <- out[[paste0(external_label, "_tested")]]
  primary_q <- out[[paste0(primary_label, "_q_value")]]
  external_q <- out[[paste0(external_label, "_q_value")]]
  comparable <- !is.na(primary_tested) & !is.na(external_tested) &
    primary_tested & external_tested & site_status == "SUPPORTED_COMMON_SITE_SET" &
    is.finite(cosine)
  primary_significant <- comparable & primary_q <= fdr
  external_significant <- comparable & external_q <= fdr
  out$category <- ifelse(
    !comparable, "UNSUPPORTED_GENE_CONTEXT",
    ifelse(primary_significant & external_significant & cosine > 0, "SHARED_EFFECT",
      ifelse(primary_significant & external_significant & cosine < 0, "OPPOSITE_EFFECT",
        ifelse(xor(primary_significant, external_significant),
               "CELL_LINE_SPECIFIC_EFFECT", "SUPPORTED_NEITHER_SIGNIFICANT")))
  )
  out
}

#' Compare APA effects across independently fitted cell contexts
#'
#' Effects are aligned only on PAS observed in both contexts. Each context keeps
#' its own guide-preserving fit and BH correction; no external-context result is
#' used to refit or tune the primary model.
#'
#' @param primary,external A `repliapa_fit`, `repliapa_fit_list`, or named list
#'   of fits for the primary and external contexts.
#' @param primary_label,external_label Safe labels used as output prefixes.
#' @param targets Optional target universe. By default, all requested targets in
#'   either context are retained, including unsupported targets.
#' @param fdr Within-context BH threshold.
#' @return A `repliapa_context_comparison` object containing gene-level effects,
#'   target summaries, and explicit target support accounting.
#' @export
repliapa_compare_contexts <- function(primary, external,
                                      primary_label = "primary",
                                      external_label = "external",
                                      targets = NULL, fdr = 0.05) {
  primary_collection <- .normalize_fit_collection(primary, "primary")
  external_collection <- .normalize_fit_collection(external, "external")
  labels <- c(primary_label, external_label)
  if (anyNA(labels) || any(labels == "") || anyDuplicated(labels) ||
      any(make.names(labels) != labels)) {
    stop("Context labels must be distinct syntactic names.", call. = FALSE)
  }
  if (!is.numeric(fdr) || length(fdr) != 1L || !is.finite(fdr) || fdr <= 0 || fdr >= 1) {
    stop("`fdr` must be between zero and one.", call. = FALSE)
  }
  if (is.null(targets)) {
    targets <- sort(unique(c(as.character(primary_collection$accounting$target),
                             as.character(external_collection$accounting$target))))
  }
  targets <- unique(as.character(targets))
  if (!length(targets) || anyNA(targets) || any(targets == "")) {
    stop("`targets` must contain non-missing target labels.", call. = FALSE)
  }

  target_accounting <- do.call(rbind, lapply(targets, function(target_value) {
    first <- .context_target_status(primary_collection, target_value)
    second <- .context_target_status(external_collection, target_value)
    first_fitted <- identical(first$status, "FITTED") &&
      target_value %in% names(primary_collection$fits)
    second_fitted <- identical(second$status, "FITTED") &&
      target_value %in% names(external_collection$fits)
    comparison_status <- if (first_fitted && second_fitted) "SUPPORTED_BOTH" else if (
      first_fitted
    ) "EXTERNAL_UNSUPPORTED" else if (second_fitted) "PRIMARY_UNSUPPORTED" else
      "BOTH_UNSUPPORTED"
    data.frame(
      target = target_value,
      primary_status = first$status, primary_reason = first$reason,
      external_status = second$status, external_reason = second$reason,
      comparison_status = comparison_status,
      stringsAsFactors = FALSE
    )
  }))

  supported <- target_accounting$target[
    target_accounting$comparison_status == "SUPPORTED_BOTH"
  ]
  results <- do.call(rbind, lapply(supported, function(target_value) {
    value <- .compare_context_target(
      primary_collection$fits[[target_value]], external_collection$fits[[target_value]],
      primary_label, external_label, fdr
    )
    value$target <- target_value
    value
  }))
  if (!is.null(results) && nrow(results)) {
    results <- results[c("target", setdiff(names(results), "target"))]
    rownames(results) <- NULL
  } else {
    results <- data.frame()
  }

  target_summary <- do.call(rbind, lapply(supported, function(target_value) {
    x <- results[results$target == target_value, , drop = FALSE]
    comparable <- x$category != "UNSUPPORTED_GENE_CONTEXT"
    primary_p <- x[[paste0(primary_label, "_p_value")]][comparable]
    external_p <- x[[paste0(external_label, "_p_value")]][comparable]
    data.frame(
      target = target_value,
      gene_n = nrow(x), comparable_gene_n = sum(comparable),
      unsupported_gene_n = sum(!comparable),
      rank_Spearman = suppressWarnings(stats::cor(
        -log10(primary_p), -log10(external_p), method = "spearman",
        use = "pairwise.complete.obs"
      )),
      median_effect_cosine = stats::median(x$effect_vector_cosine[comparable], na.rm = TRUE),
      positive_cosine_fraction = mean(x$effect_vector_cosine[comparable] > 0, na.rm = TRUE),
      shared_n = sum(x$category == "SHARED_EFFECT"),
      opposite_n = sum(x$category == "OPPOSITE_EFFECT"),
      context_specific_n = sum(x$category == "CELL_LINE_SPECIFIC_EFFECT"),
      neither_n = sum(x$category == "SUPPORTED_NEITHER_SIGNIFICANT"),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(target_summary)) target_summary <- data.frame()
  output <- list(
    results = results, target_summary = target_summary,
    target_accounting = target_accounting,
    contexts = labels, fdr = fdr
  )
  class(output) <- "repliapa_context_comparison"
  output
}

#' @export
print.repliapa_context_comparison <- function(x, ...) {
  cat("RepliAPA cross-context comparison\n")
  cat("  contexts: ", paste(x$contexts, collapse = " vs "), "\n", sep = "")
  cat("  targets supported / unsupported: ",
      sum(x$target_accounting$comparison_status == "SUPPORTED_BOTH"), " / ",
      sum(x$target_accounting$comparison_status != "SUPPORTED_BOTH"), "\n", sep = "")
  if (nrow(x$results)) {
    cat("  shared / opposite / context-specific: ",
        sum(x$results$category == "SHARED_EFFECT"), " / ",
        sum(x$results$category == "OPPOSITE_EFFECT"), " / ",
        sum(x$results$category == "CELL_LINE_SPECIFIC_EFFECT"), "\n", sep = "")
  }
  invisible(x)
}
