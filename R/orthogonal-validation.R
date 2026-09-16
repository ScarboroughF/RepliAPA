#' Derive proximal-direction effects from all-PAS fits
#'
#' The distal-bin shift is the sum of all-PAS effects assigned to the distal
#' half of the direction-aware site set. The proximal direction is its negative.
#'
#' @param fits A `repliapa_fit`, `repliapa_fit_list`, or named list of fits.
#' @param pd_data A `repliapa_pd_data` from
#'   [repliapa_collapse_proximal_distal()] using the same PAS set.
#' @return A target-by-gene table with distal-bin shift, proximal direction,
#'   within-fit q-value, and tested status.
#' @export
repliapa_proximal_directions <- function(fits, pd_data) {
  collection <- .normalize_fit_collection(fits, "fits")
  if (!inherits(pd_data, "repliapa_pd_data")) {
    stop("`pd_data` must be a repliapa_pd_data object.", call. = FALSE)
  }
  distal <- pd_data$membership[pd_data$membership$bin == "distal",
                               c("gene", "pas_id"), drop = FALSE]
  results <- do.call(rbind, lapply(collection$fits, function(fit) {
    effects <- merge(distal, fit$pas_effects[c("gene", "pas_id", "effect")],
                     by = c("gene", "pas_id"))
    shifts <- if (nrow(effects)) stats::aggregate(effect ~ gene, effects, sum) else
      data.frame(gene = character(), effect = numeric())
    names(shifts)[2] <- "distal_bin_shift"
    gene_table <- fit$gene_table[c("gene", "q_value", "tested")]
    out <- merge(gene_table, shifts, by = "gene", all.x = TRUE)
    out$target <- fit$target
    out$proximal_direction <- -out$distal_bin_shift
    out[c("target", "gene", "distal_bin_shift", "proximal_direction",
          "q_value", "tested")]
  }))
  rownames(results) <- NULL
  results
}

.orthogonal_usage <- function(data, group_columns, proximal, distal, total) {
  groups <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  pieces <- split(data, groups)
  do.call(rbind, lapply(pieces, function(x) {
    identity <- x[1, group_columns, drop = FALSE]
    proximal_count <- sum(x[[proximal]])
    distal_count <- sum(x[[distal]])
    total_count <- sum(x[[total]])
    denominator <- proximal_count + distal_count
    cbind(identity,
          proximal_count = proximal_count,
          distal_count = distal_count,
          total_count = total_count,
          row_n = nrow(x),
          proximal_usage = if (denominator > 0) proximal_count / denominator else NA_real_,
          row.names = NULL)
  }))
}

#' Validate frozen APA directions against an orthogonal count assay
#'
#' This function is deliberately post-freeze. It aggregates raw orthogonal
#' rows within construct, replicate, and perturbation; pairs each perturbation
#' with the control in the same orthogonal replicate; and compares only effect
#' directions. Orthogonal outcomes never refit the APA model.
#'
#' @param effect_table Output from [repliapa_proximal_directions()] or a data
#'   frame with `target`, `gene`, `proximal_direction`, `q_value`, and `tested`.
#' @param orthogonal_data Row-level orthogonal proximal/distal count data.
#' @param freeze Verified analysis freeze created before orthogonal outcomes.
#' @param orthogonal_registry_id Registry entry carrying the orthogonal role.
#' @param control Perturbation label identifying the orthogonal control.
#' @param targets Optional perturbations to validate. Defaults to effect targets.
#' @param construct_columns Columns uniquely defining an orthogonal construct;
#'   the first must map to the `gene` column in `effect_table`.
#' @param replicate,perturbation,proximal,distal,total Column names.
#' @param experiment Optional scalar experiment value to retain.
#' @param experiment_col Column containing the experiment label.
#' @param min_replicates Minimum paired orthogonal replicates per construct.
#' @return A `repliapa_orthogonal_validation` object with construct-level
#'   results, target summaries, and mapping/support accounting.
#' @export
repliapa_orthogonal_validation <- function(
    effect_table, orthogonal_data, freeze, orthogonal_registry_id, control = "NT",
    targets = NULL,
    construct_columns = c("gene_id", "pas_id", "aim", "subaim"),
    replicate = "replicate", perturbation = "perturbation",
    proximal = "proximal", distal = "distal", total = "total",
    experiment = NULL, experiment_col = "experiment", min_replicates = 2L) {
  repliapa_authorize(freeze, orthogonal_registry_id, "orthogonal_validation")
  if (!is.data.frame(effect_table) || !is.data.frame(orthogonal_data)) {
    stop("Inputs must be data frames.", call. = FALSE)
  }
  .require_columns(effect_table,
                   c("target", "gene", "proximal_direction", "q_value", "tested"),
                   "effect_table")
  if (anyDuplicated(effect_table[c("target", "gene")])) {
    stop("`effect_table` must have one row per target and gene.", call. = FALSE)
  }
  if (!is.logical(effect_table$tested) ||
      !is.numeric(effect_table$proximal_direction) || !is.numeric(effect_table$q_value)) {
    stop("Frozen effect columns have invalid types.", call. = FALSE)
  }
  required <- c(construct_columns, replicate, perturbation, proximal, distal, total)
  if (!is.null(experiment)) required <- c(required, experiment_col)
  .require_columns(orthogonal_data, required, "orthogonal_data")
  if (!length(construct_columns) || anyDuplicated(construct_columns)) {
    stop("`construct_columns` must contain unique column names.", call. = FALSE)
  }
  count_columns <- c(proximal, distal, total)
  valid_counts <- vapply(orthogonal_data[count_columns], function(x) {
    is.numeric(x) && all(is.finite(x)) && all(x >= 0) && all(abs(x - round(x)) < 1e-8)
  }, logical(1))
  if (!all(valid_counts)) {
    stop("Orthogonal proximal, distal, and total counts must be nonnegative integers.",
         call. = FALSE)
  }
  min_replicates <- .validate_scalar_integer(
    min_replicates, "min_replicates", minimum = 2L
  )
  data <- orthogonal_data
  if (!is.null(experiment)) {
    if (length(experiment) != 1L || is.na(experiment)) {
      stop("`experiment` must be one non-missing value.", call. = FALSE)
    }
    data <- data[as.character(data[[experiment_col]]) == as.character(experiment), , drop = FALSE]
  }
  if (!nrow(data)) stop("No orthogonal rows remain after filtering.", call. = FALSE)
  identifier_columns <- c(construct_columns, replicate, perturbation)
  if (any(vapply(data[identifier_columns], function(x) anyNA(x) || any(as.character(x) == ""),
                 logical(1)))) {
    stop("Orthogonal construct, replicate, and perturbation identifiers must be complete.",
         call. = FALSE)
  }
  if (is.null(targets)) targets <- sort(unique(as.character(effect_table$target)))
  targets <- unique(as.character(targets))
  if (!length(targets) || anyNA(targets) || any(targets == "") || control %in% targets) {
    stop("`targets` must contain non-control perturbation labels.", call. = FALSE)
  }
  data <- data[as.character(data[[perturbation]]) %in% c(control, targets), , drop = FALSE]
  group_columns <- c(construct_columns, replicate, perturbation)
  aggregated <- .orthogonal_usage(data, group_columns, proximal, distal, total)

  control_rows <- aggregated[as.character(aggregated[[perturbation]]) == control, , drop = FALSE]
  target_rows <- aggregated[as.character(aggregated[[perturbation]]) %in% targets, , drop = FALSE]
  match_columns <- c(construct_columns, replicate)
  control_keep <- control_rows[c(match_columns, "proximal_usage", "total_count")]
  names(control_keep)[names(control_keep) == "proximal_usage"] <- "control_proximal_usage"
  names(control_keep)[names(control_keep) == "total_count"] <- "control_total"
  paired <- merge(target_rows, control_keep, by = match_columns, all.x = TRUE)
  paired$orthogonal_delta <- paired$proximal_usage - paired$control_proximal_usage

  construct_key_columns <- c(construct_columns, perturbation)
  all_replicates <- sort(unique(as.character(data[[replicate]])))
  groups <- interaction(paired[construct_key_columns], drop = TRUE, lex.order = TRUE)
  constructs <- do.call(rbind, lapply(split(paired, groups), function(x) {
    identity <- x[1, construct_key_columns, drop = FALSE]
    paired_rows <- x[is.finite(x$orthogonal_delta), , drop = FALSE]
    replicate_values <- sort(unique(as.character(paired_rows[[replicate]])))
    deltas <- vapply(replicate_values, function(rep_value) {
      mean(paired_rows$orthogonal_delta[
        as.character(paired_rows[[replicate]]) == rep_value
      ])
    }, numeric(1))
    out <- cbind(
      identity,
      orthogonal_replicate_n = length(replicate_values),
      orthogonal_proximal_delta = if (length(deltas)) mean(deltas) else NA_real_,
      orthogonal_total_perturbation = sum(paired_rows$total_count),
      orthogonal_total_control = sum(paired_rows$control_total),
      orthogonal_replicate_sign_concordant = length(deltas) >= min_replicates &&
        length(unique(sign(deltas))) == 1L,
      stringsAsFactors = FALSE
    )
    for (rep_value in all_replicates) {
      out[[paste0("orthogonal_delta_", make.names(rep_value))]] <-
        if (rep_value %in% names(deltas)) deltas[[rep_value]] else NA_real_
    }
    out
  }))
  names(constructs)[names(constructs) == perturbation] <- "target"
  names(constructs)[names(constructs) == construct_columns[[1]]] <- "gene"
  mapped <- merge(constructs, effect_table, by = c("target", "gene"), all.x = TRUE,
                  suffixes = c("", "_effect"))
  mapped$mapping_status <- ifelse(
    mapped$orthogonal_replicate_n < min_replicates,
    "ORTHOGONAL_REPLICATES_INSUFFICIENT",
    ifelse(is.na(mapped$proximal_direction), "NOT_IN_FROZEN_EFFECT_SET",
           ifelse(is.na(mapped$tested) | !mapped$tested, "FROZEN_EFFECT_UNTESTED", "MAPPED"))
  )
  mapped$direction_concordant <- mapped$mapping_status == "MAPPED" &
    sign(mapped$orthogonal_proximal_delta) == sign(mapped$proximal_direction)
  mapped$direction_concordant[mapped$mapping_status != "MAPPED"] <- NA
  mapped <- mapped[order(mapped$target, mapped$gene), , drop = FALSE]
  rownames(mapped) <- NULL

  mapped_only <- mapped[mapped$mapping_status == "MAPPED", , drop = FALSE]
  target_summary <- do.call(rbind, lapply(split(mapped_only, mapped_only$target), function(x) {
    concordant_replicates <- x$orthogonal_replicate_sign_concordant
    weights <- x$orthogonal_total_perturbation + x$orthogonal_total_control
    data.frame(
      target = x$target[[1]], construct_n = nrow(x), gene_n = length(unique(x$gene)),
      replicate_concordant_construct_n = sum(concordant_replicates),
      direction_concordance_all = mean(x$direction_concordant),
      direction_concordance_orthogonal_replicate_concordant =
        if (any(concordant_replicates)) mean(x$direction_concordant[concordant_replicates]) else NA_real_,
      weighted_direction_concordance = stats::weighted.mean(x$direction_concordant, weights),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(target_summary)) target_summary <- data.frame()
  accounting <- do.call(rbind, lapply(targets, function(target_value) {
    x <- mapped[mapped$target == target_value, , drop = FALSE]
    if (!nrow(x)) {
      return(data.frame(
        target = target_value, target_status = "ORTHOGONAL_TARGET_NOT_PRESENT",
        construct_n = 0L, mapped_n = 0L, not_in_effect_set_n = 0L,
        effect_untested_n = 0L, orthogonal_unsupported_n = 0L,
        stringsAsFactors = FALSE
      ))
    }
    counts <- table(factor(x$mapping_status, levels = c(
      "MAPPED", "NOT_IN_FROZEN_EFFECT_SET", "FROZEN_EFFECT_UNTESTED",
      "ORTHOGONAL_REPLICATES_INSUFFICIENT"
    )))
    data.frame(target = target_value,
               target_status = if (unname(counts[["MAPPED"]]) > 0L)
                 "SUPPORTED_WITH_MAPPED_CONSTRUCTS" else "NO_MAPPED_CONSTRUCTS",
               construct_n = nrow(x),
               mapped_n = unname(counts[["MAPPED"]]),
               not_in_effect_set_n = unname(counts[["NOT_IN_FROZEN_EFFECT_SET"]]),
               effect_untested_n = unname(counts[["FROZEN_EFFECT_UNTESTED"]]),
               orthogonal_unsupported_n = unname(counts[["ORTHOGONAL_REPLICATES_INSUFFICIENT"]]),
               stringsAsFactors = FALSE)
  }))
  output <- list(results = mapped, target_summary = target_summary,
                 accounting = accounting, paired_replicates = paired,
                 settings = list(freeze_fingerprint = freeze$fingerprint,
                                 orthogonal_registry_id = orthogonal_registry_id,
                                 control = control,
                                 min_replicates = min_replicates,
                                 experiment = experiment))
  class(output) <- "repliapa_orthogonal_validation"
  output
}

#' @export
print.repliapa_orthogonal_validation <- function(x, ...) {
  cat("RepliAPA post-freeze orthogonal validation\n")
  cat("  mapped / total constructs: ", sum(x$accounting$mapped_n), " / ",
      sum(x$accounting$construct_n), "\n", sep = "")
  if (nrow(x$target_summary)) {
    cat("  targets: ", paste(x$target_summary$target, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}
