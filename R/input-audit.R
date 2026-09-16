.audit_positive_count <- function(value, name, allow_zero = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < if (allow_zero) 0 else 1 || abs(value - round(value)) > 1e-8) {
    qualifier <- if (allow_zero) "a nonnegative" else "a positive"
    stop("`", name, "` must be ", qualifier, " integer.", call. = FALSE)
  }
  as.integer(value)
}

#' Audit a guide-by-replicate PAS pseudobulk contract
#'
#' Performs the pre-fit data-contract and adequacy audit used by RepliAPA. The
#' function does not alter counts, choose favorable targets, or fit an outcome
#' model. PAS support tiers are computed only from predeclared development
#' replicates and qualified sample units.
#'
#' @param counts PAS-by-sample-unit pseudobulk count matrix.
#' @param sample_data Sample metadata aligned to columns of `counts`.
#' @param pas_data PAS metadata aligned to rows of `counts`.
#' @param development_replicates Replicate labels allowed to define PAS support
#'   tiers. Validation or external replicates must not be supplied here.
#' @param gene,replicate,guide,target,ntc,cell_n Column names.
#' @param min_cells,min_polya_count Predeclared sample-unit qualification rules.
#' @param min_target_guides Minimum independent qualified guides required in
#'   every biological replicate for primary target eligibility.
#' @param usable_total,usable_nonzero PAS-A2 lower bounds.
#' @param strong_total,strong_nonzero PAS-A3 lower bounds.
#' @param seqname,start,end Optional PAS coordinate column names. If all are
#'   present, exact duplicate coordinates within a gene are marked as conflicts.
#' @return A `repliapa_input_audit` containing contract checks, sample-unit,
#'   guide, target, PAS, and gene adequacy tables plus immutable settings.
#' @export
repliapa_audit_input <- function(
    counts, sample_data, pas_data, development_replicates,
    gene = "gene", replicate = "replicate", guide = "guide",
    target = "target", ntc = "ntc", cell_n = "cell_n",
    min_cells = 20L, min_polya_count = 150000L,
    min_target_guides = 2L,
    usable_total = 20L, usable_nonzero = 2L,
    strong_total = 200L, strong_nonzero = 10L,
    seqname = "seqnames", start = "start", end = "end") {
  repliapa_validate_pseudobulk(counts, sample_data, replicate, guide, target, ntc)
  .aligned_data(pas_data, rownames(counts), "pas_data")
  .require_columns(pas_data, gene, "pas_data")
  .require_columns(sample_data, cell_n, "sample_data")

  min_cells <- .audit_positive_count(min_cells, "min_cells")
  min_polya_count <- .audit_positive_count(min_polya_count, "min_polya_count")
  min_target_guides <- .audit_positive_count(min_target_guides, "min_target_guides")
  usable_total <- .audit_positive_count(usable_total, "usable_total")
  usable_nonzero <- .audit_positive_count(usable_nonzero, "usable_nonzero")
  strong_total <- .audit_positive_count(strong_total, "strong_total")
  strong_nonzero <- .audit_positive_count(strong_nonzero, "strong_nonzero")
  if (strong_total < usable_total || strong_nonzero < usable_nonzero) {
    stop("Strong PAS thresholds must not be below usable thresholds.", call. = FALSE)
  }

  development_replicates <- sort(unique(as.character(development_replicates)))
  available_replicates <- sort(unique(as.character(sample_data[[replicate]])))
  if (!length(development_replicates) || anyNA(development_replicates) ||
      any(development_replicates == "") ||
      any(!development_replicates %in% available_replicates)) {
    stop("`development_replicates` must be non-missing available replicate labels.",
         call. = FALSE)
  }
  gene_value <- as.character(pas_data[[gene]])
  if (anyNA(gene_value) || any(gene_value == "")) {
    stop("The PAS-to-gene map is incomplete.", call. = FALSE)
  }
  cell_value <- sample_data[[cell_n]]
  if (!is.numeric(cell_value) || anyNA(cell_value) || any(!is.finite(cell_value)) ||
      any(cell_value < 0) || any(abs(cell_value - round(cell_value)) > 1e-8)) {
    stop("The sample-unit cell count must contain nonnegative integers.", call. = FALSE)
  }
  guide_value <- as.character(sample_data[[guide]])
  target_value <- as.character(sample_data[[target]])
  ntc_value <- sample_data[[ntc]]
  guide_target_n <- vapply(split(target_value, guide_value),
                           function(x) length(unique(x)), integer(1))
  guide_ntc_n <- vapply(split(ntc_value, guide_value),
                        function(x) length(unique(x)), integer(1))
  if (any(guide_target_n != 1L) || any(guide_ntc_n != 1L)) {
    stop("A guide maps to multiple targets or NTC states across replicates.", call. = FALSE)
  }

  polya_count <- as.numeric(Matrix::colSums(counts))
  qualified <- cell_value >= min_cells & polya_count >= min_polya_count
  sample_units <- data.frame(
    sample_unit = rownames(sample_data),
    replicate = as.character(sample_data[[replicate]]),
    guide = guide_value,
    target = target_value,
    ntc = ntc_value,
    cell_n = as.integer(cell_value),
    polya_count = polya_count,
    qualified = qualified,
    qualification_reason = ifelse(
      qualified, "NONE",
      ifelse(cell_value < min_cells & polya_count < min_polya_count,
             "LOW_CELL_N_AND_POLYA_COUNT",
             ifelse(cell_value < min_cells, "LOW_CELL_N", "LOW_POLYA_COUNT"))
    ),
    stringsAsFactors = FALSE
  )

  guide_keys <- sort(unique(guide_value))
  guide_support <- do.call(rbind, lapply(guide_keys, function(guide_id) {
    rows <- guide_value == guide_id
    data.frame(
      guide = guide_id,
      target = unique(target_value[rows]),
      ntc = unique(ntc_value[rows]),
      replicate_n = length(unique(as.character(sample_data[[replicate]][rows]))),
      qualified_replicate_n = length(unique(sample_units$replicate[rows & qualified])),
      sample_unit_n = sum(rows),
      qualified_sample_unit_n = sum(rows & qualified),
      cell_n = sum(cell_value[rows]),
      polya_count = sum(polya_count[rows]),
      stringsAsFactors = FALSE
    )
  }))

  targets <- sort(unique(target_value))
  target_replicate_support <- do.call(rbind, lapply(targets, function(target_id) {
    do.call(rbind, lapply(available_replicates, function(rep_id) {
      rows <- target_value == target_id &
        as.character(sample_data[[replicate]]) == rep_id & qualified
      data.frame(
        target = target_id, replicate = rep_id,
        qualified_guide_n = length(unique(guide_value[rows])),
        qualified_sample_unit_n = sum(rows),
        stringsAsFactors = FALSE
      )
    }))
  }))
  target_support <- do.call(rbind, lapply(targets, function(target_id) {
    is_ntc <- all(ntc_value[target_value == target_id])
    rep_rows <- target_replicate_support$target == target_id
    replicate_counts <- target_replicate_support$qualified_guide_n[rep_rows]
    eligible <- !is_ntc && length(replicate_counts) == length(available_replicates) &&
      all(replicate_counts >= min_target_guides)
    insufficient <- target_replicate_support$replicate[rep_rows][
      replicate_counts < min_target_guides
    ]
    data.frame(
      target = target_id,
      ntc = is_ntc,
      guide_n_total = length(unique(guide_value[target_value == target_id])),
      qualified_guide_n = length(unique(guide_value[target_value == target_id & qualified])),
      replicate_n = length(available_replicates),
      min_qualified_guides_per_replicate = min(replicate_counts),
      target_eligible = eligible,
      exclusion_reason = if (is_ntc) "NEGATIVE_CONTROL" else if (eligible) "NONE" else
        paste0("INSUFFICIENT_GUIDES:", paste(insufficient, collapse = ";")),
      stringsAsFactors = FALSE
    )
  }))

  development_rows <- as.character(sample_data[[replicate]]) %in%
    development_replicates & qualified
  development_counts <- counts[, development_rows, drop = FALSE]
  pas_total <- as.numeric(Matrix::rowSums(development_counts))
  pas_nonzero <- as.integer(Matrix::rowSums(development_counts > 0))
  tier <- ifelse(
    pas_total == 0, "PAS-A0",
    ifelse(pas_total < usable_total | pas_nonzero < usable_nonzero, "PAS-A1",
           ifelse(pas_total >= strong_total & pas_nonzero >= strong_nonzero,
                  "PAS-A3", "PAS-A2"))
  )
  coordinate_conflict <- rep(FALSE, nrow(pas_data))
  coordinate_columns <- c(seqname, start, end)
  if (all(coordinate_columns %in% names(pas_data))) {
    coordinate_key <- paste(gene_value, pas_data[[seqname]], pas_data[[start]],
                            pas_data[[end]], sep = "\r")
    coordinate_conflict <- duplicated(coordinate_key) | duplicated(coordinate_key, fromLast = TRUE)
  }
  pas_adequacy <- data.frame(
    pas_id = rownames(pas_data),
    gene = gene_value,
    development_total_count = pas_total,
    development_nonzero_sample_units = pas_nonzero,
    all_sample_total_count = as.numeric(Matrix::rowSums(counts)),
    mapping_unique = TRUE,
    exact_coordinate_conflict = coordinate_conflict,
    adequacy_tier = tier,
    primary_eligible = tier %in% c("PAS-A2", "PAS-A3") & !coordinate_conflict,
    stringsAsFactors = FALSE
  )
  pas_adequacy <- pas_adequacy[order(pas_adequacy$gene, pas_adequacy$pas_id), , drop = FALSE]

  gene_indices <- split(seq_len(nrow(pas_adequacy)), pas_adequacy$gene)
  gene_adequacy <- do.call(rbind, lapply(names(gene_indices), function(gene_id) {
    rows <- gene_indices[[gene_id]]
    tiers <- pas_adequacy$adequacy_tier[rows]
    usable_n <- sum(pas_adequacy$primary_eligible[rows])
    conflict <- any(pas_adequacy$exact_coordinate_conflict[rows])
    eligible <- length(rows) >= 2L && usable_n >= 2L && !conflict
    reason <- if (eligible) "NONE" else if (conflict) "EXACT_COORDINATE_CONFLICT" else
      if (length(rows) < 2L) "SINGLE_PAS_GENE" else "FEWER_THAN_TWO_USABLE_PAS"
    data.frame(
      gene = gene_id,
      pas_n_total = length(rows),
      pas_n_A0 = sum(tiers == "PAS-A0"),
      pas_n_A1 = sum(tiers == "PAS-A1"),
      pas_n_A2 = sum(tiers == "PAS-A2"),
      pas_n_A3 = sum(tiers == "PAS-A3"),
      usable_pas_n = usable_n,
      development_total_count = sum(pas_adequacy$development_total_count[rows]),
      coordinate_conflict = conflict,
      primary_multiPAS_eligible = eligible,
      exclusion_reason = reason,
      stringsAsFactors = FALSE
    )
  }))

  contract_checks <- data.frame(
    check = c(
      "NONNEGATIVE_INTEGER_COUNTS", "SAMPLE_ALIGNMENT", "PAS_ALIGNMENT",
      "UNIQUE_PAS_IDS", "UNIQUE_SAMPLE_UNITS", "GUIDE_TARGET_MAPPING",
      "GUIDE_NTC_MAPPING", "NO_CELL_AS_REPLICATE", "DEVELOPMENT_ROLE_DECLARED"
    ),
    passed = TRUE,
    detail = c(
      nrow(counts), ncol(counts), nrow(pas_data), nrow(counts), ncol(counts),
      length(guide_keys), length(guide_keys), nrow(sample_units),
      paste(development_replicates, collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
  status <- if (any(!qualified) || any(!target_support$target_eligible & !target_support$ntc) ||
                any(!gene_adequacy$primary_multiPAS_eligible)) {
    "PASS_WITH_EXCLUSIONS"
  } else "PASS"
  summary <- data.frame(
    data_contract_status = status,
    sample_unit_n = nrow(sample_units),
    qualified_sample_unit_n = sum(qualified),
    replicate_n = length(available_replicates),
    guide_n = length(guide_keys),
    target_n_non_ntc = sum(!target_support$ntc),
    eligible_target_n = sum(target_support$target_eligible),
    pas_n = nrow(pas_adequacy),
    primary_pas_n = sum(pas_adequacy$primary_eligible),
    gene_n = nrow(gene_adequacy),
    primary_multiPAS_gene_n = sum(gene_adequacy$primary_multiPAS_eligible),
    stringsAsFactors = FALSE
  )
  settings <- list(
    development_replicates = development_replicates,
    min_cells = min_cells, min_polya_count = min_polya_count,
    min_target_guides = min_target_guides,
    usable_total = usable_total, usable_nonzero = usable_nonzero,
    strong_total = strong_total, strong_nonzero = strong_nonzero,
    outcome_used = FALSE
  )
  output <- list(
    summary = summary, contract_checks = contract_checks,
    sample_units = sample_units, guide_support = guide_support,
    target_replicate_support = target_replicate_support,
    target_support = target_support, pas_adequacy = pas_adequacy,
    gene_adequacy = gene_adequacy, settings = settings
  )
  class(output) <- "repliapa_input_audit"
  output
}

#' @export
print.repliapa_input_audit <- function(x, ...) {
  cat("RepliAPA input audit\n")
  cat("  contract: ", x$summary$data_contract_status, "\n", sep = "")
  cat("  sample units qualified: ", x$summary$qualified_sample_unit_n, " / ",
      x$summary$sample_unit_n, "\n", sep = "")
  cat("  eligible targets: ", x$summary$eligible_target_n, " / ",
      x$summary$target_n_non_ntc, "\n", sep = "")
  cat("  primary multi-PAS genes: ", x$summary$primary_multiPAS_gene_n, " / ",
      x$summary$gene_n, "\n", sep = "")
  cat("  outcome used: FALSE\n")
  invisible(x)
}
