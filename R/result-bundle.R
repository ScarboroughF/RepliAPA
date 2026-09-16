.empty_bundle_tables <- function() {
  list(
    target_effects = data.frame(
      target = character(), gene = character(), pas_n = integer(),
      total_count = numeric(), target_guide_n = integer(),
      ntc_guide_n = integer(), permutation_n = integer(),
      statistic = numeric(), p_value = numeric(), effect_norm = numeric(),
      coordinate_shift = numeric(), tested = logical(), filtered = logical(),
      abstained = logical(), reason = character(), q_value = numeric(),
      significant = logical(),
      stringsAsFactors = FALSE
    ),
    pas_effects = data.frame(
      target = character(), gene = character(), pas_id = character(),
      effect = numeric(), stringsAsFactors = FALSE
    ),
    guide_effects = data.frame(
      target = character(), gene = character(), pas_id = character(),
      replicate = character(), guide = character(), effect = numeric(),
      stringsAsFactors = FALSE
    ),
    target_heterogeneity = data.frame(
      target = character(), gene = character(), replicate = character(),
      guide_n = integer(), mean_pairwise_cosine = numeric(),
      effect_dispersion = numeric(), stringsAsFactors = FALSE
    )
  )
}

.bind_fit_table <- function(fits, component) {
  values <- lapply(names(fits), function(target_value) {
    table <- fits[[target_value]][[component]]
    if (is.null(table) || !nrow(table)) return(NULL)
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    table$target <- target_value
    table[c("target", setdiff(names(table), "target"))]
  })
  values <- Filter(Negate(is.null), values)
  if (!length(values)) return(NULL)
  output <- do.call(rbind, values)
  rownames(output) <- NULL
  output
}

.sort_bundle_table <- function(table, keys) {
  if (!nrow(table)) return(table)
  keys <- intersect(keys, names(table))
  if (!length(keys)) return(table)
  key_values <- unname(lapply(table[keys], function(x) as.character(x)))
  order_value <- do.call(order, c(key_values, list(na.last = TRUE, method = "radix")))
  output <- table[order_value, , drop = FALSE]
  rownames(output) <- NULL
  output
}

#' Construct a canonical auditable result bundle
#'
#' Standardizes fitted results without rerunning inference. Every tested,
#' filtered, abstained, and unsupported gene or target remains in an accounting
#' table. Optional input-audit tables and freeze provenance are included without
#' changing any result.
#'
#' @param x A `repliapa_fit` or `repliapa_fit_list` object.
#' @param freeze Optional verified `repliapa_freeze` object.
#' @param input_audit Optional `repliapa_input_audit` object.
#' @param compute_benchmark Optional externally instrumented benchmark table.
#'   It must report unique methods, repetitions, median elapsed seconds, median
#'   peak RSS, input SHA256, and the elapsed/memory measurement scopes. RepliAPA
#'   validates and records these measurements but does not infer missing values.
#' @return A deterministic `repliapa_result_bundle` containing named data-frame
#'   tables and provenance. Runtime and peak memory are not invented; absent
#'   instrumentation is explicitly recorded.
#' @export
repliapa_result_bundle <- function(x, freeze = NULL, input_audit = NULL,
                                   compute_benchmark = NULL) {
  if (inherits(x, "repliapa_fit")) {
    fits <- stats::setNames(list(x), x$target)
    accounting <- data.frame(
      target = x$target, status = "FITTED", reason = "NONE",
      target_guide_n = length(x$target_guides),
      ntc_guide_n = length(x$ntc_guides),
      tested_gene_n = sum(x$gene_table$tested),
      filtered_gene_n = sum(x$gene_table$filtered),
      abstained_gene_n = sum(x$gene_table$abstained),
      discovery_n = sum(x$gene_table$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    requested_targets <- x$target
  } else if (inherits(x, "repliapa_fit_list")) {
    fits <- x$fits
    accounting <- as.data.frame(x$accounting, stringsAsFactors = FALSE)
    requested_targets <- as.character(x$targets)
  } else {
    stop("`x` must be a repliapa_fit or repliapa_fit_list object.", call. = FALSE)
  }
  if (anyDuplicated(names(fits)) || anyDuplicated(as.character(accounting$target))) {
    stop("Target identifiers in fitted results must be unique.", call. = FALSE)
  }

  tables <- .empty_bundle_tables()
  collected <- list(
    target_effects = .bind_fit_table(fits, "gene_table"),
    pas_effects = .bind_fit_table(fits, "pas_effects"),
    guide_effects = .bind_fit_table(fits, "guide_effects"),
    target_heterogeneity = .bind_fit_table(fits, "heterogeneity")
  )
  for (name in names(collected)) {
    if (!is.null(collected[[name]])) tables[[name]] <- collected[[name]]
  }
  tables$target_accounting <- accounting

  target_summary <- do.call(rbind, lapply(requested_targets, function(target_value) {
    row <- accounting[as.character(accounting$target) == target_value, , drop = FALSE]
    fit <- fits[[target_value]]
    if (is.null(fit)) {
      return(data.frame(
        target = target_value, status = "UNSUPPORTED",
        native_gene_n = 0L, tested_gene_n = 0L, filtered_gene_n = 0L,
        abstained_gene_n = 0L, discovery_n = 0L,
        tested_fraction = NA_real_, discovery_fraction_of_tested = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    native_n <- nrow(fit$gene_table)
    tested_n <- sum(fit$gene_table$tested)
    discovery_n <- sum(fit$gene_table$significant, na.rm = TRUE)
    data.frame(
      target = target_value,
      status = if (nrow(row)) as.character(row$status[[1]]) else "FITTED",
      native_gene_n = native_n,
      tested_gene_n = tested_n,
      filtered_gene_n = sum(fit$gene_table$filtered),
      abstained_gene_n = sum(fit$gene_table$abstained),
      discovery_n = discovery_n,
      tested_fraction = if (native_n) tested_n / native_n else NA_real_,
      discovery_fraction_of_tested = if (tested_n) discovery_n / tested_n else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  tables$tested_gene_accounting <- target_summary

  audit_included <- !is.null(input_audit)
  if (audit_included) {
    if (!inherits(input_audit, "repliapa_input_audit")) {
      stop("`input_audit` must be a repliapa_input_audit object.", call. = FALSE)
    }
    tables$input_audit_summary <- as.data.frame(input_audit$summary)
    tables$input_contract_checks <- as.data.frame(input_audit$contract_checks)
    tables$sample_unit_manifest <- as.data.frame(input_audit$sample_units)
    tables$guide_support <- as.data.frame(input_audit$guide_support)
    tables$target_guide_support <- as.data.frame(input_audit$target_support)
    tables$pas_adequacy <- as.data.frame(input_audit$pas_adequacy)
    tables$gene_pas_adequacy <- as.data.frame(input_audit$gene_adequacy)
  }

  benchmark_included <- !is.null(compute_benchmark)
  if (benchmark_included) {
    if (!is.data.frame(compute_benchmark)) {
      stop("`compute_benchmark` must be a data frame.", call. = FALSE)
    }
    required_benchmark <- c(
      "method", "repetitions", "elapsed_seconds_median", "peak_rss_mb_median",
      "input_sha256", "elapsed_scope", "peak_rss_scope", "peak_rss_source"
    )
    .require_columns(compute_benchmark, required_benchmark, "compute_benchmark")
    benchmark <- as.data.frame(compute_benchmark, stringsAsFactors = FALSE)
    method <- as.character(benchmark$method)
    if (!nrow(benchmark) || anyNA(method) || any(method == "") || anyDuplicated(method)) {
      stop("Compute benchmark methods must be non-missing and unique.", call. = FALSE)
    }
    numeric_positive <- c(
      "repetitions", "elapsed_seconds_median", "peak_rss_mb_median"
    )
    if (any(vapply(numeric_positive, function(column) {
      value <- benchmark[[column]]
      !is.numeric(value) || anyNA(value) || any(!is.finite(value)) || any(value <= 0)
    }, logical(1)))) {
      stop("Compute benchmark repetitions, elapsed time, and peak RSS must be positive numbers.",
           call. = FALSE)
    }
    if (any(abs(benchmark$repetitions - round(benchmark$repetitions)) > 1e-8)) {
      stop("Compute benchmark repetitions must be integers.", call. = FALSE)
    }
    character_required <- c(
      "input_sha256", "elapsed_scope", "peak_rss_scope", "peak_rss_source"
    )
    if (any(vapply(character_required, function(column) {
      value <- as.character(benchmark[[column]])
      anyNA(value) || any(value == "")
    }, logical(1)))) {
      stop("Compute benchmark provenance fields must be complete.", call. = FALSE)
    }
    if (any(!grepl("^[0-9a-f]{64}$", as.character(benchmark$input_sha256)))) {
      stop("Compute benchmark input SHA256 values are invalid.", call. = FALSE)
    }
    tables$compute_benchmark <- .sort_bundle_table(benchmark, "method")
  }

  freeze_fingerprint <- "NOT_PROVIDED"
  if (!is.null(freeze)) {
    repliapa_verify_freeze(freeze)
    freeze_fingerprint <- freeze$fingerprint
  }
  package_version <- tryCatch(
    as.character(utils::packageVersion("RepliAPA")),
    error = function(e) "DEVELOPMENT_SOURCE"
  )
  provenance <- data.frame(
    bundle_schema_version = 1L,
    package_version = package_version,
    target_n_requested = length(requested_targets),
    target_n_fitted = length(fits),
    freeze_fingerprint = freeze_fingerprint,
    input_audit_included = audit_included,
    compute_benchmark_included = benchmark_included,
    inference_changed = FALSE,
    runtime_status = if (benchmark_included) {
      "MEASURED_ISOLATED_PROCESS_MEDIAN"
    } else "NOT_INSTRUMENTED",
    peak_memory_status = if (benchmark_included) {
      "MEASURED_KERNEL_VM_HWM"
    } else "NOT_INSTRUMENTED",
    stringsAsFactors = FALSE
  )
  tables$provenance <- provenance

  sort_keys <- list(
    target_effects = c("target", "gene"),
    pas_effects = c("target", "gene", "pas_id"),
    guide_effects = c("target", "replicate", "guide", "gene", "pas_id"),
    target_heterogeneity = c("target", "replicate", "gene"),
    target_accounting = "target",
    tested_gene_accounting = "target",
    sample_unit_manifest = c("replicate", "target", "guide", "sample_unit"),
    guide_support = c("target", "guide"),
    target_guide_support = "target",
    pas_adequacy = c("gene", "pas_id"),
    gene_pas_adequacy = "gene",
    input_contract_checks = "check",
    compute_benchmark = "method"
  )
  for (name in intersect(names(sort_keys), names(tables))) {
    tables[[name]] <- .sort_bundle_table(tables[[name]], sort_keys[[name]])
  }
  output <- list(tables = tables, schema_version = 1L)
  class(output) <- "repliapa_result_bundle"
  output
}

#' @export
print.repliapa_result_bundle <- function(x, ...) {
  accounting <- x$tables$tested_gene_accounting
  cat("RepliAPA result bundle\n")
  cat("  fitted / requested targets: ", sum(accounting$status == "FITTED"), " / ",
      nrow(accounting), "\n", sep = "")
  cat("  tested genes: ", sum(accounting$tested_gene_n), "\n", sep = "")
  cat("  discoveries: ", sum(accounting$discovery_n), "\n", sep = "")
  provenance <- x$tables$provenance
  cat("  runtime / peak memory: ", provenance$runtime_status, " / ",
      provenance$peak_memory_status, "\n", sep = "")
  invisible(x)
}

.bundle_file_names <- function(table_names) {
  stats::setNames(paste0(toupper(table_names), ".tsv"), table_names)
}

#' Write an auditable result bundle
#'
#' Writes deterministic tab-separated tables, a table manifest, and a SHA256
#' checksum file. Existing output directories fail closed unless `overwrite` is
#' explicitly enabled.
#'
#' @param bundle A `repliapa_result_bundle` object.
#' @param path Output directory.
#' @param overwrite Whether to replace an existing non-empty output directory.
#' @return The normalized output path, invisibly.
#' @export
repliapa_write_bundle <- function(bundle, path, overwrite = FALSE) {
  if (!inherits(bundle, "repliapa_result_bundle")) {
    stop("`bundle` must be a repliapa_result_bundle object.", call. = FALSE)
  }
  if (length(path) != 1L || is.na(path) || path == "") {
    stop("`path` must be one non-missing directory path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  path <- path.expand(path)
  if (dir.exists(path)) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) && !overwrite) {
      stop("Output directory is not empty; set `overwrite = TRUE` explicitly.",
           call. = FALSE)
    }
    if (length(existing) && overwrite) {
      expected <- c(unname(.bundle_file_names(names(bundle$tables))),
                    "BUNDLE_MANIFEST.tsv", "CHECKSUMS.sha256")
      unexpected <- setdiff(existing, expected)
      if (length(unexpected)) {
        stop("Refusing to overwrite a directory containing unrelated files.", call. = FALSE)
      }
      unlink(file.path(path, existing), recursive = FALSE, force = TRUE)
    }
  } else if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create the output directory.", call. = FALSE)
  }

  file_names <- .bundle_file_names(names(bundle$tables))
  manifest_rows <- lapply(names(bundle$tables), function(name) {
    table <- bundle$tables[[name]]
    if (!is.data.frame(table)) stop("Every bundle table must be a data frame.", call. = FALSE)
    destination <- file.path(path, file_names[[name]])
    utils::write.table(table, destination, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE, na = "")
    data.frame(
      table = name,
      file = file_names[[name]],
      row_n = nrow(table),
      column_n = ncol(table),
      sha256 = digest::digest(destination, file = TRUE, algo = "sha256",
                              serialize = FALSE),
      stringsAsFactors = FALSE
    )
  })
  manifest <- do.call(rbind, manifest_rows)
  manifest <- manifest[order(manifest$table, method = "radix"), , drop = FALSE]
  manifest_path <- file.path(path, "BUNDLE_MANIFEST.tsv")
  utils::write.table(manifest, manifest_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "")

  checksum_files <- c(manifest$file, "BUNDLE_MANIFEST.tsv")
  checksum_files <- sort(checksum_files, method = "radix")
  hashes <- vapply(file.path(path, checksum_files), digest::digest, character(1),
                   file = TRUE, algo = "sha256", serialize = FALSE)
  writeLines(paste(hashes, checksum_files, sep = "  "),
             file.path(path, "CHECKSUMS.sha256"), useBytes = TRUE)
  repliapa_verify_bundle(path)
  invisible(normalizePath(path, mustWork = TRUE))
}

#' Verify an auditable result bundle
#'
#' @param path Result-bundle directory created by [repliapa_write_bundle()].
#' @return `TRUE`; otherwise stops when a file is missing, unsafe, duplicated,
#'   or differs from its recorded SHA256.
#' @export
repliapa_verify_bundle <- function(path) {
  if (length(path) != 1L || is.na(path) || !dir.exists(path)) {
    stop("`path` must be an existing result-bundle directory.", call. = FALSE)
  }
  checksum_path <- file.path(path, "CHECKSUMS.sha256")
  if (!file.exists(checksum_path)) stop("Bundle checksum file is missing.", call. = FALSE)
  lines <- readLines(checksum_path, warn = FALSE)
  matched <- regexec("^([0-9a-f]{64})  ([^/\\\\]+)$", lines)
  pieces <- regmatches(lines, matched)
  if (!length(lines) || any(lengths(pieces) != 3L)) {
    stop("Bundle checksum file has an invalid or unsafe entry.", call. = FALSE)
  }
  hashes <- vapply(pieces, `[[`, character(1), 2L)
  files <- vapply(pieces, `[[`, character(1), 3L)
  if (anyDuplicated(files) || !"BUNDLE_MANIFEST.tsv" %in% files) {
    stop("Bundle checksum entries are duplicated or incomplete.", call. = FALSE)
  }
  paths <- file.path(path, files)
  if (any(!file.exists(paths))) stop("A recorded bundle file is missing.", call. = FALSE)
  observed <- vapply(paths, digest::digest, character(1), file = TRUE,
                     algo = "sha256", serialize = FALSE)
  if (!identical(unname(observed), unname(hashes))) {
    stop("Bundle SHA256 verification failed.", call. = FALSE)
  }
  TRUE
}
