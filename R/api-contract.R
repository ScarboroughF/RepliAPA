.canonical_formals <- function(fun) {
  values <- vapply(formals(fun), function(value) {
    if (identical(value, quote(expr = ))) return("<required>")
    paste(deparse(value, width.cutoff = 500L), collapse = " ")
  }, character(1))
  paste(paste(names(values), values, sep = "="), collapse = "|")
}

.contract_table <- function(required_columns, required = TRUE,
                            optional_columns = character(),
                            extension_patterns = character(),
                            extension_type = NULL) {
  list(
    required = required,
    required_columns = required_columns,
    optional_columns = optional_columns,
    extension_patterns = extension_patterns,
    extension_type = extension_type,
    column_types = character()
  )
}

.contract_types <- function(character = character(), integer = character(),
                            double = character(), numeric = character(),
                            logical = character()) {
  output <- c(
    stats::setNames(rep("character", length(character)), character),
    stats::setNames(rep("integer", length(integer)), integer),
    stats::setNames(rep("double", length(double)), double),
    stats::setNames(rep("numeric", length(numeric)), numeric),
    stats::setNames(rep("logical", length(logical)), logical)
  )
  output
}

.repliapa_contract_payload <- function() {
  functions <- c(
    repliapa_api_contract = "",
    repliapa_apply_external = 'primary=<required>|external=<required>|freeze=<required>|external_registry_id=<required>|primary_label="primary"|external_label="external"|targets=NULL|fdr=0.05',
    repliapa_audit_input = 'counts=<required>|sample_data=<required>|pas_data=<required>|development_replicates=<required>|gene="gene"|replicate="replicate"|guide="guide"|target="target"|ntc="ntc"|cell_n="cell_n"|min_cells=20L|min_polya_count=150000L|min_target_guides=2L|usable_total=20L|usable_nonzero=2L|strong_total=200L|strong_nonzero=10L|seqname="seqnames"|start="start"|end="end"',
    repliapa_authorize = "freeze=<required>|registry_id=<required>|action=<required>",
    repliapa_classify_effects = "fit=<required>|gene_expression=<required>|apa_fdr=0.05|ge_abs_log2fc=0.25|ge_significant=NULL",
    repliapa_collapse_proximal_distal = 'counts=<required>|pas_data=<required>|gene="gene"|position="position"|strand="strand"',
    repliapa_collect_guide_effects = "fits=<required>",
    repliapa_compare_allpas_pd = "all_pas=<required>|proximal_distal=<required>|fdr=0.05",
    repliapa_compare_baseline = 'repliapa=<required>|baseline=<required>|baseline_gene="gene_id"|baseline_p="pvalue"|baseline_q=NULL|baseline_tested=NULL|baseline_name="DRIMSeq"|alpha=0.05',
    repliapa_compare_contexts = 'primary=<required>|external=<required>|primary_label="primary"|external_label="external"|targets=NULL|fdr=0.05',
    repliapa_create_freeze = 'config=<required>|pas_ids=<required>|primary_metrics=<required>|role_registry=<required>|source_files=NULL|frozen_on=Sys.Date()|status="FROZEN"',
    repliapa_drimseq_input = 'counts=<required>|sample_data=<required>|pas_data=<required>|target=<required>|gene="gene"|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|replicates=NULL|pas_ids=NULL',
    repliapa_example = "",
    repliapa_fit = 'counts=<required>|sample_data=<required>|pas_data=<required>|target=<required>|gene="gene"|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|replicates=NULL|pseudocount=0.5|min_gene_total=20|max_combinations=1e+06',
    repliapa_fit_proximal_distal = 'counts=<required>|sample_data=<required>|pas_data=<required>|target=<required>|gene="gene"|position="position"|strand="strand"|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|replicates=NULL|pseudocount=0.5|min_gene_total=20|max_combinations=1e+06',
    repliapa_fit_targets = 'counts=<required>|sample_data=<required>|pas_data=<required>|targets=NULL|gene="gene"|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|replicates=NULL|pseudocount=0.5|min_gene_total=20|max_combinations=1e+06|strict=FALSE',
    repliapa_gene_expression_effect = 'rna_counts=<required>|sample_data=<required>|target=<required>|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|cell_n="cell_n"|replicates=NULL|pseudocount=0.5',
    repliapa_guide_reproducibility = 'guide_effects=<required>|guide_data=<required>|match_columns=<required>|replicate="replicate"|target="target"|guide="guide"|n_boot=2000L|seed=20260812L',
    repliapa_ntc_calibration = 'counts=<required>|sample_data=<required>|pas_data=<required>|contrasts=NULL|guide_sizes=c(2L, 3L, 4L)|n_contrasts=500L|seed=1L|gene="gene"|replicate="replicate"|guide="guide"|target_col="target"|ntc="ntc"|replicates=NULL|pseudocount=0.5|min_gene_total=20|max_combinations=1e+06',
    repliapa_orthogonal_validation = 'effect_table=<required>|orthogonal_data=<required>|freeze=<required>|orthogonal_registry_id=<required>|control="NT"|targets=NULL|construct_columns=c("gene_id", "pas_id", "aim", "subaim")|replicate="replicate"|perturbation="perturbation"|proximal="proximal"|distal="distal"|total="total"|experiment=NULL|experiment_col="experiment"|min_replicates=2L',
    repliapa_palette = 'type=c("method", "effect", "context", "pas")|n=8L',
    repliapa_plot_allpas = "comparison=<required>",
    repliapa_plot_apa_ge = "classification=<required>|ge_threshold=0.25|label_n=8L",
    repliapa_plot_calibration = "calibration=<required>|max_alpha=0.1",
    repliapa_plot_context = "comparison=<required>",
    repliapa_plot_guide_evidence = "fit=<required>|gene=<required>|limits=NULL",
    repliapa_plot_pas_effect = "fit=<required>|gene=<required>|pas_data=NULL|show_replicates=TRUE",
    repliapa_plot_reproducibility = "reproducibility=<required>|label_n=6L",
    repliapa_proximal_directions = "fits=<required>|pd_data=<required>",
    repliapa_pseudobulk = 'counts=<required>|cell_data=<required>|cell_line="cell_line"|replicate="replicate"|guide="guide"|target="target"|ntc="ntc"',
    repliapa_read_pas_counts = 'x=<required>|cell_data=NULL|pas_data=NULL|assay=NULL|layer="counts"',
    repliapa_reproducibility = "development=<required>|validation=<required>|fdr=0.05",
    repliapa_result_bundle = "x=<required>|freeze=NULL|input_audit=NULL|compute_benchmark=NULL",
    repliapa_theme = 'base_size=11|base_family=""|grid=c("y", "x", "xy", "none")',
    repliapa_validate_counts = "counts=<required>",
    repliapa_validate_pseudobulk = 'counts=<required>|sample_data=<required>|replicate="replicate"|guide="guide"|target="target"|ntc="ntc"',
    repliapa_verify_bundle = "path=<required>",
    repliapa_verify_contract = "contract=repliapa_api_contract()|bundle=NULL",
    repliapa_verify_freeze = "freeze=<required>|check_source_files=FALSE",
    repliapa_write_bundle = "bundle=<required>|path=<required>|overwrite=FALSE"
  )

  classes <- list(
    repliapa_api_contract = c("contract_version", "candidate_api_version", "functions", "s3_methods", "classes", "bundle_tables", "enums", "compatibility_policy", "fingerprint"),
    repliapa_baseline_comparison = c("gene_table", "native_summary", "matched_summary", "agreement", "baseline_name", "alpha"),
    repliapa_context_comparison = c("results", "target_summary", "target_accounting", "contexts", "fdr"),
    repliapa_drimseq_input = c("counts", "samples", "design", "accounting", "target", "replicates", "target_guides", "ntc_guides"),
    repliapa_fit = c("target", "replicates", "target_guides", "ntc_guides", "gene_table", "pas_effects", "guide_effects", "heterogeneity", "settings"),
    repliapa_fit_list = c("fits", "accounting", "targets"),
    repliapa_freeze = c("schema_version", "status", "frozen_on", "config", "pas_ids", "primary_metrics", "role_registry", "source_hashes", "r_version", "fingerprint"),
    repliapa_guide_reproducibility = c("results", "pair_summary", "target_summary", "matching", "bootstrap", "settings"),
    repliapa_input_audit = c("summary", "contract_checks", "sample_units", "guide_support", "target_replicate_support", "target_support", "pas_adequacy", "gene_adequacy", "settings"),
    repliapa_input = c("counts", "cell_data", "pas_data", "source", "assay", "layer"),
    repliapa_null_calibration = c("contrasts", "results", "summary", "ntc_guides", "replicates"),
    repliapa_orthogonal_validation = c("results", "target_summary", "accounting", "paired_replicates", "settings"),
    repliapa_pd_data = c("counts", "pas_data", "membership", "accounting", "settings"),
    repliapa_pd_fit = c("target", "replicates", "target_guides", "ntc_guides", "gene_table", "pas_effects", "guide_effects", "heterogeneity", "settings", "proximal_distal"),
    repliapa_pseudobulk = c("counts", "sample_data"),
    repliapa_result_bundle = c("tables", "schema_version"),
    repliapa_contract_verification = c("passed", "checks", "contract_version", "fingerprint")
  )

  tables <- list(
    target_effects = .contract_table(c("target", "gene", "pas_n", "total_count", "target_guide_n", "ntc_guide_n", "permutation_n", "statistic", "p_value", "effect_norm", "coordinate_shift", "tested", "filtered", "abstained", "reason", "q_value", "significant")),
    pas_effects = .contract_table(c("target", "gene", "pas_id", "effect"), extension_patterns = "^effect_.+$", extension_type = "double"),
    guide_effects = .contract_table(c("target", "gene", "pas_id", "replicate", "guide", "effect")),
    target_heterogeneity = .contract_table(c("target", "gene", "replicate", "guide_n", "mean_pairwise_cosine", "effect_dispersion")),
    target_accounting = .contract_table(c("target", "status", "reason", "target_guide_n", "ntc_guide_n", "tested_gene_n", "filtered_gene_n", "abstained_gene_n", "discovery_n")),
    tested_gene_accounting = .contract_table(c("target", "status", "native_gene_n", "tested_gene_n", "filtered_gene_n", "abstained_gene_n", "discovery_n", "tested_fraction", "discovery_fraction_of_tested")),
    provenance = .contract_table(c("bundle_schema_version", "package_version", "target_n_requested", "target_n_fitted", "freeze_fingerprint", "input_audit_included", "compute_benchmark_included", "inference_changed", "runtime_status", "peak_memory_status")),
    input_audit_summary = .contract_table(c("data_contract_status", "sample_unit_n", "qualified_sample_unit_n", "replicate_n", "guide_n", "target_n_non_ntc", "eligible_target_n", "pas_n", "primary_pas_n", "gene_n", "primary_multiPAS_gene_n"), required = FALSE),
    input_contract_checks = .contract_table(c("check", "passed", "detail"), required = FALSE),
    sample_unit_manifest = .contract_table(c("sample_unit", "replicate", "guide", "target", "ntc", "cell_n", "polya_count", "qualified", "qualification_reason"), required = FALSE),
    guide_support = .contract_table(c("guide", "target", "ntc", "replicate_n", "qualified_replicate_n", "sample_unit_n", "qualified_sample_unit_n", "cell_n", "polya_count"), required = FALSE),
    target_guide_support = .contract_table(c("target", "ntc", "guide_n_total", "qualified_guide_n", "replicate_n", "min_qualified_guides_per_replicate", "target_eligible", "exclusion_reason"), required = FALSE),
    pas_adequacy = .contract_table(c("pas_id", "gene", "development_total_count", "development_nonzero_sample_units", "all_sample_total_count", "mapping_unique", "exact_coordinate_conflict", "adequacy_tier", "primary_eligible"), required = FALSE),
    gene_pas_adequacy = .contract_table(c("gene", "pas_n_total", "pas_n_A0", "pas_n_A1", "pas_n_A2", "pas_n_A3", "usable_pas_n", "development_total_count", "coordinate_conflict", "primary_multiPAS_eligible", "exclusion_reason"), required = FALSE),
    compute_benchmark = .contract_table(
      c("method", "repetitions", "elapsed_seconds_median", "peak_rss_mb_median", "input_sha256", "elapsed_scope", "peak_rss_scope", "peak_rss_source"),
      required = FALSE,
      optional_columns = c("elapsed_seconds_min", "elapsed_seconds_max", "peak_rss_mb_min", "peak_rss_mb_max", "peak_above_baseline_rss_mb_median", "tested_gene_n", "reference_status", "package_version", "r_version")
    )
  )

  tables$target_effects$column_types <- .contract_types(
    character = c("target", "gene", "reason"),
    integer = c("pas_n", "target_guide_n", "ntc_guide_n"),
    double = c("total_count", "statistic", "p_value", "effect_norm", "coordinate_shift", "q_value"),
    numeric = "permutation_n",
    logical = c("tested", "filtered", "abstained", "significant")
  )
  tables$pas_effects$column_types <- .contract_types(
    character = c("target", "gene", "pas_id"), double = "effect"
  )
  tables$guide_effects$column_types <- .contract_types(
    character = c("target", "gene", "pas_id", "replicate", "guide"),
    double = "effect"
  )
  tables$target_heterogeneity$column_types <- .contract_types(
    character = c("target", "gene", "replicate"), integer = "guide_n",
    double = c("mean_pairwise_cosine", "effect_dispersion")
  )
  tables$target_accounting$column_types <- .contract_types(
    character = c("target", "status", "reason"),
    integer = c("target_guide_n", "ntc_guide_n", "tested_gene_n", "filtered_gene_n", "abstained_gene_n", "discovery_n")
  )
  tables$tested_gene_accounting$column_types <- .contract_types(
    character = c("target", "status"),
    integer = c("native_gene_n", "tested_gene_n", "filtered_gene_n", "abstained_gene_n", "discovery_n"),
    double = c("tested_fraction", "discovery_fraction_of_tested")
  )
  tables$provenance$column_types <- .contract_types(
    character = c("package_version", "freeze_fingerprint", "runtime_status", "peak_memory_status"),
    integer = c("bundle_schema_version", "target_n_requested", "target_n_fitted"),
    logical = c("input_audit_included", "compute_benchmark_included", "inference_changed")
  )
  tables$input_audit_summary$column_types <- .contract_types(
    character = "data_contract_status",
    integer = c("sample_unit_n", "qualified_sample_unit_n", "replicate_n", "guide_n", "target_n_non_ntc", "eligible_target_n", "pas_n", "primary_pas_n", "gene_n", "primary_multiPAS_gene_n")
  )
  tables$input_contract_checks$column_types <- .contract_types(
    character = c("check", "detail"), logical = "passed"
  )
  tables$sample_unit_manifest$column_types <- .contract_types(
    character = c("sample_unit", "replicate", "guide", "target", "qualification_reason"),
    integer = c("cell_n", "polya_count"), logical = c("ntc", "qualified")
  )
  tables$guide_support$column_types <- .contract_types(
    character = c("guide", "target"),
    integer = c("replicate_n", "qualified_replicate_n", "sample_unit_n", "qualified_sample_unit_n", "cell_n", "polya_count"),
    logical = "ntc"
  )
  tables$target_guide_support$column_types <- .contract_types(
    character = c("target", "exclusion_reason"),
    integer = c("guide_n_total", "qualified_guide_n", "replicate_n", "min_qualified_guides_per_replicate"),
    logical = c("ntc", "target_eligible")
  )
  tables$pas_adequacy$column_types <- .contract_types(
    character = c("pas_id", "gene", "adequacy_tier"),
    integer = c("development_total_count", "development_nonzero_sample_units", "all_sample_total_count"),
    logical = c("mapping_unique", "exact_coordinate_conflict", "primary_eligible")
  )
  tables$gene_pas_adequacy$column_types <- .contract_types(
    character = c("gene", "exclusion_reason"),
    integer = c("pas_n_total", "pas_n_A0", "pas_n_A1", "pas_n_A2", "pas_n_A3", "usable_pas_n", "development_total_count"),
    logical = c("coordinate_conflict", "primary_multiPAS_eligible")
  )
  tables$compute_benchmark$column_types <- .contract_types(
    character = c("method", "reference_status", "input_sha256", "elapsed_scope", "peak_rss_scope", "peak_rss_source", "package_version", "r_version"),
    numeric = c("repetitions", "tested_gene_n"),
    double = c("elapsed_seconds_median", "elapsed_seconds_min", "elapsed_seconds_max", "peak_rss_mb_median", "peak_rss_mb_min", "peak_rss_mb_max", "peak_above_baseline_rss_mb_median")
  )

  list(
    contract_version = "1.0.0-rc4",
    candidate_api_version = "1.0.0",
    package_version_floor = "0.99.0",
    bundle_schema_version = 1L,
    functions = functions,
    s3_methods = paste0("print.", c(
      "repliapa_api_contract", "repliapa_baseline_comparison",
      "repliapa_context_comparison", "repliapa_contract_verification",
      "repliapa_drimseq_input", "repliapa_fit", "repliapa_fit_list",
      "repliapa_freeze", "repliapa_guide_reproducibility",
      "repliapa_input", "repliapa_input_audit", "repliapa_null_calibration",
      "repliapa_orthogonal_validation", "repliapa_result_bundle"
    )),
    classes = classes,
    bundle_tables = tables,
    enums = list(
      fit_status = c("FITTED", "UNSUPPORTED"),
      fit_reason = c(
        "NONE", "TOO_MANY_ADMISSIBLE_PERMUTATIONS",
        "TOO_MANY_EXACT_PERMUTATIONS", # accepted legacy machine code; deprecated
        "SINGLE_PAS_GENE", "LOW_GENE_TOTAL"
      ),
      input_contract_status = c("PASS", "PASS_WITH_EXCLUSIONS"),
      pas_adequacy_tier = paste0("PAS-A", 0:3),
      effect_class = c("APA_AND_GE", "APA_ONLY", "GE_ONLY", "NEITHER"),
      context_category = c("UNSUPPORTED_GENE_CONTEXT", "SHARED_EFFECT", "OPPOSITE_EFFECT", "CELL_LINE_SPECIFIC_EFFECT", "SUPPORTED_NEITHER_SIGNIFICANT"),
      context_comparison_status = c("SUPPORTED_BOTH", "EXTERNAL_UNSUPPORTED", "PRIMARY_UNSUPPORTED", "BOTH_UNSUPPORTED"),
      runtime_status = c("NOT_INSTRUMENTED", "MEASURED_ISOLATED_PROCESS_MEDIAN"),
      peak_memory_status = c("NOT_INSTRUMENTED", "MEASURED_KERNEL_VM_HWM")
    ),
    compatibility_policy = list(
      public_function_signature = "exact until an explicitly versioned breaking release",
      s3_components = "required components may not be removed or renamed",
      bundle_columns = "required columns may not be removed or renamed; declared extension patterns are additive",
      enum_values = "existing meanings are stable; new values require release-note disclosure",
      inference = "API/schema maintenance must not alter scientific estimands, frozen thresholds, or historical results"
    )
  )
}

.repliapa_contract_fingerprint <- "29467b553b7df6ed26a21a8dac19eca3f223f6afac5f249d65eeeb3ed68addff"

#' Inspect the frozen v1-candidate API and result-schema contract
#'
#' Returns the machine-readable contract for exported function signatures, S3
#' result components, bundle tables, enumerations, and compatibility policy.
#' The fingerprint fails closed if the package-internal contract is edited
#' without deliberately updating the candidate contract.
#'
#' @return A `repliapa_api_contract` object.
#' @export
repliapa_api_contract <- function() {
  payload <- .repliapa_contract_payload()
  fingerprint <- digest::digest(payload, algo = "sha256", serialize = TRUE,
                                serializeVersion = 2L)
  if (!identical(fingerprint, .repliapa_contract_fingerprint)) {
    stop("The internal API contract does not match its frozen fingerprint.",
         call. = FALSE)
  }
  payload$fingerprint <- fingerprint
  class(payload) <- "repliapa_api_contract"
  unserialize(serialize(payload, NULL, version = 2L))
}

#' @export
print.repliapa_api_contract <- function(x, ...) {
  cat("RepliAPA v1-candidate API contract\n")
  cat("  contract / candidate: ", x$contract_version, " / ",
      x$candidate_api_version, "\n", sep = "")
  cat("  exported functions / S3 classes / bundle tables: ",
      length(x$functions), " / ", length(x$classes), " / ",
      length(x$bundle_tables), "\n", sep = "")
  cat("  SHA256: ", x$fingerprint, "\n", sep = "")
  invisible(x)
}

.contract_check_row <- function(check, passed, detail) {
  data.frame(check = check, passed = isTRUE(passed), detail = detail,
             stringsAsFactors = FALSE)
}

#' Verify the v1-candidate API and an optional result bundle
#'
#' Checks the supplied contract fingerprint against the package's frozen
#' contract, verifies exported function signatures and registered print
#' methods, and optionally validates a result bundle's required tables and
#' columns. The function reports every check and does not rerun inference.
#'
#' @param contract A contract returned by `repliapa_api_contract()`.
#' @param bundle Optional `repliapa_result_bundle` to validate.
#' @return A `repliapa_contract_verification` object.
#' @export
repliapa_verify_contract <- function(contract = repliapa_api_contract(),
                                     bundle = NULL) {
  if (!inherits(contract, "repliapa_api_contract")) {
    stop("`contract` must be a repliapa_api_contract object.", call. = FALSE)
  }
  supplied_fingerprint <- contract$fingerprint
  payload <- unclass(contract)
  payload$fingerprint <- NULL
  calculated <- digest::digest(payload, algo = "sha256", serialize = TRUE,
                               serializeVersion = 2L)
  checks <- list(.contract_check_row(
    "CONTRACT_FINGERPRINT",
    identical(supplied_fingerprint, .repliapa_contract_fingerprint) &&
      identical(calculated, supplied_fingerprint),
    calculated
  ))

  namespace <- asNamespace("RepliAPA")
  actual_exports <- sort(getNamespaceExports(namespace))
  expected_exports <- sort(names(contract$functions))
  checks[[length(checks) + 1L]] <- .contract_check_row(
    "EXPORTED_FUNCTION_SET", identical(actual_exports, expected_exports),
    paste(actual_exports, collapse = ",")
  )
  for (function_name in intersect(expected_exports, actual_exports)) {
    actual <- .canonical_formals(get(function_name, envir = namespace,
                                     inherits = FALSE))
    expected <- unname(contract$functions[[function_name]])
    checks[[length(checks) + 1L]] <- .contract_check_row(
      paste0("SIGNATURE:", function_name), identical(actual, expected), actual
    )
  }
  for (method_name in contract$s3_methods) {
    present <- exists(method_name, envir = namespace, inherits = FALSE) &&
      is.function(get(method_name, envir = namespace, inherits = FALSE))
    checks[[length(checks) + 1L]] <- .contract_check_row(
      paste0("S3_METHOD:", method_name), present,
      if (present) "PRESENT" else "MISSING"
    )
  }

  if (!is.null(bundle)) {
    if (!inherits(bundle, "repliapa_result_bundle")) {
      stop("`bundle` must be a repliapa_result_bundle object.", call. = FALSE)
    }
    expected_tables <- names(contract$bundle_tables)
    actual_tables <- names(bundle$tables)
    required_tables <- expected_tables[vapply(contract$bundle_tables, `[[`,
                                               logical(1), "required")]
    checks[[length(checks) + 1L]] <- .contract_check_row(
      "BUNDLE_REQUIRED_TABLES",
      all(required_tables %in% actual_tables),
      paste(actual_tables, collapse = ",")
    )
    checks[[length(checks) + 1L]] <- .contract_check_row(
      "BUNDLE_DECLARED_TABLES", all(actual_tables %in% expected_tables),
      paste(setdiff(actual_tables, expected_tables), collapse = ",")
    )
    for (table_name in intersect(expected_tables, actual_tables)) {
      specification <- contract$bundle_tables[[table_name]]
      actual_columns <- names(bundle$tables[[table_name]])
      missing <- setdiff(specification$required_columns, actual_columns)
      declared <- c(specification$required_columns, specification$optional_columns)
      extra <- setdiff(actual_columns, declared)
      allowed_extra <- if (!length(extra)) logical() else vapply(extra, function(column) {
        any(vapply(specification$extension_patterns, grepl, logical(1), x = column))
      }, logical(1))
      undeclared <- extra[!allowed_extra]
      checks[[length(checks) + 1L]] <- .contract_check_row(
        paste0("TABLE_SCHEMA:", table_name),
        !length(missing) && !length(undeclared),
        paste0("missing=", paste(missing, collapse = ","),
               ";undeclared=", paste(undeclared, collapse = ","))
      )
      expected_types <- specification$column_types[
        intersect(names(specification$column_types), actual_columns)
      ]
      actual_types <- vapply(
        bundle$tables[[table_name]][names(expected_types)], typeof, character(1)
      )
      extension_columns <- setdiff(extra[allowed_extra], names(expected_types))
      extension_passed <- !length(extension_columns) ||
        is.null(specification$extension_type) ||
        all(vapply(bundle$tables[[table_name]][extension_columns], typeof,
                   character(1)) == specification$extension_type)
      type_match <- actual_types == expected_types |
        (expected_types == "numeric" & actual_types %in% c("integer", "double"))
      type_passed <- all(type_match) && extension_passed
      checks[[length(checks) + 1L]] <- .contract_check_row(
        paste0("TABLE_TYPES:", table_name), type_passed,
        paste(paste(names(actual_types), actual_types, sep = "="), collapse = ",")
      )
    }
  }

  checks <- do.call(rbind, checks)
  rownames(checks) <- NULL
  output <- list(
    passed = all(checks$passed), checks = checks,
    contract_version = contract$contract_version,
    fingerprint = supplied_fingerprint
  )
  class(output) <- "repliapa_contract_verification"
  output
}

#' @export
print.repliapa_contract_verification <- function(x, ...) {
  cat("RepliAPA API/schema contract verification\n")
  cat("  status: ", if (x$passed) "PASS" else "FAIL", "\n", sep = "")
  cat("  checks / failures: ", nrow(x$checks), " / ",
      sum(!x$checks$passed), "\n", sep = "")
  invisible(x)
}
