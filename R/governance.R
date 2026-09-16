.freeze_payload <- function(x) {
  x[c("schema_version", "status", "frozen_on", "config", "pas_ids",
      "primary_metrics", "role_registry", "source_hashes", "r_version")]
}

.freeze_fingerprint <- function(payload) {
  digest::digest(payload, algo = "sha256", serialize = TRUE, serializeVersion = 2L)
}

.validate_role_registry <- function(role_registry) {
  if (!is.data.frame(role_registry)) {
    stop("`role_registry` must be a data frame.", call. = FALSE)
  }
  .require_columns(role_registry, c("dataset", "role"), "role_registry")
  registry <- as.data.frame(role_registry, stringsAsFactors = FALSE)
  if (!"registry_id" %in% names(registry)) {
    registry$registry_id <- sprintf("ROLE_%03d", seq_len(nrow(registry)))
  }
  if (!nrow(registry) || anyNA(registry$registry_id) || any(registry$registry_id == "") ||
      anyDuplicated(registry$registry_id)) {
    stop("Role registry IDs must be non-missing and unique.", call. = FALSE)
  }
  if (anyNA(registry$dataset) || any(registry$dataset == "") ||
      anyNA(registry$role) || any(registry$role == "")) {
    stop("Dataset and role assignments must be complete.", call. = FALSE)
  }
  registry
}

#' Create an immutable analysis freeze contract
#'
#' The returned object records the statistical configuration, fixed PAS set,
#' primary metrics, data roles, and optional source-file SHA256 checksums. A
#' second SHA256 fingerprint protects the complete payload from silent changes.
#'
#' @param config Named list containing the frozen model and filter settings.
#' @param pas_ids Character vector of frozen PAS identifiers.
#' @param primary_metrics Named list or character vector of primary metrics.
#' @param role_registry Data frame with at least `dataset` and `role` columns.
#' @param source_files Optional named character vector of governance files to
#'   checksum. Names are used as stable source labels.
#' @param frozen_on Explicit freeze date or timestamp.
#' @param status Freeze status label.
#' @return A `repliapa_freeze` object.
#' @export
repliapa_create_freeze <- function(
    config, pas_ids, primary_metrics, role_registry, source_files = NULL,
    frozen_on = Sys.Date(), status = "FROZEN") {
  if (!is.list(config) || is.null(names(config)) || any(names(config) == "")) {
    stop("`config` must be a named list.", call. = FALSE)
  }
  pas_ids <- as.character(pas_ids)
  if (!length(pas_ids) || anyNA(pas_ids) || any(pas_ids == "") || anyDuplicated(pas_ids)) {
    stop("`pas_ids` must be non-missing and unique.", call. = FALSE)
  }
  if (!(is.list(primary_metrics) || is.character(primary_metrics)) ||
      !length(primary_metrics)) {
    stop("`primary_metrics` must be a non-empty list or character vector.", call. = FALSE)
  }
  registry <- .validate_role_registry(role_registry)
  if (length(frozen_on) != 1L || is.na(frozen_on) || as.character(frozen_on) == "") {
    stop("`frozen_on` must be one explicit date or timestamp.", call. = FALSE)
  }
  if (length(status) != 1L || is.na(status) || status == "") {
    stop("`status` must be one non-missing label.", call. = FALSE)
  }
  source_hashes <- data.frame(
    source = character(), path = character(), sha256 = character(),
    stringsAsFactors = FALSE
  )
  if (!is.null(source_files)) {
    source_names <- names(source_files)
    source_files <- stats::setNames(as.character(source_files), source_names)
    if (!length(source_files) || anyNA(source_files) || any(source_files == "") ||
        is.null(names(source_files)) || any(names(source_files) == "") ||
        anyDuplicated(names(source_files)) || any(!file.exists(source_files))) {
      stop("`source_files` must be a named vector of existing files.", call. = FALSE)
    }
    source_hashes <- data.frame(
      source = names(source_files), path = normalizePath(source_files, mustWork = TRUE),
      sha256 = vapply(source_files, digest::digest, character(1),
                      file = TRUE, algo = "sha256", serialize = FALSE),
      stringsAsFactors = FALSE
    )
  }
  output <- list(
    schema_version = 1L, status = status, frozen_on = as.character(frozen_on),
    config = config, pas_ids = pas_ids, primary_metrics = primary_metrics,
    role_registry = registry, source_hashes = source_hashes,
    r_version = paste(R.version$major, R.version$minor, sep = ".")
  )
  output$fingerprint <- .freeze_fingerprint(.freeze_payload(output))
  class(output) <- "repliapa_freeze"
  output
}

#' Verify an analysis freeze and optional source files
#'
#' @param freeze A `repliapa_freeze` object.
#' @param check_source_files If `TRUE`, require every recorded source path to
#'   exist and still match its recorded SHA256.
#' @return `TRUE`; otherwise stops on any integrity failure.
#' @export
repliapa_verify_freeze <- function(freeze, check_source_files = FALSE) {
  if (!inherits(freeze, "repliapa_freeze")) {
    stop("`freeze` must be a repliapa_freeze object.", call. = FALSE)
  }
  required <- c("schema_version", "status", "frozen_on", "config", "pas_ids",
                "primary_metrics", "role_registry", "source_hashes", "r_version",
                "fingerprint")
  if (!all(required %in% names(freeze))) stop("Freeze payload is incomplete.", call. = FALSE)
  observed <- .freeze_fingerprint(.freeze_payload(freeze))
  if (!identical(observed, freeze$fingerprint)) {
    stop("Freeze fingerprint mismatch: payload was modified.", call. = FALSE)
  }
  .validate_role_registry(freeze$role_registry)
  if (isTRUE(check_source_files) && nrow(freeze$source_hashes)) {
    missing <- !file.exists(freeze$source_hashes$path)
    if (any(missing)) stop("A frozen governance source file is missing.", call. = FALSE)
    observed_hashes <- vapply(
      freeze$source_hashes$path, digest::digest, character(1),
      file = TRUE, algo = "sha256", serialize = FALSE
    )
    if (!identical(unname(observed_hashes), unname(freeze$source_hashes$sha256))) {
      stop("A frozen governance source file changed.", call. = FALSE)
    }
  }
  TRUE
}

.role_action_policy <- list(
  DEVELOPMENT = c("schema_audit", "debug", "tune_model", "tune_filter",
                  "change_pas_set", "change_metric", "fit_model", "read_outcome"),
  FROZEN_INTERNAL_VALIDATION = c("schema_audit", "apply_frozen_model", "read_outcome"),
  EXTERNAL_CELL_CONTEXT_VALIDATION = c("schema_audit", "apply_frozen_model", "read_outcome"),
  NEGATIVE_CONTROL = c("schema_audit", "null_calibration", "read_outcome"),
  NULL_CALIBRATION = c("schema_audit", "null_calibration", "read_outcome"),
  OPTIONAL_ORTHOGONAL_VALIDATION = c("schema_audit", "orthogonal_validation", "read_outcome"),
  METADATA_ONLY_APPLICATION_CANDIDATE = c("metadata_audit")
)

.governance_actions <- c(
  "schema_audit", "metadata_audit", "debug", "tune_model", "tune_filter",
  "change_pas_set", "change_metric", "fit_model", "apply_frozen_model",
  "read_outcome", "null_calibration", "split_guide", "orthogonal_validation",
  "download_raw", "claim_independent_validation"
)

#' Authorize an action against a frozen data-role registry
#'
#' Composite roles separated by semicolons are allowed when at least one role
#' permits the action. Unknown roles and actions fail closed.
#'
#' @param freeze A verified `repliapa_freeze` object.
#' @param registry_id Stable registry ID stored in `freeze$role_registry`.
#' @param action Requested action label.
#' @return `TRUE` invisibly; otherwise stops with the role violation.
#' @export
repliapa_authorize <- function(freeze, registry_id, action) {
  repliapa_verify_freeze(freeze)
  if (length(registry_id) != 1L || is.na(registry_id) || registry_id == "") {
    stop("`registry_id` must be one non-missing label.", call. = FALSE)
  }
  if (length(action) != 1L || is.na(action) || action == "") {
    stop("`action` must be one non-missing label.", call. = FALSE)
  }
  if (!action %in% .governance_actions) stop("Unknown governance action.", call. = FALSE)
  row <- freeze$role_registry[freeze$role_registry$registry_id == registry_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Registry ID is unavailable or ambiguous.", call. = FALSE)
  roles <- trimws(strsplit(as.character(row$role[[1]]), ";", fixed = TRUE)[[1]])
  if (any(!roles %in% names(.role_action_policy))) {
    stop("Unknown data role in frozen registry.", call. = FALSE)
  }
  allowed <- any(vapply(roles, function(role) action %in% .role_action_policy[[role]], logical(1)))
  if (!allowed) {
    stop("Action `", action, "` is prohibited for role `",
         paste(roles, collapse = ";"), "`.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Apply the frozen method to an external context and compare effects
#'
#' @inheritParams repliapa_compare_contexts
#' @param freeze Verified analysis freeze.
#' @param external_registry_id Registry entry carrying the external-validation role.
#' @return A `repliapa_context_comparison` with freeze provenance attached.
#' @export
repliapa_apply_external <- function(primary, external, freeze, external_registry_id,
                                    primary_label = "primary", external_label = "external",
                                    targets = NULL, fdr = 0.05) {
  repliapa_authorize(freeze, external_registry_id, "apply_frozen_model")
  output <- repliapa_compare_contexts(
    primary, external, primary_label = primary_label,
    external_label = external_label, targets = targets, fdr = fdr
  )
  output$freeze_fingerprint <- freeze$fingerprint
  output$external_registry_id <- external_registry_id
  output
}

#' @export
print.repliapa_freeze <- function(x, ...) {
  cat("RepliAPA analysis freeze\n")
  cat("  status / frozen on: ", x$status, " / ", x$frozen_on, "\n", sep = "")
  cat("  PAS / role records: ", length(x$pas_ids), " / ", nrow(x$role_registry), "\n", sep = "")
  cat("  SHA256 fingerprint: ", x$fingerprint, "\n", sep = "")
  invisible(x)
}
