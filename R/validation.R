#' Validate a PAS count matrix
#'
#' Checks that an input is a matrix-like object containing finite,
#' nonnegative integer counts with unique feature and sample identifiers.
#'
#' @param counts A feature-by-cell or feature-by-sample matrix.
#' @return `counts`, invisibly.
#' @export
repliapa_validate_counts <- function(counts) {
  if (!(is.matrix(counts) || inherits(counts, "Matrix"))) {
    stop("`counts` must be a matrix or Matrix object.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("PAS feature identifiers must be present and unique.", call. = FALSE)
  }
  if (is.null(colnames(counts)) || anyDuplicated(colnames(counts))) {
    stop("Cell or sample identifiers must be present and unique.", call. = FALSE)
  }
  value <- if (inherits(counts, "sparseMatrix")) methods::slot(counts, "x") else as.vector(counts)
  if (any(!is.finite(value)) || any(value < 0) || any(abs(value - round(value)) > 1e-8)) {
    stop("`counts` must contain finite nonnegative integers.", call. = FALSE)
  }
  invisible(counts)
}

.require_columns <- function(x, required, object = "data") {
  missing <- setdiff(required, colnames(x))
  if (length(missing)) {
    stop(object, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_scalar_number <- function(value, name, lower = -Inf,
                                    lower_inclusive = TRUE, upper = Inf,
                                    upper_inclusive = TRUE) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value)
  if (valid) {
    valid <- if (lower_inclusive) value >= lower else value > lower
    valid <- valid && if (upper_inclusive) value <= upper else value < upper
  }
  if (!valid) stop("`", name, "` is outside its allowed finite numeric range.",
                   call. = FALSE)
  value
}

.validate_scalar_integer <- function(value, name, minimum = 0L) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value >= minimum &&
    abs(value - round(value)) <= 1e-8 && value <= .Machine$integer.max
  if (!valid) {
    stop("`", name, "` must be one integer of at least ", minimum, ".",
         call. = FALSE)
  }
  as.integer(value)
}

.validate_target_label <- function(value, name = "target") {
  if (length(value) != 1L || is.na(value) || as.character(value) == "") {
    stop("`", name, "` must be one non-missing label.", call. = FALSE)
  }
  as.character(value)
}

.validate_replicates <- function(sample_data, replicate, replicates = NULL) {
  available <- sort(unique(as.character(sample_data[[replicate]])))
  if (is.null(replicates)) replicates <- available
  replicates <- sort(unique(as.character(replicates)))
  if (!length(replicates) || anyNA(replicates) || any(replicates == "") ||
      any(!replicates %in% available)) {
    stop("Requested replicates must be non-missing available labels.", call. = FALSE)
  }
  replicates
}

.aligned_data <- function(data, ids, object) {
  if (!is.data.frame(data)) stop("`", object, "` must be a data.frame.", call. = FALSE)
  if (is.null(rownames(data)) || anyDuplicated(rownames(data))) {
    stop("`", object, "` row names must be unique identifiers.", call. = FALSE)
  }
  if (!identical(rownames(data), ids)) {
    stop("`", object, "` row names must be identical to and ordered like matrix columns/rows.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate guide-by-replicate pseudobulk input
#'
#' @param counts PAS-by-sample-unit count matrix.
#' @param sample_data Data frame aligned to columns of `counts`.
#' @param replicate,guide,target,ntc Column names in `sample_data`.
#' @return `TRUE`, invisibly.
#' @export
repliapa_validate_pseudobulk <- function(counts, sample_data,
                                         replicate = "replicate", guide = "guide",
                                         target = "target", ntc = "ntc") {
  repliapa_validate_counts(counts)
  .aligned_data(sample_data, colnames(counts), "sample_data")
  .require_columns(sample_data, c(replicate, guide, target, ntc), "sample_data")
  for (column in c(replicate, guide, target)) {
    value <- sample_data[[column]]
    if (anyNA(value) || any(as.character(value) == "")) {
      stop("`", column, "` has missing assignments.", call. = FALSE)
    }
  }
  ntc_value <- sample_data[[ntc]]
  if (!is.logical(ntc_value) || anyNA(ntc_value)) {
    stop("`", ntc, "` must be a non-missing logical vector.", call. = FALSE)
  }
  key <- paste(sample_data[[replicate]], sample_data[[guide]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("Each guide may occur at most once in each biological replicate.", call. = FALSE)
  }
  invisible(TRUE)
}

.composition <- function(counts, pseudocount = 0.5) {
  sweep(counts + pseudocount, 2L,
        colSums(counts) + pseudocount * nrow(counts), "/")
}

.cosine <- function(a, b) {
  denominator <- sqrt(sum(a^2) * sum(b^2))
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  sum(a * b) / denominator
}
