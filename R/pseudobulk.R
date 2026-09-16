#' Aggregate cell-level PAS counts to biological sample units
#'
#' Cells are aggregated into cell-line by biological-replicate by guide units.
#' The function never creates cell-level inferential replicates.
#'
#' @param counts PAS-by-cell integer count matrix.
#' @param cell_data Cell metadata aligned to columns of `counts`.
#' @param cell_line,replicate,guide,target,ntc Column names in `cell_data`.
#' @return A list containing `counts` and standardized `sample_data`.
#' @export
repliapa_pseudobulk <- function(counts, cell_data,
                                cell_line = "cell_line", replicate = "replicate",
                                guide = "guide", target = "target", ntc = "ntc") {
  repliapa_validate_counts(counts)
  .aligned_data(cell_data, colnames(counts), "cell_data")
  .require_columns(cell_data, c(cell_line, replicate, guide, target, ntc), "cell_data")
  if (!is.logical(cell_data[[ntc]]) || anyNA(cell_data[[ntc]])) {
    stop("The NTC indicator must be a non-missing logical column.", call. = FALSE)
  }
  assignments <- cell_data[c(cell_line, replicate, guide, target, ntc)]
  for (column in c(cell_line, replicate, guide, target)) {
    if (anyNA(assignments[[column]]) || any(as.character(assignments[[column]]) == "")) {
      stop("`", column, "` has missing assignments.", call. = FALSE)
    }
  }
  unit <- paste(assignments[[cell_line]], assignments[[replicate]], assignments[[guide]], sep = "|")
  unit_levels <- sort(unique(unit))
  unit_index <- match(unit, unit_levels)
  design <- Matrix::sparseMatrix(i = seq_along(unit_index), j = unit_index, x = 1,
                                 dims = c(length(unit_index), length(unit_levels)),
                                 dimnames = list(colnames(counts), unit_levels))
  aggregated <- counts %*% design
  if (!isTRUE(all.equal(as.numeric(Matrix::rowSums(aggregated)), as.numeric(Matrix::rowSums(counts)),
                        tolerance = 0))) {
    stop("Pseudobulk aggregation did not conserve PAS counts.", call. = FALSE)
  }
  first <- match(unit_levels, unit)
  sample_data <- data.frame(
    cell_line = as.character(assignments[[cell_line]][first]),
    replicate = as.character(assignments[[replicate]][first]),
    guide = as.character(assignments[[guide]][first]),
    target = as.character(assignments[[target]][first]),
    ntc = assignments[[ntc]][first],
    cell_n = tabulate(unit_index, nbins = length(unit_levels)),
    row.names = unit_levels,
    stringsAsFactors = FALSE
  )
  # A guide must map to one target and one NTC state across all of its cells.
  consistent <- vapply(unit_levels, function(z) {
    rows <- unit == z
    length(unique(assignments[[target]][rows])) == 1L &&
      length(unique(assignments[[ntc]][rows])) == 1L
  }, logical(1))
  if (!all(consistent)) stop("A sample unit maps to multiple target or NTC states.", call. = FALSE)
  repliapa_validate_pseudobulk(aggregated, sample_data)
  structure(list(counts = aggregated, sample_data = sample_data),
            class = "repliapa_pseudobulk")
}
