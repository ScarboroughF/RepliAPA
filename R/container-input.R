.repliapa_align_metadata <- function(data, ids, name) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (is.null(rownames(data)) || anyDuplicated(rownames(data))) {
    stop("`", name, "` row names must be unique identifiers.", call. = FALSE)
  }
  if (!setequal(rownames(data), ids)) {
    stop("`", name, "` identifiers must match the extracted count matrix.",
         call. = FALSE)
  }
  data[ids, , drop = FALSE]
}

.repliapa_require_assay <- function(assay, object) {
  if (is.null(assay) || length(assay) != 1L || is.na(assay) || assay == "") {
    stop("`assay` must be supplied explicitly for a ", object,
         " input so that an RNA assay is not mistaken for PAS counts.",
         call. = FALSE)
  }
  as.character(assay)
}

#' Read PAS counts from matrices or single-cell containers
#'
#' Standardizes a PAS-by-cell count matrix, cell metadata, and PAS metadata
#' without detecting or modifying PASs. For container inputs the PAS assay must
#' be named explicitly; this fail-closed rule prevents an RNA assay from being
#' analysed as poly(A)-site counts.
#'
#' @param x A base matrix, `Matrix`, Seurat object,
#'   `SingleCellExperiment`, or `SummarizedExperiment`.
#' @param cell_data Cell metadata for matrix input. For container input, an
#'   optional replacement for the object's column metadata.
#' @param pas_data PAS metadata for matrix input. For container input, an
#'   optional replacement for the object's feature metadata.
#' @param assay Explicit PAS assay name for container inputs. Ignored for a
#'   matrix input.
#' @param layer Seurat count layer or slot. The default is `"counts"`.
#' @return A `repliapa_input` containing aligned `counts`, `cell_data`, and
#'   `pas_data`, plus input provenance.
#' @export
repliapa_read_pas_counts <- function(x, cell_data = NULL, pas_data = NULL,
                                     assay = NULL, layer = "counts") {
  source <- NULL
  assay_used <- NA_character_
  layer_used <- NA_character_

  if (is.matrix(x) || inherits(x, "Matrix")) {
    counts <- x
    source <- "MATRIX"
    if (is.null(cell_data) || is.null(pas_data)) {
      stop("Matrix input requires both `cell_data` and `pas_data`.",
           call. = FALSE)
    }
  } else if (inherits(x, "Seurat")) {
    assay <- .repliapa_require_assay(assay, "Seurat")
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Reading a Seurat object requires the optional SeuratObject package.",
           call. = FALSE)
    }
    available <- SeuratObject::Assays(x)
    if (!assay %in% available) {
      stop("Seurat assay `", assay, "` is unavailable. Available assays: ",
           paste(available, collapse = ", "), ".", call. = FALSE)
    }
    if (length(layer) != 1L || is.na(layer) || layer == "") {
      stop("`layer` must be one non-empty Seurat count layer or slot.",
           call. = FALSE)
    }
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      counts <- SeuratObject::GetAssayData(x, assay = assay, layer = layer)
    } else {
      counts <- SeuratObject::GetAssayData(x, assay = assay, slot = layer)
    }
    if (is.null(cell_data)) cell_data <- x[[]]
    if (is.null(pas_data)) pas_data <- x[[assay]][[]]
    source <- "SEURAT"
    assay_used <- assay
    layer_used <- layer
  } else if (inherits(x, "SummarizedExperiment")) {
    assay <- .repliapa_require_assay(assay, class(x)[1L])
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("Reading this object requires the optional SummarizedExperiment package.",
           call. = FALSE)
    }
    available <- SummarizedExperiment::assayNames(x)
    if (!assay %in% available) {
      stop("Assay `", assay, "` is unavailable. Available assays: ",
           paste(available, collapse = ", "), ".", call. = FALSE)
    }
    counts <- SummarizedExperiment::assay(x, assay)
    if (is.null(cell_data)) {
      cell_data <- as.data.frame(SummarizedExperiment::colData(x),
                                 stringsAsFactors = FALSE)
    }
    if (is.null(pas_data)) {
      pas_data <- as.data.frame(SummarizedExperiment::rowData(x),
                                stringsAsFactors = FALSE)
    }
    source <- if (inherits(x, "SingleCellExperiment")) {
      "SINGLE_CELL_EXPERIMENT"
    } else {
      "SUMMARIZED_EXPERIMENT"
    }
    assay_used <- assay
  } else {
    stop("`x` must be a matrix, Matrix, Seurat, SingleCellExperiment, or SummarizedExperiment object.",
         call. = FALSE)
  }

  repliapa_validate_counts(counts)
  cell_data <- .repliapa_align_metadata(cell_data, colnames(counts), "cell_data")
  pas_data <- .repliapa_align_metadata(pas_data, rownames(counts), "pas_data")

  output <- list(
    counts = counts,
    cell_data = cell_data,
    pas_data = pas_data,
    source = source,
    assay = assay_used,
    layer = layer_used
  )
  class(output) <- "repliapa_input"
  output
}

#' @export
print.repliapa_input <- function(x, ...) {
  cat("RepliAPA standardized PAS input\n")
  cat("  source: ", x$source, "\n", sep = "")
  cat("  PAS / cells or sample units: ", nrow(x$counts), " / ",
      ncol(x$counts), "\n", sep = "")
  if (!is.na(x$assay)) cat("  assay: ", x$assay, "\n", sep = "")
  invisible(x)
}
