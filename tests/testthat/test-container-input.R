test_that("matrix input is standardized and metadata is aligned", {
  x <- repliapa_example()
  cell_data <- x$sample_data[rev(rownames(x$sample_data)), , drop = FALSE]
  pas_data <- x$pas_data[rev(rownames(x$pas_data)), , drop = FALSE]

  input <- repliapa_read_pas_counts(x$counts, cell_data, pas_data)

  expect_s3_class(input, "repliapa_input")
  expect_identical(input$source, "MATRIX")
  expect_identical(rownames(input$cell_data), colnames(input$counts))
  expect_identical(rownames(input$pas_data), rownames(input$counts))
  expect_identical(input$counts, x$counts)
  expect_output(print(input), "RepliAPA standardized PAS input")
})

test_that("matrix input fails closed on incomplete or mismatched metadata", {
  x <- repliapa_example()
  expect_error(repliapa_read_pas_counts(x$counts), "requires both")

  bad <- x$sample_data[-1, , drop = FALSE]
  expect_error(
    repliapa_read_pas_counts(x$counts, bad, x$pas_data),
    "identifiers must match"
  )
})

test_that("Seurat input uses an explicit PAS assay", {
  skip_if_not_installed("SeuratObject")
  x <- repliapa_example()
  object <- SeuratObject::CreateSeuratObject(
    counts = methods::as(x$counts, "dgCMatrix"),
    assay = "polyA",
    meta.data = x$sample_data
  )
  object[["polyA"]][["gene"]] <- x$pas_data$gene

  expect_error(repliapa_read_pas_counts(object), "must be supplied explicitly")
  input <- repliapa_read_pas_counts(object, assay = "polyA")

  expect_s3_class(input, "repliapa_input")
  expect_identical(input$source, "SEURAT")
  expect_identical(input$assay, "polyA")
  expect_identical(input$layer, "counts")
  expect_identical(input$pas_data$gene, x$pas_data$gene)
  expect_identical(colnames(input$counts), rownames(input$cell_data))
})

test_that("SummarizedExperiment and SingleCellExperiment inputs are supported", {
  skip_if_not_installed("SummarizedExperiment")
  x <- repliapa_example()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(polyA = x$counts),
    colData = x$sample_data,
    rowData = x$pas_data
  )
  se_input <- repliapa_read_pas_counts(se, assay = "polyA")
  expect_identical(se_input$source, "SUMMARIZED_EXPERIMENT")
  expect_identical(se_input$counts, x$counts)
  expect_identical(se_input$pas_data$gene, x$pas_data$gene)

  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(polyA = x$counts),
    colData = x$sample_data,
    rowData = x$pas_data
  )
  sce_input <- repliapa_read_pas_counts(sce, assay = "polyA")
  expect_identical(sce_input$source, "SINGLE_CELL_EXPERIMENT")
  expect_identical(sce_input$counts, x$counts)
  expect_identical(rownames(sce_input$cell_data), colnames(x$counts))
})

test_that("container reader rejects absent assays and non-count values", {
  skip_if_not_installed("SummarizedExperiment")
  x <- repliapa_example()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(polyA = x$counts),
    colData = x$sample_data,
    rowData = x$pas_data
  )
  expect_error(repliapa_read_pas_counts(se, assay = "RNA"), "unavailable")

  bad <- x$counts
  bad[1, 1] <- 0.5
  expect_error(
    repliapa_read_pas_counts(bad, x$sample_data, x$pas_data),
    "nonnegative integers"
  )
})
