test_that("count validation rejects invalid inputs", {
  x <- matrix(c(0, 1, 2, 3), 2, dimnames = list(c("p1", "p2"), c("s1", "s2")))
  expect_invisible(repliapa_validate_counts(x))
  x[1, 1] <- -1
  expect_error(repliapa_validate_counts(x), "nonnegative integers")
  x[1, 1] <- 0.5
  expect_error(repliapa_validate_counts(x), "nonnegative integers")
  rownames(x) <- NULL
  expect_error(repliapa_validate_counts(x), "feature identifiers")
})

test_that("pseudobulk contract requires one guide-replicate row", {
  x <- repliapa_example()
  expect_invisible(repliapa_validate_pseudobulk(x$counts, x$sample_data))
  duplicated_data <- rbind(x$sample_data, x$sample_data[1, , drop = FALSE])
  rownames(duplicated_data)[nrow(duplicated_data)] <- "duplicate"
  duplicated_counts <- cbind(x$counts, x$counts[, 1])
  colnames(duplicated_counts)[ncol(duplicated_counts)] <- "duplicate"
  expect_error(repliapa_validate_pseudobulk(duplicated_counts, duplicated_data),
               "at most once")
})

test_that("pseudobulk aggregation conserves every PAS", {
  counts <- matrix(c(1, 2, 3, 4, 0, 1, 0, 1), nrow = 2,
                   dimnames = list(c("p1", "p2"), paste0("c", 1:4)))
  meta <- data.frame(
    cell_line = "L", replicate = c("R1", "R1", "R2", "R2"),
    guide = c("g1", "g1", "n1", "n1"), target = c("T", "T", "NT", "NT"),
    ntc = c(FALSE, FALSE, TRUE, TRUE), row.names = colnames(counts))
  out <- repliapa_pseudobulk(counts, meta)
  expect_equal(as.numeric(Matrix::rowSums(out$counts)), as.numeric(rowSums(counts)))
  expect_equal(out$sample_data$cell_n, c(2L, 2L))
  expect_equal(ncol(out$counts), 2L)
})
