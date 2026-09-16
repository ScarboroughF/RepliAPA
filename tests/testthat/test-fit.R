test_that("admissible permutations preserve guide identity across replicates", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  expect_s3_class(fit, "repliapa_fit")
  expect_equal(unique(fit$gene_table$permutation_n[fit$gene_table$tested]), 35)
  expect_equal(fit$gene_table$p_value[fit$gene_table$gene == "GENE_A"], 1 / 35)
  expect_identical(fit$replicates, c("R1", "R2"))
  expect_identical(fit$target_guides, paste0("T1_g", 1:3))
})

test_that("PAS effects are compositional and guide evidence is retained", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  sums <- tapply(fit$pas_effects$effect, fit$pas_effects$gene, sum)
  expect_equal(as.numeric(sums), rep(0, length(sums)), tolerance = 1e-12)
  expect_true(all(c("replicate", "guide", "effect") %in% names(fit$guide_effects)))
  expect_equal(unique(fit$heterogeneity$guide_n), 3)
})

test_that("filtering and abstention are explicit", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  single <- fit$gene_table[fit$gene_table$gene == "GENE_C", ]
  expect_false(single$tested)
  expect_true(single$filtered)
  expect_false(single$abstained)
  expect_identical(single$reason, "SINGLE_PAS_GENE")

  limited <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1",
                          max_combinations = 20)
  expect_true(all(limited$gene_table$abstained[limited$gene_table$pas_n >= 2]))
  expect_true(all(limited$gene_table$reason[limited$gene_table$pas_n >= 2] ==
                    "TOO_MANY_ADMISSIBLE_PERMUTATIONS"))
})

test_that("sample units, not cells, are required by fit", {
  x <- repliapa_example()
  cell_like <- rbind(x$sample_data, x$sample_data[1, , drop = FALSE])
  rownames(cell_like)[nrow(cell_like)] <- "extra_cell"
  counts <- cbind(x$counts, x$counts[, 1])
  colnames(counts)[ncol(counts)] <- "extra_cell"
  expect_error(repliapa_fit(counts, cell_like, x$pas_data, "T1"), "at most once")
})

test_that("fitting is deterministic and does not split guide labels", {
  x <- repliapa_example()
  a <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  b <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  expect_identical(a$gene_table, b$gene_table)
  expect_identical(a$pas_effects, b$pas_effects)
})
