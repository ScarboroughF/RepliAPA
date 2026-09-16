test_that("independent replicate fits can be compared without cell splitting", {
  x <- repliapa_example()
  r1 <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1", replicates = "R1")
  r2 <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1", replicates = "R2")
  out <- repliapa_reproducibility(r1, r2)
  expect_true(out$effect_vector_cosine[out$gene == "GENE_A"] > 0.99)
  expect_equal(attr(out, "summary")$tested_in_both, 2)
  expect_equal(attr(out, "summary")$rank_spearman, 1)
})

test_that("different targets cannot be compared as replicate validation", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  other <- fit
  other$target <- "T2"
  expect_error(repliapa_reproducibility(fit, other), "same target")
})
