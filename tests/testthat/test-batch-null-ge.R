test_that("batch fitting retains unsupported-target accounting", {
  x <- repliapa_example()
  batch <- repliapa_fit_targets(x$counts, x$sample_data, x$pas_data,
                                targets = c("T1", "MISSING"))
  expect_s3_class(batch, "repliapa_fit_list")
  expect_identical(batch$accounting$status, c("FITTED", "UNSUPPORTED"))
  expect_named(batch$fits, "T1")
  expect_equal(batch$accounting$tested_gene_n[1], 2)
  expect_error(repliapa_fit_targets(x$counts, x$sample_data, x$pas_data,
                                    targets = "MISSING", strict = TRUE),
               "At least two target guides")
})

test_that("predeclared NTC contrasts preserve whole guides", {
  x <- repliapa_example()
  contrasts <- data.frame(
    contrast_id = c("C1", "C2"),
    pseudo_target_guides = c("NT_g1;NT_g2", "NT_g3;NT_g4"),
    stringsAsFactors = FALSE
  )
  out <- repliapa_ntc_calibration(x$counts, x$sample_data, x$pas_data,
                                  contrasts = contrasts)
  expect_s3_class(out, "repliapa_null_calibration")
  expect_equal(out$summary$contrasts, 2)
  expect_equal(unique(out$results$target_guide_n[out$results$tested]), 2)
  expect_true(all(out$results$permutation_n[out$results$tested] == 6))
  expect_identical(out$contrasts$contrast_id, c("C1", "C2"))
})

test_that("generated NTC contrasts are deterministic without changing caller RNG", {
  x <- repliapa_example()
  set.seed(99)
  before <- .Random.seed
  a <- repliapa_ntc_calibration(x$counts, x$sample_data, x$pas_data,
                                guide_sizes = 2, n_contrasts = 4, seed = 12)
  after <- .Random.seed
  b <- repliapa_ntc_calibration(x$counts, x$sample_data, x$pas_data,
                                guide_sizes = 2, n_contrasts = 4, seed = 12)
  expect_identical(before, after)
  expect_identical(a$contrasts, b$contrasts)
  expect_identical(a$results, b$results)
})

test_that("gene expression remains separate from APA composition", {
  x <- repliapa_example()
  sample_data <- x$sample_data
  sample_data$cell_n <- 10
  rna <- matrix(10L, nrow = 3, ncol = nrow(sample_data),
                dimnames = list(c("GENE_A", "GENE_B", "GENE_C"), rownames(sample_data)))
  rna["GENE_B", !sample_data$ntc] <- 20L
  ge <- repliapa_gene_expression_effect(rna, sample_data, "T1")
  expect_equal(ge$log2FC[ge$gene == "GENE_A"], 0)
  expect_equal(ge$log2FC[ge$gene == "GENE_B"], log2(2.5 / 1.5))

  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  fit$gene_table$q_value[fit$gene_table$gene == "GENE_A"] <- 0.01
  classified <- repliapa_classify_effects(fit, ge, ge_abs_log2fc = 0.25)
  expect_identical(classified$effect_class[classified$gene == "GENE_A"], "APA_ONLY")
  expect_identical(classified$effect_class[classified$gene == "GENE_B"], "GE_ONLY")
  expect_identical(classified$effect_class[classified$gene == "GENE_C"], "NEITHER")
})
