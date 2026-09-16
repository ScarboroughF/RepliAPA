test_that("replicate reproducibility rejects overlapping biological units", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1", replicates = "R1")
  expect_error(
    repliapa_reproducibility(fit, fit),
    "disjoint biological replicates"
  )
})

test_that("batch fitting fails on global configuration errors", {
  x <- repliapa_example()
  expect_error(
    repliapa_fit_targets(
      x$counts, x$sample_data, x$pas_data,
      targets = c("T1", "missing"), pseudocount = NA_real_
    ),
    "pseudocount"
  )
  expect_error(
    repliapa_fit_targets(
      x$counts, x$sample_data, x$pas_data, targets = "T1", strict = NA
    ),
    "strict"
  )
})

test_that("integer controls are not silently truncated", {
  x <- repliapa_example()
  expect_error(
    repliapa_ntc_calibration(
      x$counts, x$sample_data, x$pas_data,
      n_contrasts = 1.5, guide_sizes = 2L
    ),
    "n_contrasts"
  )

  effects <- data.frame(
    target = rep(c("A", "B", "C"), each = 4),
    replicate = "R1", guide = rep(paste0("g", 1:6), each = 2),
    gene = rep("G", 12), pas_id = rep(c("p1", "p2"), 6),
    effect = rep(c(-0.1, 0.1), 6), stringsAsFactors = FALSE
  )
  guide_data <- data.frame(
    replicate = "R1", target = rep(c("A", "B", "C"), each = 2),
    guide = paste0("g", 1:6), depth = seq_len(6), stringsAsFactors = FALSE
  )
  expect_error(
    repliapa_guide_reproducibility(
      effects, guide_data, "depth", n_boot = 1.5
    ),
    "n_boot"
  )
})

test_that("named fit collections cannot misroute target effects", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  expect_error(
    repliapa_compare_contexts(
      list(wrong_name = fit), list(T1 = fit)
    ),
    "names identical"
  )
})

test_that("gene-expression classification validates thresholds and identities", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  ge <- data.frame(gene = fit$gene_table$gene, log2FC = 0)
  expect_error(repliapa_classify_effects(fit, ge, apa_fdr = NA_real_), "apa_fdr")
  expect_error(repliapa_classify_effects(fit, ge, ge_abs_log2fc = -1), "ge_abs_log2fc")
  duplicate <- rbind(ge, ge[1, ])
  expect_error(repliapa_classify_effects(fit, duplicate), "unique")
})

test_that("replicate-effect schema accepts realistic labels", {
  x <- repliapa_example()
  x$sample_data$replicate <- sub("R", "HEK_rep_", x$sample_data$replicate)
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  bundle <- repliapa_result_bundle(fit)
  expect_true(all(c("effect_HEK_rep_1", "effect_HEK_rep_2") %in%
                    names(bundle$tables$pas_effects)))
  expect_true(repliapa_verify_contract(bundle = bundle)$passed)
})
