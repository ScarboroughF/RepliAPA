test_that("DRIMSeq export is deterministic and conserves selected counts", {
  x <- repliapa_example()
  out <- repliapa_drimseq_input(
    x$counts, x$sample_data, x$pas_data, target = "T1",
    pas_ids = c("pas1", "pas2", "pas3", "pas4", "pas5")
  )

  expect_s3_class(out, "repliapa_drimseq_input")
  expect_identical(names(out$counts)[1:2], c("gene_id", "feature_id"))
  expect_identical(out$counts$feature_id, paste0("pas", 1:5))
  expect_identical(out$samples$sample_id, colnames(out$counts)[-(1:2)])
  expect_identical(out$samples$group, rep(c(rep(1L, 3), rep(0L, 4)), 2))
  expect_equal(sum(as.matrix(out$counts[, -(1:2)])), sum(x$counts[1:5, ]))
  expect_true(out$accounting$count_conserved)
  expect_identical(out$accounting$export_filtering, "NONE")
  expect_equal(out$accounting$exported_pas_n, 5L)
  expect_equal(as.numeric(out$design[, "group"]), out$samples$group)

  again <- repliapa_drimseq_input(
    x$counts, x$sample_data, x$pas_data, target = "T1",
    pas_ids = c("pas1", "pas2", "pas3", "pas4", "pas5")
  )
  expect_identical(out, again)
})

test_that("DRIMSeq export keeps only whole guides complete across replicates", {
  x <- repliapa_example()
  missing_unit <- rownames(x$sample_data)[
    x$sample_data$replicate == "R2" & x$sample_data$guide == "T1_g3"
  ]
  keep <- setdiff(colnames(x$counts), missing_unit)
  out <- repliapa_drimseq_input(
    x$counts[, keep], x$sample_data[keep, ], x$pas_data, target = "T1"
  )

  expect_identical(out$target_guides, c("T1_g1", "T1_g2"))
  expect_false(any(out$samples$guide == "T1_g3"))
  expect_equal(as.vector(table(out$samples$guide)), rep(2L, 6))
  expect_error(
    repliapa_drimseq_input(x$counts, x$sample_data, x$pas_data, "T1",
                           pas_ids = c("pas1", "missing")),
    "unavailable"
  )
})

test_that("baseline comparison reports native and matched tested sets", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  fit$gene_table$q_value[fit$gene_table$gene == "GENE_A"] <- 0.01
  fit$gene_table$q_value[fit$gene_table$gene == "GENE_B"] <- 0.5

  baseline <- data.frame(
    gene_id = c("GENE_A", "GENE_B", "GENE_D"),
    pvalue = c(0.01, NA, 0.03),
    adj = c(0.02, NA, 0.04),
    tested = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  comparison <- repliapa_compare_baseline(
    fit, baseline, baseline_q = "adj", baseline_tested = "tested"
  )

  expect_s3_class(comparison, "repliapa_baseline_comparison")
  expect_identical(comparison$native_summary$tested_gene_n, c(2L, 2L))
  expect_identical(comparison$native_summary$untested_gene_n, c(1L, 1L))
  expect_equal(comparison$native_summary$discovery_n, c(1, 2))
  expect_identical(comparison$matched_summary$matched_tested_gene_n, c(1L, 1L))
  expect_equal(comparison$matched_summary$discovery_n, c(1, 1))
  expect_equal(comparison$agreement$discovery_overlap_n, 1)
  expect_true(comparison$gene_table$baseline_present[
    comparison$gene_table$gene == "GENE_D"
  ])
  expect_false(comparison$gene_table$repliapa_present[
    comparison$gene_table$gene == "GENE_D"
  ])
})

test_that("baseline BH correction uses only native tested genes", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  baseline <- data.frame(
    gene_id = c("GENE_A", "GENE_B", "GENE_C"),
    pvalue = c(0.01, NA, 0.04),
    tested = c(TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  comparison <- repliapa_compare_baseline(
    fit, baseline, baseline_tested = "tested"
  )
  observed <- comparison$gene_table$baseline_q_value[
    match(baseline$gene_id, comparison$gene_table$gene)
  ]
  expect_equal(observed, c(0.02, NA, 0.04))

  duplicate <- rbind(baseline, baseline[1, ])
  expect_error(repliapa_compare_baseline(fit, duplicate), "unique")
  nonnumeric <- baseline
  nonnumeric$pvalue <- as.character(nonnumeric$pvalue)
  expect_error(repliapa_compare_baseline(fit, nonnumeric), "must be numeric")
})
