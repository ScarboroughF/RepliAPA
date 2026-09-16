test_that("input audit separates contract, support, and adequacy", {
  x <- repliapa_example()
  x$sample_data$cell_n <- 30L
  audit <- repliapa_audit_input(
    x$counts, x$sample_data, x$pas_data,
    development_replicates = "R1", min_polya_count = 1L
  )

  expect_s3_class(audit, "repliapa_input_audit")
  expect_true(all(audit$contract_checks$passed))
  expect_equal(audit$summary$sample_unit_n, ncol(x$counts))
  expect_equal(audit$summary$qualified_sample_unit_n, ncol(x$counts))
  expect_equal(audit$summary$eligible_target_n, 1L)
  expect_true(audit$target_support$target_eligible[audit$target_support$target == "T1"])
  expect_identical(audit$target_support$exclusion_reason[
    audit$target_support$target == "NT"
  ], "NEGATIVE_CONTROL")
  expect_equal(nrow(audit$pas_adequacy), nrow(x$counts))
  expect_equal(nrow(audit$gene_adequacy), 3L)
  expect_identical(audit$gene_adequacy$exclusion_reason[
    audit$gene_adequacy$gene == "GENE_C"
  ], "SINGLE_PAS_GENE")
  expect_false(audit$settings$outcome_used)
})

test_that("PAS tiers use only predeclared development replicates", {
  x <- repliapa_example()
  x$sample_data$cell_n <- 30L
  a <- repliapa_audit_input(
    x$counts, x$sample_data, x$pas_data, "R1", min_polya_count = 1L
  )
  changed <- x$counts
  changed[, x$sample_data$replicate == "R2"] <-
    changed[, x$sample_data$replicate == "R2"] * 100L
  b <- repliapa_audit_input(
    changed, x$sample_data, x$pas_data, "R1", min_polya_count = 1L
  )

  expect_identical(a$pas_adequacy$adequacy_tier, b$pas_adequacy$adequacy_tier)
  expect_identical(a$pas_adequacy$development_total_count,
                   b$pas_adequacy$development_total_count)
  expect_false(identical(a$pas_adequacy$all_sample_total_count,
                         b$pas_adequacy$all_sample_total_count))
})

test_that("target eligibility requires independent guides in every replicate", {
  x <- repliapa_example()
  x$sample_data$cell_n <- 30L
  low <- x$sample_data$replicate == "R2" &
    x$sample_data$guide %in% c("T1_g2", "T1_g3")
  x$sample_data$cell_n[low] <- 1L
  audit <- repliapa_audit_input(
    x$counts, x$sample_data, x$pas_data, "R1", min_polya_count = 1L
  )

  t1 <- audit$target_support[audit$target_support$target == "T1", ]
  expect_false(t1$target_eligible)
  expect_match(t1$exclusion_reason, "R2")
  expect_equal(audit$target_replicate_support$qualified_guide_n[
    audit$target_replicate_support$target == "T1" &
      audit$target_replicate_support$replicate == "R2"
  ], 1L)
})

test_that("exact coordinate conflicts and mapping violations are explicit", {
  x <- repliapa_example()
  x$sample_data$cell_n <- 30L
  x$pas_data$seqnames <- "chr1"
  x$pas_data$start <- x$pas_data$position
  x$pas_data$end <- x$pas_data$position
  x$pas_data[c("pas1", "pas2"), c("start", "end")] <- 100
  audit <- repliapa_audit_input(
    x$counts, x$sample_data, x$pas_data, "R1", min_polya_count = 1L
  )
  expect_true(audit$gene_adequacy$coordinate_conflict[
    audit$gene_adequacy$gene == "GENE_A"
  ])
  expect_identical(audit$gene_adequacy$exclusion_reason[
    audit$gene_adequacy$gene == "GENE_A"
  ], "EXACT_COORDINATE_CONFLICT")

  inconsistent <- x$sample_data
  inconsistent$target[inconsistent$replicate == "R2" &
                        inconsistent$guide == "T1_g1"] <- "OTHER"
  expect_error(
    repliapa_audit_input(x$counts, inconsistent, x$pas_data, "R1",
                         min_polya_count = 1L),
    "multiple targets"
  )
  expect_error(
    repliapa_audit_input(x$counts, x$sample_data, x$pas_data, "missing",
                         min_polya_count = 1L),
    "available replicate"
  )
})
