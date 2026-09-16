test_that("freeze fingerprints are deterministic and detect mutation", {
  first <- make_test_freeze()
  second <- make_test_freeze()
  expect_s3_class(first, "repliapa_freeze")
  expect_identical(first$fingerprint, second$fingerprint)
  expect_true(repliapa_verify_freeze(first))
  first$config$pseudocount <- 1
  expect_error(repliapa_verify_freeze(first), "fingerprint mismatch")
})

test_that("recorded governance source files are checked", {
  source <- tempfile()
  writeLines("frozen: true", source)
  freeze <- repliapa_create_freeze(
    config = list(model = "TEST"), pas_ids = c("p1", "p2"),
    primary_metrics = "replication",
    role_registry = data.frame(dataset = "development", role = "DEVELOPMENT"),
    source_files = c(config = source), frozen_on = "2026-08-11"
  )
  expect_true(repliapa_verify_freeze(freeze, check_source_files = TRUE))
  writeLines("frozen: false", source)
  expect_error(repliapa_verify_freeze(freeze, check_source_files = TRUE),
               "source file changed")
})

test_that("data roles authorize allowed actions and fail closed", {
  freeze <- make_test_freeze()
  expect_true(repliapa_authorize(freeze, "DEV", "tune_model"))
  expect_true(repliapa_authorize(freeze, "REP2", "apply_frozen_model"))
  expect_true(repliapa_authorize(freeze, "K562", "read_outcome"))
  expect_true(repliapa_authorize(freeze, "MPRA", "orthogonal_validation"))
  expect_true(repliapa_authorize(freeze, "PATIENT", "metadata_audit"))
  expect_true(repliapa_authorize(freeze, "NTC", "null_calibration"))
  expect_error(repliapa_authorize(freeze, "REP2", "tune_filter"), "prohibited")
  expect_error(repliapa_authorize(freeze, "K562", "change_metric"), "prohibited")
  expect_error(repliapa_authorize(freeze, "PATIENT", "fit_model"), "prohibited")
  expect_error(repliapa_authorize(freeze, "NTC", "split_guide"), "prohibited")
})

test_that("external application records freeze provenance", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  out <- repliapa_apply_external(
    fit, fit, make_test_freeze(), "K562",
    primary_label = "HEK", external_label = "K562"
  )
  expect_s3_class(out, "repliapa_context_comparison")
  expect_identical(out$external_registry_id, "K562")
  expect_identical(out$freeze_fingerprint, make_test_freeze()$fingerprint)
  expect_error(repliapa_apply_external(fit, fit, make_test_freeze(), "PATIENT"),
               "prohibited")
})
