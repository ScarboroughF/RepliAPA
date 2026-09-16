test_that("v1-candidate API contract is fingerprinted and self-consistent", {
  contract <- repliapa_api_contract()
  verification <- repliapa_verify_contract(contract)

  expect_s3_class(contract, "repliapa_api_contract")
  expect_identical(contract$contract_version, "1.0.0-rc4")
  expect_identical(contract$candidate_api_version, "1.0.0")
  expect_length(contract$functions, 40L)
  expect_length(contract$fingerprint, 1L)
  expect_match(contract$fingerprint, "^[0-9a-f]{64}$")
  expect_s3_class(verification, "repliapa_contract_verification")
  expect_true(verification$passed)
  expect_true(all(verification$checks$passed))
})

test_that("contract verification detects payload and signature tampering", {
  contract <- repliapa_api_contract()
  changed <- contract
  changed$compatibility_policy$inference <- "results may change"
  verification <- repliapa_verify_contract(changed)
  expect_false(verification$passed)
  expect_false(verification$checks$passed[
    verification$checks$check == "CONTRACT_FINGERPRINT"
  ])

  changed <- contract
  changed$functions[["repliapa_fit"]] <- paste0(
    changed$functions[["repliapa_fit"]], "|new_argument=TRUE"
  )
  verification <- repliapa_verify_contract(changed)
  expect_false(verification$passed)
  expect_false(verification$checks$passed[
    verification$checks$check == "SIGNATURE:repliapa_fit"
  ])
})

test_that("contract validates canonical and optional bundle schemas", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  benchmark <- data.frame(
    method = "REPLIAPA_ALL_PAS", repetitions = 1L,
    elapsed_seconds_median = 1, peak_rss_mb_median = 100,
    input_sha256 = paste(rep("a", 64), collapse = ""),
    elapsed_scope = "MODEL_EXCLUDES_INPUT_LOAD",
    peak_rss_scope = "WHOLE_ISOLATED_PROCESS",
    peak_rss_source = "LINUX_PROC_VM_HWM",
    stringsAsFactors = FALSE
  )
  bundle <- repliapa_result_bundle(fit, compute_benchmark = benchmark)
  verification <- repliapa_verify_contract(bundle = bundle)
  expect_true(verification$passed)

  missing <- bundle
  missing$tables$target_effects$tested <- NULL
  verification <- repliapa_verify_contract(bundle = missing)
  expect_false(verification$passed)
  expect_false(verification$checks$passed[
    verification$checks$check == "TABLE_SCHEMA:target_effects"
  ])

  unknown <- bundle
  unknown$tables$target_effects$unfrozen_result <- 1
  verification <- repliapa_verify_contract(bundle = unknown)
  expect_false(verification$passed)

  declared_extension <- bundle
  declared_extension$tables$pas_effects$effect_R3 <- 0
  expect_true(repliapa_verify_contract(bundle = declared_extension)$passed)
})

test_that("contract records required components of principal result objects", {
  contract <- repliapa_api_contract()
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  pd <- repliapa_collapse_proximal_distal(x$counts, x$pas_data)
  batch <- repliapa_fit_targets(
    x$counts, x$sample_data, x$pas_data, targets = c("T1", "missing")
  )
  drim <- repliapa_drimseq_input(
    x$counts, x$sample_data, x$pas_data, target = "T1"
  )
  objects <- list(
    repliapa_input = repliapa_read_pas_counts(
      x$counts, x$sample_data, x$pas_data
    ),
    repliapa_fit = fit,
    repliapa_pd_data = pd,
    repliapa_fit_list = batch,
    repliapa_drimseq_input = drim,
    repliapa_freeze = make_test_freeze(),
    repliapa_result_bundle = repliapa_result_bundle(fit)
  )
  for (class_name in names(objects)) {
    expect_true(all(contract$classes[[class_name]] %in% names(objects[[class_name]])),
                info = class_name)
  }
})
