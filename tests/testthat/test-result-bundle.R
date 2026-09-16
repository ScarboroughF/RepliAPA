test_that("result bundle preserves fitted evidence and accounting", {
  x <- repliapa_example()
  batch <- repliapa_fit_targets(
    x$counts, x$sample_data, x$pas_data,
    targets = c("T1", "missing_target")
  )
  bundle <- repliapa_result_bundle(batch)

  expect_s3_class(bundle, "repliapa_result_bundle")
  expect_identical(bundle$tables$tested_gene_accounting$target,
                   c("T1", "missing_target"))
  expect_identical(bundle$tables$tested_gene_accounting$status,
                   c("FITTED", "UNSUPPORTED"))
  expect_equal(bundle$tables$tested_gene_accounting$native_gene_n, c(3L, 0L))
  expect_equal(bundle$tables$tested_gene_accounting$tested_gene_n, c(2L, 0L))
  expect_true(all(bundle$tables$target_effects$target == "T1"))
  expect_true(all(bundle$tables$pas_effects$target == "T1"))
  expect_true(all(bundle$tables$guide_effects$target == "T1"))
  expect_false(bundle$tables$provenance$inference_changed)
  expect_identical(bundle$tables$provenance$peak_memory_status,
                   "NOT_INSTRUMENTED")
})

test_that("bundle can include verified freeze and outcome-free input audit", {
  x <- repliapa_example()
  x$sample_data$cell_n <- 30L
  audit <- repliapa_audit_input(
    x$counts, x$sample_data, x$pas_data, "R1", min_polya_count = 1L
  )
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  freeze <- make_test_freeze()
  bundle <- repliapa_result_bundle(fit, freeze = freeze, input_audit = audit)

  expect_identical(bundle$tables$provenance$freeze_fingerprint,
                   freeze$fingerprint)
  expect_true(bundle$tables$provenance$input_audit_included)
  expect_equal(nrow(bundle$tables$sample_unit_manifest), ncol(x$counts))
  expect_equal(nrow(bundle$tables$pas_adequacy), nrow(x$counts))
  expect_false(audit$settings$outcome_used)

  changed <- freeze
  changed$config$pseudocount <- 1
  expect_error(repliapa_result_bundle(fit, freeze = changed), "fingerprint")
})

test_that("written bundles are deterministic and tamper evident", {
  x <- repliapa_example()
  bundle <- repliapa_result_bundle(
    repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  )
  first <- tempfile("repliapa_bundle_a_")
  second <- tempfile("repliapa_bundle_b_")
  on.exit(unlink(c(first, second), recursive = TRUE, force = TRUE), add = TRUE)
  expect_invisible(repliapa_write_bundle(bundle, first))
  expect_invisible(repliapa_write_bundle(bundle, second))
  expect_true(repliapa_verify_bundle(first))

  checks_a <- readLines(file.path(first, "CHECKSUMS.sha256"))
  checks_b <- readLines(file.path(second, "CHECKSUMS.sha256"))
  expect_identical(checks_a, checks_b)
  expect_true(file.exists(file.path(first, "TARGET_EFFECTS.tsv")))
  expect_true(file.exists(file.path(first, "TESTED_GENE_ACCOUNTING.tsv")))
  expect_error(repliapa_write_bundle(bundle, first), "not empty")

  write("tampered", file.path(first, "TARGET_EFFECTS.tsv"), append = TRUE)
  expect_error(repliapa_verify_bundle(first), "SHA256")
})

test_that("overwrite fails closed around unrelated files", {
  x <- repliapa_example()
  unsupported <- repliapa_fit_targets(
    x$counts, x$sample_data, x$pas_data, targets = "missing_target"
  )
  bundle <- repliapa_result_bundle(unsupported)
  path <- tempfile("repliapa_bundle_unrelated_")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  writeLines("user file", file.path(path, "KEEP.txt"))
  expect_error(repliapa_write_bundle(bundle, path, overwrite = TRUE),
               "unrelated files")

  clean <- tempfile("repliapa_bundle_emptyfit_")
  on.exit(unlink(clean, recursive = TRUE, force = TRUE), add = TRUE)
  expect_invisible(repliapa_write_bundle(bundle, clean))
  expect_true(repliapa_verify_bundle(clean))
  target_effects <- utils::read.delim(
    file.path(clean, "TARGET_EFFECTS.tsv"), check.names = FALSE
  )
  expect_equal(nrow(target_effects), 0L)
})

test_that("measured compute provenance is validated and retained", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  benchmark <- data.frame(
    method = c("REPLIAPA_ALL_PAS", "DRIMSEQ_ALL_PAS"),
    repetitions = c(3L, 3L),
    elapsed_seconds_median = c(1.2, 2.4),
    peak_rss_mb_median = c(100, 300),
    input_sha256 = rep(paste(rep("a", 64), collapse = ""), 2),
    elapsed_scope = "MODEL_EXCLUDES_INPUT_LOAD",
    peak_rss_scope = "WHOLE_ISOLATED_PROCESS",
    peak_rss_source = "LINUX_PROC_VM_HWM",
    stringsAsFactors = FALSE
  )
  bundle <- repliapa_result_bundle(fit, compute_benchmark = benchmark)

  expect_identical(bundle$tables$compute_benchmark$method,
                   c("DRIMSEQ_ALL_PAS", "REPLIAPA_ALL_PAS"))
  expect_true(bundle$tables$provenance$compute_benchmark_included)
  expect_identical(bundle$tables$provenance$runtime_status,
                   "MEASURED_ISOLATED_PROCESS_MEDIAN")
  expect_identical(bundle$tables$provenance$peak_memory_status,
                   "MEASURED_KERNEL_VM_HWM")

  duplicate <- rbind(benchmark, benchmark[1, ])
  expect_error(repliapa_result_bundle(fit, compute_benchmark = duplicate), "unique")
  invalid <- benchmark
  invalid$peak_rss_mb_median[1] <- NA_real_
  expect_error(repliapa_result_bundle(fit, compute_benchmark = invalid),
               "positive numbers")
  invalid_hash <- benchmark
  invalid_hash$input_sha256 <- "not-a-hash"
  expect_error(repliapa_result_bundle(fit, compute_benchmark = invalid_hash),
               "SHA256")
})
