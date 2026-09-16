test_that("cross-context categories use independently fitted effects", {
  x <- repliapa_example()
  primary <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  external <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  primary$gene_table$q_value[primary$gene_table$gene %in% c("GENE_A", "GENE_B")] <- 0.01
  external$gene_table$q_value[external$gene_table$gene %in% c("GENE_A", "GENE_B")] <- 0.01
  external$pas_effects$effect[external$pas_effects$gene == "GENE_B"] <-
    -external$pas_effects$effect[external$pas_effects$gene == "GENE_B"]

  comparison <- repliapa_compare_contexts(primary, external,
                                           primary_label = "HEK",
                                           external_label = "K562")
  expect_s3_class(comparison, "repliapa_context_comparison")
  expect_identical(comparison$results$category[comparison$results$gene == "GENE_A"],
                   "SHARED_EFFECT")
  expect_identical(comparison$results$category[comparison$results$gene == "GENE_B"],
                   "OPPOSITE_EFFECT")
  expect_identical(comparison$results$category[comparison$results$gene == "GENE_C"],
                   "UNSUPPORTED_GENE_CONTEXT")
  expect_equal(comparison$target_summary$comparable_gene_n, 2)
})

test_that("context-specific effects are not called shared", {
  x <- repliapa_example()
  primary <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  external <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  primary$gene_table$q_value[primary$gene_table$gene == "GENE_A"] <- 0.01
  external$gene_table$q_value[external$gene_table$gene == "GENE_A"] <- 0.5
  comparison <- repliapa_compare_contexts(primary, external)
  expect_identical(comparison$results$category[comparison$results$gene == "GENE_A"],
                   "CELL_LINE_SPECIFIC_EFFECT")
})

test_that("unsupported targets remain explicit", {
  x <- repliapa_example()
  fitted <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  primary_only <- fitted
  primary_only$target <- "P_ONLY"
  external_only <- fitted
  external_only$target <- "E_ONLY"
  primary <- structure(list(
    fits = list(T1 = fitted, P_ONLY = primary_only),
    accounting = data.frame(target = c("T1", "P_ONLY", "E_ONLY", "BOTH_MISSING"),
      status = c("FITTED", "FITTED", "UNSUPPORTED", "UNSUPPORTED"),
      reason = c("NONE", "NONE", "LOW_GUIDE_SUPPORT", "LOW_GUIDE_SUPPORT"))
  ), class = "repliapa_fit_list")
  external <- structure(list(
    fits = list(T1 = fitted, E_ONLY = external_only),
    accounting = data.frame(target = c("T1", "P_ONLY", "E_ONLY", "BOTH_MISSING"),
      status = c("FITTED", "UNSUPPORTED", "FITTED", "UNSUPPORTED"),
      reason = c("NONE", "LOW_GUIDE_SUPPORT", "NONE", "LOW_GUIDE_SUPPORT"))
  ), class = "repliapa_fit_list")
  comparison <- repliapa_compare_contexts(
    primary, external, targets = c("T1", "P_ONLY", "E_ONLY", "BOTH_MISSING")
  )
  expect_identical(comparison$target_accounting$comparison_status,
    c("SUPPORTED_BOTH", "EXTERNAL_UNSUPPORTED", "PRIMARY_UNSUPPORTED", "BOTH_UNSUPPORTED"))
  expect_equal(nrow(comparison$results), 3)
})

test_that("context labels and target identities are guarded", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  expect_error(repliapa_compare_contexts(fit, fit, "same", "same"),
               "distinct syntactic")
  duplicate <- list(first = fit, second = fit)
  expect_error(repliapa_compare_contexts(duplicate, fit), "unique")
})
