test_that("publication theme and palettes are stable ggplot components", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(repliapa_theme(), "theme")
  expect_identical(names(repliapa_palette("effect")),
                   c("Increased usage", "Decreased usage", "No shift"))
  expect_length(repliapa_palette("pas", 5), 5L)
})

test_that("fit-level visualizations preserve PAS, guide, and replicate evidence", {
  skip_if_not_installed("ggplot2")
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")

  effect <- repliapa_plot_pas_effect(fit, "GENE_A", x$pas_data)
  guide <- repliapa_plot_guide_evidence(fit, "GENE_A")

  expect_s3_class(effect, "ggplot")
  expect_s3_class(guide, "ggplot")
  expect_equal(nrow(effect$data), 3L)
  expect_equal(nrow(guide$data), 18L)
  expect_silent(ggplot2::ggplot_build(effect))
  expect_silent(ggplot2::ggplot_build(guide))
})

test_that("validation visualizations accept native RepliAPA results", {
  skip_if_not_installed("ggplot2")
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  development <- repliapa_fit(
    x$counts, x$sample_data, x$pas_data, "T1", replicates = "R1"
  )
  validation <- repliapa_fit(
    x$counts, x$sample_data, x$pas_data, "T1", replicates = "R2"
  )
  reproducibility <- repliapa_reproducibility(development, validation)
  null <- repliapa_ntc_calibration(
    x$counts, x$sample_data, x$pas_data, guide_sizes = 2L, n_contrasts = 2L
  )
  pd <- repliapa_fit_proximal_distal(x$counts, x$sample_data, x$pas_data, "T1")
  allpas <- repliapa_compare_allpas_pd(fit, pd)
  context <- repliapa_compare_contexts(fit, fit, primary_label = "HEK",
                                       external_label = "K562")
  ge <- data.frame(gene = fit$gene_table$gene,
                   log2FC = seq(-0.4, 0.4, length.out = nrow(fit$gene_table)))
  classification <- repliapa_classify_effects(fit, ge)

  plots <- list(
    repliapa_plot_calibration(null),
    repliapa_plot_reproducibility(reproducibility, label_n = 0L),
    repliapa_plot_allpas(allpas),
    repliapa_plot_context(context),
    repliapa_plot_apa_ge(classification, label_n = 0L)
  )
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(plots, function(plot) {
    !inherits(try(ggplot2::ggplot_build(plot), silent = TRUE), "try-error")
  }, logical(1))))
})

test_that("visualizations fail explicitly for unsupported inputs", {
  skip_if_not_installed("ggplot2")
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  expect_error(repliapa_plot_pas_effect(fit, "MISSING"), "no tested PAS")
  expect_error(repliapa_plot_guide_evidence(fit, "MISSING"), "no tested PAS")
  expect_error(repliapa_plot_calibration(list()), "null_calibration")
  expect_error(repliapa_palette("pas", 0), "n")
})
