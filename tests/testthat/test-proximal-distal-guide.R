test_that("proximal/distal collapse follows transcript direction and conserves counts", {
  x <- repliapa_example()
  x$pas_data$strand[x$pas_data$gene == "GENE_A"] <- "-"
  pd <- repliapa_collapse_proximal_distal(x$counts, x$pas_data)
  expect_s3_class(pd, "repliapa_pd_data")
  expect_equal(nrow(pd$counts), 4)
  expect_equal(Matrix::colSums(pd$counts),
               Matrix::colSums(x$counts[x$pas_data$gene %in% c("GENE_A", "GENE_B"), ]))
  a <- pd$membership[pd$membership$gene == "GENE_A", ]
  expect_identical(a$pas_id, c("pas3", "pas2", "pas1"))
  expect_identical(a$bin, c("proximal", "proximal", "distal"))
  expect_identical(pd$accounting$reason[pd$accounting$gene == "GENE_C"],
                   "SINGLE_PAS_GENE")
})

test_that("two-bin fitting is directly comparable with all-PAS fitting", {
  x <- repliapa_example()
  all_pas <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  pd <- repliapa_fit_proximal_distal(x$counts, x$sample_data, x$pas_data, "T1")
  comparison <- repliapa_compare_allpas_pd(all_pas, pd)
  expect_s3_class(pd, "repliapa_pd_fit")
  expect_true(all(comparison$matched_tested[comparison$gene != "GENE_C"]))
  expect_equal(attr(comparison, "accounting")$matched_tested_gene_n, 2)
  expect_error(repliapa_compare_allpas_pd(all_pas,
    repliapa_fit_proximal_distal(x$counts, x$sample_data, x$pas_data, "T1",
                                replicates = "R1")), "same target and biological replicates")
})

make_guide_reproducibility_example <- function() {
  guide_data <- do.call(rbind, lapply(c("R1", "R2"), function(rep_value) {
    data.frame(
      replicate = rep_value,
      target = rep(c("A", "B", "C"), each = 2),
      guide = paste0(rep(c("A", "B", "C"), each = 2), "_g", rep(1:2, 3)),
      cell_n = c(100, 102, 80, 83, 120, 118),
      depth = c(1000, 1010, 800, 810, 1200, 1190),
      stringsAsFactors = FALSE
    )
  }))
  vectors <- list(
    A_g1 = c(.20, -.10, -.10), A_g2 = c(.18, -.08, -.10),
    B_g1 = c(-.20, .10, .10), B_g2 = c(-.18, .08, .10),
    C_g1 = c(.05, .10, -.15), C_g2 = c(-.05, -.10, .15)
  )
  effects <- do.call(rbind, lapply(seq_len(nrow(guide_data)), function(i) {
    row <- guide_data[i, ]
    data.frame(target = row$target, replicate = row$replicate, guide = row$guide,
               gene = "G", pas_id = paste0("p", 1:3),
               effect = vectors[[row$guide]], stringsAsFactors = FALSE)
  }))
  list(effects = effects, guide_data = guide_data)
}

test_that("matched guide reproducibility preserves replicate and different targets", {
  x <- make_guide_reproducibility_example()
  out <- repliapa_guide_reproducibility(
    x$effects, x$guide_data, c("cell_n", "depth"), n_boot = 20, seed = 7
  )
  expect_s3_class(out, "repliapa_guide_reproducibility")
  expect_equal(nrow(out$pair_summary), 6)
  expect_true(all(out$matching$null_target_1 != out$matching$target))
  expect_true(all(out$matching$null_target_2 != out$matching$target))
  expect_true(all(out$matching$null_target_1 != out$matching$null_target_2))
  expect_true(mean(out$pair_summary$mean_difference) > 0)
  expect_equal(out$bootstrap$resamples, rep(20L, 2))
})

test_that("guide effects are collected losslessly from batch fits", {
  x <- repliapa_example()
  batch <- repliapa_fit_targets(x$counts, x$sample_data, x$pas_data, "T1")
  effects <- repliapa_collect_guide_effects(batch)
  expect_named(effects, c("target", "replicate", "guide", "gene", "pas_id", "effect"))
  expect_identical(unique(effects$target), "T1")
  expect_equal(nrow(effects), nrow(batch$fits$T1$guide_effects))
})

test_that("guide matching and bootstrap are deterministic without changing RNG", {
  x <- make_guide_reproducibility_example()
  set.seed(91)
  before <- .Random.seed
  first <- repliapa_guide_reproducibility(
    x$effects, x$guide_data, c("cell_n", "depth"), n_boot = 10, seed = 4
  )
  after <- .Random.seed
  second <- repliapa_guide_reproducibility(
    x$effects, x$guide_data, c("cell_n", "depth"), n_boot = 10, seed = 4
  )
  expect_identical(before, after)
  expect_identical(first$matching, second$matching)
  expect_identical(first$bootstrap, second$bootstrap)
})
