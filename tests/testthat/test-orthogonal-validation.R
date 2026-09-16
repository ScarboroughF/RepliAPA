make_orthogonal_example <- function() {
  rows <- expand.grid(
    gene_id = c("GENE_A", "UNMAPPED"),
    replicate = c("rep1", "rep2"),
    perturbation = c("NT", "T1"),
    barcode = 1:2,
    stringsAsFactors = FALSE
  )
  rows$pas_id <- paste0(rows$gene_id, "_construct")
  rows$aim <- "validation"
  rows$subaim <- "wt_site"
  rows$experiment <- "wt"
  rows$total <- 100L
  target_positive <- rows$perturbation == "T1"
  rows$proximal <- ifelse(target_positive, 70L, 40L)
  rows$distal <- rows$total - rows$proximal
  rows
}

test_that("proximal direction is derived from the frozen all-PAS effect", {
  x <- repliapa_example()
  fit <- repliapa_fit(x$counts, x$sample_data, x$pas_data, "T1")
  pd <- repliapa_collapse_proximal_distal(x$counts, x$pas_data)
  directions <- repliapa_proximal_directions(fit, pd)
  expect_true(directions$proximal_direction[directions$gene == "GENE_A"] > 0)
  expect_equal(nrow(directions), 3)
  expect_true(is.na(directions$proximal_direction[directions$gene == "GENE_C"]))
})

test_that("orthogonal outcomes cannot be read before model freeze", {
  effects <- data.frame(target = "T1", gene = "GENE_A", proximal_direction = 0.2,
                        q_value = 0.01, tested = TRUE)
  expect_error(repliapa_orthogonal_validation(
    effects, make_orthogonal_example(), list(), "MPRA"
  ), "repliapa_freeze")
})

test_that("orthogonal validation aggregates rows and reports mapping", {
  effects <- data.frame(target = "T1", gene = "GENE_A", proximal_direction = 0.2,
                        q_value = 0.01, tested = TRUE)
  out <- repliapa_orthogonal_validation(
    effects, make_orthogonal_example(), make_test_freeze(), "MPRA",
    targets = "T1", experiment = "wt"
  )
  expect_s3_class(out, "repliapa_orthogonal_validation")
  expect_equal(out$results$orthogonal_proximal_delta[out$results$gene == "GENE_A"], 0.3)
  expect_true(out$results$direction_concordant[out$results$gene == "GENE_A"])
  expect_identical(out$results$mapping_status[out$results$gene == "UNMAPPED"],
                   "NOT_IN_FROZEN_EFFECT_SET")
  expect_equal(out$accounting$mapped_n, 1)
  expect_equal(out$accounting$not_in_effect_set_n, 1)
  expect_equal(out$target_summary$direction_concordance_all, 1)
})

test_that("insufficient orthogonal replicates abstain", {
  effects <- data.frame(target = "T1", gene = "GENE_A", proximal_direction = 0.2,
                        q_value = 0.01, tested = TRUE)
  data <- make_orthogonal_example()
  data <- data[data$replicate == "rep1" & data$gene_id == "GENE_A", ]
  out <- repliapa_orthogonal_validation(
    effects, data, make_test_freeze(), "MPRA", targets = "T1"
  )
  expect_identical(out$results$mapping_status,
                   "ORTHOGONAL_REPLICATES_INSUFFICIENT")
  expect_true(is.na(out$results$direction_concordant))
  expect_equal(out$accounting$orthogonal_unsupported_n, 1)
})

test_that("constructs with missing replicates retain a stable schema", {
  effects <- data.frame(target = "T1", gene = c("GENE_A", "UNMAPPED"),
                        proximal_direction = c(0.2, 0.1), q_value = 0.01,
                        tested = TRUE)
  data <- make_orthogonal_example()
  data <- data[!(data$gene_id == "UNMAPPED" & data$replicate == "rep2"), ]
  out <- repliapa_orthogonal_validation(
    effects, data, make_test_freeze(), "MPRA", targets = "T1"
  )
  expect_true(all(c("orthogonal_delta_rep1", "orthogonal_delta_rep2") %in%
                    names(out$results)))
  expect_true(is.na(out$results$orthogonal_delta_rep2[
    out$results$gene == "UNMAPPED"
  ]))
  expect_identical(out$results$mapping_status[out$results$gene == "UNMAPPED"],
                   "ORTHOGONAL_REPLICATES_INSUFFICIENT")
})

test_that("requested targets absent from the orthogonal assay remain explicit", {
  effects <- data.frame(target = c("T1", "ABSENT"), gene = "GENE_A",
                        proximal_direction = 0.2, q_value = 0.01, tested = TRUE)
  out <- repliapa_orthogonal_validation(
    effects, make_orthogonal_example(), make_test_freeze(), "MPRA",
    targets = c("T1", "ABSENT")
  )
  expect_identical(out$accounting$target_status,
    c("SUPPORTED_WITH_MAPPED_CONSTRUCTS", "ORTHOGONAL_TARGET_NOT_PRESENT"))
  expect_equal(out$accounting$construct_n, c(2, 0))
})
