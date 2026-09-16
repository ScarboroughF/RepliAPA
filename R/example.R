#' Construct a deterministic RepliAPA example
#'
#' Returns a small synthetic count matrix with two biological replicates,
#' three guides targeting `T1`, four non-targeting guides, and three genes.
#' It is intended for interface demonstrations, not performance claims.
#'
#' @return A list with `counts`, `sample_data`, and `pas_data`.
#' @export
repliapa_example <- function() {
  replicates <- c("R1", "R2")
  target_guides <- paste0("T1_g", 1:3)
  ntc_guides <- paste0("NT_g", 1:4)
  sample_data <- do.call(rbind, lapply(replicates, function(rep_value) {
    guides <- c(target_guides, ntc_guides)
    data.frame(replicate = rep_value, guide = guides,
               target = ifelse(guides %in% target_guides, "T1", "NT"),
               ntc = guides %in% ntc_guides, stringsAsFactors = FALSE)
  }))
  rownames(sample_data) <- paste("example", sample_data$replicate, sample_data$guide, sep = "|")
  pas_data <- data.frame(
    gene = c(rep("GENE_A", 3), rep("GENE_B", 2), "GENE_C"),
    position = c(100, 200, 300, 100, 250, 100), strand = "+",
    row.names = paste0("pas", 1:6), stringsAsFactors = FALSE
  )
  counts <- matrix(0, nrow = nrow(pas_data), ncol = nrow(sample_data),
                   dimnames = list(rownames(pas_data), rownames(sample_data)))
  for (j in seq_len(ncol(counts))) {
    is_target <- !sample_data$ntc[j]
    total_a <- 90 + 7 * j
    prop_a <- if (is_target) c(.45, .25, .30) else c(.20, .30, .50)
    value_a <- round(total_a * prop_a)
    value_a[3] <- total_a - sum(value_a[1:2])
    total_b <- 65 + 3 * j
    value_b <- c(round(total_b * .4), 0)
    value_b[2] <- total_b - value_b[1]
    counts[, j] <- c(value_a, value_b, 20 + j)
  }
  list(counts = counts, sample_data = sample_data, pas_data = pas_data)
}
