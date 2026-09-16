#' Collapse all PAS into a direction-aware proximal/distal baseline
#'
#' Within each multi-PAS gene, sites are ordered in transcript direction. The
#' first `ceiling(K / 2)` sites form the proximal bin and the remaining sites
#' form the distal bin. Gene totals are conserved exactly.
#'
#' @param counts PAS-by-sample pseudobulk count matrix.
#' @param pas_data PAS metadata aligned to rows of `counts`.
#' @param gene,position,strand Column names in `pas_data`.
#' @return A `repliapa_pd_data` object containing collapsed counts, two-bin
#'   feature metadata, site membership, and gene accounting.
#' @export
repliapa_collapse_proximal_distal <- function(counts, pas_data, gene = "gene",
                                               position = "position",
                                               strand = "strand") {
  repliapa_validate_counts(counts)
  .aligned_data(pas_data, rownames(counts), "pas_data")
  .require_columns(pas_data, c(gene, position, strand), "pas_data")
  genes <- as.character(pas_data[[gene]])
  positions <- pas_data[[position]]
  strands <- as.character(pas_data[[strand]])
  if (anyNA(genes) || any(genes == "")) {
    stop("The PAS-to-gene map is incomplete.", call. = FALSE)
  }
  if (!is.numeric(positions) || any(!is.finite(positions))) {
    stop("PAS positions must be finite numbers.", call. = FALSE)
  }
  indices <- split(seq_along(genes), genes)
  accounting <- do.call(rbind, lapply(names(indices), function(gene_value) {
    idx <- indices[[gene_value]]
    gene_strand <- unique(strands[idx])
    reason <- if (length(idx) < 2L) "SINGLE_PAS_GENE" else if (
      length(gene_strand) != 1L || !gene_strand %in% c("+", "-")
    ) "INVALID_GENE_STRAND" else "NONE"
    data.frame(gene = gene_value, pas_n = length(idx), eligible = reason == "NONE",
               reason = reason, stringsAsFactors = FALSE)
  }))
  eligible <- accounting$gene[accounting$eligible]
  if (!length(eligible)) stop("No eligible multi-PAS genes are available.", call. = FALSE)

  collapsed <- vector("list", length(eligible))
  membership <- vector("list", length(eligible))
  for (i in seq_along(eligible)) {
    gene_value <- eligible[[i]]
    idx <- indices[[gene_value]]
    gene_strand <- unique(strands[idx])
    ordered <- idx[order(
      if (gene_strand == "-") -positions[idx] else positions[idx],
      rownames(pas_data)[idx]
    )]
    proximal_n <- ceiling(length(ordered) / 2)
    bins <- rep(c("proximal", "distal"), c(proximal_n, length(ordered) - proximal_n))
    membership[[i]] <- data.frame(
      gene = gene_value, pas_id = rownames(pas_data)[ordered],
      transcript_order = seq_along(ordered), bin = bins,
      stringsAsFactors = FALSE
    )
    collapsed[[i]] <- rbind(
      Matrix::colSums(counts[ordered[bins == "proximal"], , drop = FALSE]),
      Matrix::colSums(counts[ordered[bins == "distal"], , drop = FALSE])
    )
  }
  collapsed_counts <- do.call(rbind, collapsed)
  feature_ids <- as.vector(rbind(paste0(eligible, "__proximal_bin"),
                                 paste0(eligible, "__distal_bin")))
  rownames(collapsed_counts) <- feature_ids
  colnames(collapsed_counts) <- colnames(counts)
  collapsed_pas_data <- data.frame(
    gene = rep(eligible, each = 2L),
    position = rep(c(0, 1), times = length(eligible)),
    strand = "+", bin = rep(c("proximal", "distal"), times = length(eligible)),
    row.names = feature_ids, stringsAsFactors = FALSE
  )
  if (!isTRUE(all.equal(as.numeric(Matrix::colSums(collapsed_counts)),
                        as.numeric(Matrix::colSums(counts[genes %in% eligible, , drop = FALSE]))))) {
    stop("Proximal/distal collapse did not conserve eligible-gene counts.", call. = FALSE)
  }
  output <- list(
    counts = collapsed_counts,
    pas_data = collapsed_pas_data,
    membership = do.call(rbind, membership),
    accounting = accounting,
    settings = list(gene = gene, position = position, strand = strand,
                    split_rule = "transcript_order_ceiling_half")
  )
  class(output) <- "repliapa_pd_data"
  output
}

#' Fit the direction-aware proximal/distal baseline
#'
#' This uses the same guide-identity exact test as [repliapa_fit()], changing
#' only the feature representation from all PAS to two bins. It therefore
#' isolates the information added by retaining the full PAS composition.
#'
#' @inheritParams repliapa_fit
#' @param position,strand Column names in `pas_data`.
#' @return A `repliapa_pd_fit`, inheriting from `repliapa_fit`.
#' @export
repliapa_fit_proximal_distal <- function(counts, sample_data, pas_data, target,
                                         gene = "gene", position = "position",
                                         strand = "strand", replicate = "replicate",
                                         guide = "guide", target_col = "target",
                                         ntc = "ntc", replicates = NULL,
                                         pseudocount = 0.5, min_gene_total = 20,
                                         max_combinations = 1000000) {
  pd <- repliapa_collapse_proximal_distal(counts, pas_data, gene, position, strand)
  fit <- repliapa_fit(
    pd$counts, sample_data, pd$pas_data, target = target,
    gene = "gene", replicate = replicate, guide = guide,
    target_col = target_col, ntc = ntc, replicates = replicates,
    pseudocount = pseudocount, min_gene_total = min_gene_total,
    max_combinations = max_combinations
  )
  fit$proximal_distal <- pd
  fit$settings$site_representation <- "PROXIMAL_DISTAL_TWO_BIN"
  class(fit) <- c("repliapa_pd_fit", "repliapa_fit")
  fit
}

#' Compare all-PAS and proximal/distal fits on the same target
#'
#' @param all_pas A standard `repliapa_fit` using all PAS.
#' @param proximal_distal A `repliapa_pd_fit` using the same target and samples.
#' @param fdr BH threshold used for the comparison labels.
#' @return A gene-level matched comparison with tested-set accounting.
#' @export
repliapa_compare_allpas_pd <- function(all_pas, proximal_distal, fdr = 0.05) {
  if (!inherits(all_pas, "repliapa_fit") || !inherits(proximal_distal, "repliapa_pd_fit")) {
    stop("Inputs must be all-PAS and proximal/distal RepliAPA fits.", call. = FALSE)
  }
  if (!identical(all_pas$target, proximal_distal$target) ||
      !identical(all_pas$replicates, proximal_distal$replicates) ||
      !identical(all_pas$target_guides, proximal_distal$target_guides) ||
      !identical(all_pas$ntc_guides, proximal_distal$ntc_guides)) {
    stop("Fits must use the same target and biological replicates, with identical guide units.",
         call. = FALSE)
  }
  fdr <- .validate_scalar_number(
    fdr, "fdr", lower = 0, lower_inclusive = FALSE,
    upper = 1, upper_inclusive = FALSE
  )
  all_table <- all_pas$gene_table[c("gene", "pas_n", "tested", "filtered", "abstained",
                                     "statistic", "p_value", "q_value", "effect_norm")]
  names(all_table)[-1] <- paste0("allPAS_", names(all_table)[-1])
  pd_table <- proximal_distal$gene_table[c("gene", "tested", "filtered", "abstained",
                                           "statistic", "p_value", "q_value", "effect_norm")]
  names(pd_table)[-1] <- paste0("PD_", names(pd_table)[-1])
  out <- merge(all_table, pd_table, by = "gene", all = TRUE)
  out$matched_tested <- out$allPAS_tested & out$PD_tested
  out$allPAS_significant <- out$allPAS_tested & out$allPAS_q_value <= fdr
  out$PD_significant <- out$PD_tested & out$PD_q_value <= fdr
  out$comparison_class <- ifelse(
    out$allPAS_significant & out$PD_significant, "BOTH",
    ifelse(out$allPAS_significant, "ALLPAS_ONLY",
           ifelse(out$PD_significant, "PD_ONLY", "NEITHER"))
  )
  attr(out, "accounting") <- data.frame(
    target = all_pas$target,
    universe_gene_n = nrow(out),
    matched_tested_gene_n = sum(out$matched_tested, na.rm = TRUE),
    allPAS_tested_gene_n = sum(out$allPAS_tested, na.rm = TRUE),
    PD_tested_gene_n = sum(out$PD_tested, na.rm = TRUE),
    allPAS_only_discovery_n = sum(out$comparison_class == "ALLPAS_ONLY", na.rm = TRUE),
    PD_only_discovery_n = sum(out$comparison_class == "PD_ONLY", na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  out
}
