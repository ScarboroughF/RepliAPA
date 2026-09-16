.repliapa_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "RepliAPA plotting requires the suggested package `ggplot2`. ",
      "Install it with install.packages(\"ggplot2\").",
      call. = FALSE
    )
  }
}

.repliapa_plot_gene <- function(fit, gene) {
  if (!inherits(fit, "repliapa_fit")) {
    stop("`fit` must be a repliapa_fit object.", call. = FALSE)
  }
  if (!is.character(gene) || length(gene) != 1L || is.na(gene) || gene == "") {
    stop("`gene` must be one non-missing gene identifier.", call. = FALSE)
  }
  result <- fit$pas_effects[as.character(fit$pas_effects$gene) == gene, , drop = FALSE]
  if (!nrow(result)) stop("The requested gene has no tested PAS effects.", call. = FALSE)
  result
}

.repliapa_named_colours <- c(
  navy = "#17324D", blue = "#0077B6", cyan = "#00A6A6",
  coral = "#E76F51", gold = "#E9C46A", violet = "#6C5CE7",
  slate = "#7A8793", ink = "#17212B", mist = "#E8EEF3"
)

#' RepliAPA publication palette
#'
#' A restrained, colour-vision-deficiency-aware palette used consistently by
#' the package plots. Named semantic palettes prevent biological categories
#' from changing colour between figures.
#'
#' @param type Palette family: `method`, `effect`, `context`, or `pas`.
#' @param n Number of colours for the `pas` palette.
#' @return A named character vector of hexadecimal colours.
#' @export
repliapa_palette <- function(type = c("method", "effect", "context", "pas"), n = 8L) {
  type <- match.arg(type)
  colours <- .repliapa_named_colours
  if (type == "method") {
    return(c(RepliAPA = colours[["blue"]], DRIMSeq = colours[["coral"]],
             PASTA = colours[["violet"]], `Proximal/distal` = colours[["cyan"]],
             `Cell-level naive` = colours[["gold"]]))
  }
  if (type == "effect") {
    return(c(`Increased usage` = colours[["coral"]],
             `Decreased usage` = colours[["blue"]],
             `No shift` = colours[["slate"]]))
  }
  if (type == "context") {
    return(c(
      SHARED_EFFECT = colours[["cyan"]], OPPOSITE_EFFECT = colours[["coral"]],
      CELL_LINE_SPECIFIC_EFFECT = colours[["violet"]],
      SUPPORTED_NEITHER_SIGNIFICANT = colours[["slate"]],
      UNSUPPORTED_GENE_CONTEXT = colours[["mist"]]
    ))
  }
  n <- .validate_scalar_integer(n, "n", minimum = 1L)
  grDevices::hcl.colors(n, palette = "Dark 3")
}

#' RepliAPA publication theme
#'
#' @param base_size Base font size in points.
#' @param base_family Font family. The empty string uses the graphics-device
#'   default and remains portable across operating systems.
#' @param grid Display subtle horizontal (`y`), vertical (`x`), both (`xy`),
#'   or no (`none`) major grid lines.
#' @return A complete ggplot2 theme.
#' @export
repliapa_theme <- function(base_size = 11, base_family = "",
                           grid = c("y", "x", "xy", "none")) {
  .repliapa_require_ggplot2()
  grid <- match.arg(grid)
  base_size <- .validate_scalar_number(base_size, "base_size", lower = 6)
  value <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = .repliapa_named_colours[["ink"]]),
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.34,
                                         margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(colour = "#52606D", size = base_size,
                                            margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(colour = "#6B7785", size = base_size * 0.78,
                                           hjust = 0, margin = ggplot2::margin(t = 8)),
      axis.title = ggplot2::element_text(face = "bold", size = base_size * 0.93),
      axis.text = ggplot2::element_text(colour = "#334155", size = base_size * 0.84),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.key = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E8EDF2", linewidth = 0.35),
      panel.border = ggplot2::element_rect(colour = "#CBD5DF", fill = NA, linewidth = 0.45),
      strip.text = ggplot2::element_text(face = "bold", colour = .repliapa_named_colours[["navy"]]),
      strip.background = ggplot2::element_rect(fill = "#EEF4F7", colour = NA),
      plot.margin = ggplot2::margin(10, 14, 10, 10)
    )
  if (grid == "y") value <- value + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (grid == "x") value <- value + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
  if (grid == "none") value <- value + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  value
}

#' Plot a target-level all-PAS effect signature
#'
#' The target estimate is shown as a lollipop signature; independently fitted
#' biological-replicate estimates remain visible as smaller points. The plot
#' represents compositional change (delta PAS usage), not gene expression.
#'
#' @param fit A `repliapa_fit` object.
#' @param gene Gene identifier to display.
#' @param pas_data Optional PAS metadata with row names matching PAS IDs and
#'   optional `position` and `strand` columns.
#' @param show_replicates Show the replicate-specific effect estimates retained
#'   in the fit.
#' @return A ggplot object.
#' @export
repliapa_plot_pas_effect <- function(fit, gene, pas_data = NULL,
                                     show_replicates = TRUE) {
  .repliapa_require_ggplot2()
  x <- .repliapa_plot_gene(fit, gene)
  if (!is.logical(show_replicates) || length(show_replicates) != 1L || is.na(show_replicates)) {
    stop("`show_replicates` must be TRUE or FALSE.", call. = FALSE)
  }
  x$site_order <- seq_len(nrow(x))
  if (!is.null(pas_data)) {
    if (!is.data.frame(pas_data) || is.null(rownames(pas_data))) {
      stop("`pas_data` must be a data frame with PAS IDs as row names.", call. = FALSE)
    }
    index <- match(x$pas_id, rownames(pas_data))
    if (anyNA(index)) stop("`pas_data` does not contain every plotted PAS ID.", call. = FALSE)
    if ("position" %in% names(pas_data)) {
      positions <- pas_data$position[index]
      if (!is.numeric(positions) || any(!is.finite(positions))) {
        stop("PAS positions must be finite numbers.", call. = FALSE)
      }
      direction <- 1
      if ("strand" %in% names(pas_data)) {
        strands <- unique(as.character(pas_data$strand[index]))
        if (length(strands) == 1L && strands == "-") direction <- -1
      }
      x$site_order <- rank(direction * positions, ties.method = "first")
      x$position <- positions
    }
  }
  x <- x[order(x$site_order, x$pas_id), , drop = FALSE]
  x$site_index <- seq_len(nrow(x))
  x$site_label <- paste0("PAS ", x$site_index, "  \u00b7  ", x$pas_id)
  x$direction <- factor(ifelse(x$effect > 0, "Increased usage",
                               ifelse(x$effect < 0, "Decreased usage", "No shift")),
                        levels = names(repliapa_palette("effect")))
  gene_row <- fit$gene_table[as.character(fit$gene_table$gene) == gene, , drop = FALSE]
  subtitle <- paste0(
    fit$target, " vs NTC  \u00b7  ", length(fit$target_guides), " guides  \u00b7  ",
    length(fit$replicates), " biological replicate",
    if (length(fit$replicates) == 1L) "" else "s",
    if (nrow(gene_row) && is.finite(gene_row$q_value[[1]]))
      paste0("  \u00b7  BH q = ", format(gene_row$q_value[[1]], digits = 2)) else ""
  )
  plot <- ggplot2::ggplot(x, ggplot2::aes(x = effect, y = stats::reorder(site_label, site_index))) +
    ggplot2::geom_vline(xintercept = 0, colour = "#8A98A6", linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = effect, yend = stats::reorder(site_label, site_index),
                                      colour = direction), linewidth = 2.6,
                          lineend = "round", alpha = 0.82) +
    ggplot2::geom_point(ggplot2::aes(colour = direction), size = 4.1) +
    ggplot2::scale_colour_manual(values = repliapa_palette("effect"), drop = TRUE) +
    ggplot2::labs(
      title = paste0(gene, " all-PAS composition signature"), subtitle = subtitle,
      x = expression(Delta * " PAS usage"), y = NULL, colour = NULL,
      caption = "Cells contribute measurement precision; guide \u00d7 biological replicate remains the evidence unit."
    ) +
    repliapa_theme(grid = "x")

  effect_columns <- grep("^effect_", names(x), value = TRUE)
  if (show_replicates && length(effect_columns)) {
    replicate_data <- do.call(rbind, lapply(effect_columns, function(column) {
      data.frame(
        site_label = x$site_label, site_index = x$site_index,
        replicate = sub("^effect_", "", column), effect = x[[column]],
        stringsAsFactors = FALSE
      )
    }))
    plot <- plot + ggplot2::geom_point(
      data = replicate_data,
      ggplot2::aes(x = effect, y = stats::reorder(site_label, site_index), shape = replicate),
      inherit.aes = FALSE, size = 2.15, stroke = 0.7, fill = "white",
      colour = .repliapa_named_colours[["ink"]]
    ) + ggplot2::guides(shape = ggplot2::guide_legend(title = "Biological replicate"))
  }
  plot
}

#' Plot retained guide-by-replicate PAS evidence
#'
#' @param fit A `repliapa_fit` object.
#' @param gene Gene identifier to display.
#' @param limits Optional symmetric fill-scale limits. By default a robust limit
#'   is estimated from the displayed guide effects.
#' @return A ggplot object.
#' @export
repliapa_plot_guide_evidence <- function(fit, gene, limits = NULL) {
  .repliapa_require_ggplot2()
  .repliapa_plot_gene(fit, gene)
  x <- fit$guide_effects[as.character(fit$guide_effects$gene) == gene, , drop = FALSE]
  if (!nrow(x)) stop("The requested gene has no retained guide effects.", call. = FALSE)
  if (is.null(limits)) {
    bound <- as.numeric(stats::quantile(abs(x$effect), 0.98, na.rm = TRUE))
    if (!is.finite(bound) || bound == 0) bound <- max(abs(x$effect), 0.01)
    limits <- c(-bound, bound)
  }
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[[1]] >= 0 || limits[[2]] <= 0) {
    stop("`limits` must contain finite negative and positive bounds.", call. = FALSE)
  }
  x$guide <- factor(x$guide, levels = rev(unique(x$guide)))
  x$pas_id <- factor(x$pas_id, levels = unique(x$pas_id))
  ggplot2::ggplot(x, ggplot2::aes(x = pas_id, y = guide, fill = effect)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.75) +
    ggplot2::facet_grid(rows = ggplot2::vars(replicate), scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(
      low = .repliapa_named_colours[["blue"]], mid = "#F7F9FB",
      high = .repliapa_named_colours[["coral"]], midpoint = 0,
      limits = limits,
      oob = function(value, range) pmin(pmax(value, range[[1]]), range[[2]])
    ) +
    ggplot2::labs(
      title = paste0(gene, " guide-level evidence"),
      subtitle = paste0(fit$target, " guides relative to replicate-matched NTC composition"),
      x = "Poly(A) site", y = "Independent guide", fill = expression(Delta * " usage"),
      caption = "Each row is one guide; cells are not plotted as independent observations."
    ) +
    repliapa_theme(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Plot guide-unit NTC calibration
#'
#' @param calibration A `repliapa_null_calibration` object.
#' @param max_alpha Largest nominal threshold displayed.
#' @return A ggplot object.
#' @export
repliapa_plot_calibration <- function(calibration, max_alpha = 0.1) {
  .repliapa_require_ggplot2()
  if (!inherits(calibration, "repliapa_null_calibration")) {
    stop("`calibration` must be a repliapa_null_calibration object.", call. = FALSE)
  }
  max_alpha <- .validate_scalar_number(max_alpha, "max_alpha", lower = 0,
                                       lower_inclusive = FALSE, upper = 1,
                                       upper_inclusive = FALSE)
  p <- calibration$results$p_value[calibration$results$tested]
  p <- p[is.finite(p)]
  if (!length(p)) stop("The calibration object contains no finite tested P values.", call. = FALSE)
  alpha <- sort(unique(c(seq(0, max_alpha, length.out = 101), 0.01, 0.025, 0.05)))
  alpha <- alpha[alpha <= max_alpha]
  empirical <- vapply(alpha, function(value) mean(p <= value), numeric(1))
  se <- sqrt(pmax(empirical * (1 - empirical) / length(p), 0))
  curve <- data.frame(
    nominal = alpha, empirical = empirical,
    low = pmax(0, empirical - 1.96 * se), high = pmin(1, empirical + 1.96 * se)
  )
  at_five <- curve[which.min(abs(curve$nominal - 0.05)), , drop = FALSE]
  ggplot2::ggplot(curve, ggplot2::aes(x = nominal, y = empirical)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = low, ymax = high),
                         fill = .repliapa_named_colours[["blue"]], alpha = 0.13) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#7A8793",
                         linetype = "22", linewidth = 0.7) +
    ggplot2::geom_line(colour = .repliapa_named_colours[["blue"]], linewidth = 1.2) +
    ggplot2::geom_point(data = at_five, colour = .repliapa_named_colours[["coral"]],
                        size = 3.4) +
    ggplot2::annotate(
      "label", x = at_five$nominal, y = at_five$empirical,
      label = paste0("FPR@0.05 = ", format(at_five$empirical, digits = 3)),
      hjust = -0.08, vjust = 1.45, size = 3.3, linewidth = 0,
      fill = grDevices::adjustcolor("white", alpha.f = 0.88)
    ) +
    ggplot2::coord_equal(xlim = c(0, max_alpha), ylim = c(0, max_alpha), expand = FALSE) +
    ggplot2::labs(
      title = "Guide-unit null calibration",
      subtitle = paste0(format(length(p), big.mark = ","),
                        " tested gene\u2013contrast P values; whole guides remain intact"),
      x = "Nominal type-I error", y = "Empirical rejection rate",
      caption = "Ribbon: pointwise normal approximation for the empirical rejection fraction."
    ) +
    repliapa_theme(grid = "xy")
}

#' Plot frozen biological-replicate reproducibility
#'
#' @param reproducibility Output from [repliapa_reproducibility()].
#' @param label_n Number of most significant development genes to label.
#' @return A ggplot object.
#' @export
repliapa_plot_reproducibility <- function(reproducibility, label_n = 6L) {
  .repliapa_require_ggplot2()
  if (!is.data.frame(reproducibility)) {
    stop("`reproducibility` must be the data frame returned by repliapa_reproducibility().",
         call. = FALSE)
  }
  .require_columns(reproducibility,
                   c("gene", "p_value_development", "p_value_validation",
                     "effect_vector_cosine", "development_hit", "validation_hit"),
                   "reproducibility")
  label_n <- .validate_scalar_integer(label_n, "label_n", minimum = 0L)
  x <- reproducibility[is.finite(reproducibility$p_value_development) &
                         is.finite(reproducibility$p_value_validation), , drop = FALSE]
  if (!nrow(x)) stop("No genes have finite P values in both replicates.", call. = FALSE)
  floor_value <- 1 / (10 * nrow(x))
  x$development_score <- -log10(pmax(x$p_value_development, floor_value))
  x$validation_score <- -log10(pmax(x$p_value_validation, floor_value))
  x$evidence <- factor(
    ifelse(x$development_hit & x$validation_hit & x$effect_vector_cosine > 0,
           "Replicated discovery",
           ifelse(x$development_hit, "Development discovery", "Tested in both")),
    levels = c("Tested in both", "Development discovery", "Replicated discovery")
  )
  colours <- c(`Tested in both` = "#A8B3BD", `Development discovery` = "#E9C46A",
               `Replicated discovery` = "#00A6A6")
  summary <- attr(reproducibility, "summary")
  subtitle <- if (is.data.frame(summary) && nrow(summary)) {
    paste0(summary$target[[1]], "  \u00b7  rank Spearman = ",
           format(summary$rank_spearman[[1]], digits = 3),
           "  \u00b7  directional replication = ",
           format(100 * summary$directional_replication_rate[[1]], digits = 3), "%")
  } else "Independent biological replicates"
  labels <- x[order(x$p_value_development), , drop = FALSE]
  labels <- labels[seq_len(min(label_n, nrow(labels))), , drop = FALSE]
  ggplot2::ggplot(x, ggplot2::aes(x = development_score, y = validation_score)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#8D99AE",
                         linetype = "22", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = evidence, size = abs(effect_vector_cosine)),
                        alpha = 0.62) +
    ggplot2::geom_text(data = labels, ggplot2::aes(label = gene), colour = "#17212B",
                       size = 3, vjust = -0.8, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE) +
    ggplot2::scale_size_continuous(range = c(1.4, 4), limits = c(0, 1)) +
    ggplot2::labs(
      title = "Frozen biological-replicate validation", subtitle = subtitle,
      x = expression(-log[10] * " development P"),
      y = expression(-log[10] * " validation P"),
      colour = NULL, size = "Effect-vector |cosine|",
      caption = "The diagonal indicates equal evidence strength, not a fitted trend."
    ) +
    repliapa_theme(grid = "xy")
}

#' Plot all-PAS added value against a proximal/distal summary
#'
#' @param comparison Output from [repliapa_compare_allpas_pd()].
#' @return A ggplot object.
#' @export
repliapa_plot_allpas <- function(comparison) {
  .repliapa_require_ggplot2()
  if (!is.data.frame(comparison)) stop("`comparison` must be a data frame.", call. = FALSE)
  .require_columns(comparison,
                   c("gene", "allPAS_p_value", "PD_p_value", "comparison_class"),
                   "comparison")
  x <- comparison[comparison$matched_tested & is.finite(comparison$allPAS_p_value) &
                    is.finite(comparison$PD_p_value), , drop = FALSE]
  if (!nrow(x)) stop("No genes were tested by both representations.", call. = FALSE)
  floor_value <- 1 / (10 * nrow(x))
  x$PD_score <- -log10(pmax(x$PD_p_value, floor_value))
  x$allPAS_score <- -log10(pmax(x$allPAS_p_value, floor_value))
  x$comparison_class <- factor(x$comparison_class,
    levels = c("NEITHER", "PD_ONLY", "BOTH", "ALLPAS_ONLY"))
  colours <- c(NEITHER = "#B8C1CA", PD_ONLY = "#E9C46A",
               BOTH = "#00A6A6", ALLPAS_ONLY = "#6C5CE7")
  ggplot2::ggplot(x, ggplot2::aes(x = PD_score, y = allPAS_score,
                                  colour = comparison_class)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#8D99AE",
                         linetype = "22", linewidth = 0.6) +
    ggplot2::geom_point(alpha = 0.58, size = 2) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE,
      labels = c(NEITHER = "Neither", PD_ONLY = "Proximal/distal only",
                 BOTH = "Both", ALLPAS_ONLY = "All-PAS only")) +
    ggplot2::labs(
      title = "Information retained by the all-PAS representation",
      subtitle = paste0(format(nrow(x), big.mark = ","), " matched tested genes"),
      x = expression(-log[10] * " proximal/distal P"),
      y = expression(-log[10] * " all-PAS P"), colour = NULL,
      caption = "Points above the diagonal have stronger evidence in the all-PAS representation."
    ) +
    repliapa_theme(grid = "xy")
}

#' Plot cross-context support and effect direction
#'
#' @param comparison A `repliapa_context_comparison` object.
#' @return A ggplot object.
#' @export
repliapa_plot_context <- function(comparison) {
  .repliapa_require_ggplot2()
  if (!inherits(comparison, "repliapa_context_comparison")) {
    stop("`comparison` must be a repliapa_context_comparison object.", call. = FALSE)
  }
  x <- comparison$results
  if (!nrow(x)) stop("The context comparison contains no supported results.", call. = FALSE)
  primary_q <- paste0(comparison$contexts[[1]], "_q_value")
  external_q <- paste0(comparison$contexts[[2]], "_q_value")
  .require_columns(x, c("target", "gene", "effect_vector_cosine", "category",
                        primary_q, external_q), "comparison$results")
  x$joint_q <- pmax(x[[primary_q]], x[[external_q]], na.rm = TRUE)
  x$joint_q[!is.finite(x$joint_q)] <- 1
  x$joint_evidence <- -log10(pmax(x$joint_q, 1e-300))
  x$category <- factor(x$category, levels = names(repliapa_palette("context")))
  ggplot2::ggplot(x, ggplot2::aes(x = effect_vector_cosine, y = joint_evidence,
                                  colour = category)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#8D99AE", linetype = "22",
                        linewidth = 0.6) +
    ggplot2::geom_point(alpha = 0.58, size = 2) +
    ggplot2::scale_colour_manual(values = repliapa_palette("context"), drop = FALSE,
      labels = c(SHARED_EFFECT = "Shared", OPPOSITE_EFFECT = "Opposite",
                 CELL_LINE_SPECIFIC_EFFECT = "Context-specific",
                 SUPPORTED_NEITHER_SIGNIFICANT = "Supported, neither significant",
                 UNSUPPORTED_GENE_CONTEXT = "Unsupported")) +
    ggplot2::labs(
      title = paste(comparison$contexts, collapse = " \u2194 "),
      subtitle = "Direction similarity and evidence jointly supported across contexts",
      x = "All-PAS effect-vector cosine", y = expression(-log[10] * " worst-context BH q"),
      colour = NULL,
      caption = "A positive cosine supports a shared direction; unsupported genes remain explicitly accounted for."
    ) +
    repliapa_theme(grid = "xy")
}

#' Plot APA and gene-expression effects as distinct estimands
#'
#' @param classification Output from [repliapa_classify_effects()].
#' @param ge_threshold Absolute log2-fold-change reference line.
#' @param label_n Number of strongest APA effects to label.
#' @return A ggplot object.
#' @export
repliapa_plot_apa_ge <- function(classification, ge_threshold = 0.25,
                                 label_n = 8L) {
  .repliapa_require_ggplot2()
  if (!is.data.frame(classification)) stop("`classification` must be a data frame.", call. = FALSE)
  .require_columns(classification, c("gene", "log2FC", "effect_norm", "effect_class"),
                   "classification")
  ge_threshold <- .validate_scalar_number(ge_threshold, "ge_threshold", lower = 0)
  label_n <- .validate_scalar_integer(label_n, "label_n", minimum = 0L)
  x <- classification[is.finite(classification$log2FC) &
                        is.finite(classification$effect_norm), , drop = FALSE]
  if (!nrow(x)) stop("No finite APA and gene-expression effects are available.", call. = FALSE)
  x$effect_class <- factor(x$effect_class,
    levels = c("NEITHER", "GE_ONLY", "APA_ONLY", "APA_AND_GE"))
  colours <- c(NEITHER = "#B8C1CA", GE_ONLY = "#E9C46A",
               APA_ONLY = "#6C5CE7", APA_AND_GE = "#00A6A6")
  labels <- x[order(-x$effect_norm), , drop = FALSE]
  labels <- labels[seq_len(min(label_n, nrow(labels))), , drop = FALSE]
  ggplot2::ggplot(x, ggplot2::aes(x = log2FC, y = effect_norm, colour = effect_class)) +
    ggplot2::geom_vline(xintercept = c(-ge_threshold, ge_threshold), colour = "#8D99AE",
                        linetype = "22", linewidth = 0.55) +
    ggplot2::geom_point(alpha = 0.58, size = 2) +
    ggplot2::geom_text(data = labels, ggplot2::aes(label = gene), size = 3,
                       vjust = -0.8, colour = "#17212B", check_overlap = TRUE,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = colours, drop = FALSE,
      labels = c(NEITHER = "Neither", GE_ONLY = "GE only", APA_ONLY = "APA only",
                 APA_AND_GE = "APA + GE")) +
    ggplot2::labs(
      title = "APA composition and gene expression are separate estimands",
      x = expression(Delta * " gene expression (log"[2] * " fold change)"),
      y = "All-PAS composition effect norm", colour = NULL,
      caption = "Horizontal position represents total expression; vertical position represents PAS redistribution."
    ) +
    repliapa_theme(grid = "xy")
}

utils::globalVariables(c(
  "effect", "site_label", "site_index", "direction", "replicate", "guide",
  "pas_id", "nominal", "empirical", "low", "high", "development_score",
  "validation_score", "evidence", "effect_vector_cosine", "gene", "PD_score",
  "allPAS_score", "comparison_class", "joint_evidence", "category", "log2FC",
  "effect_norm", "effect_class"
))
