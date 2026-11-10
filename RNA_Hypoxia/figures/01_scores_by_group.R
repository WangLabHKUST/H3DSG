# Figure 01 - hypoxia score by group (PL vs TL), one figure per signature.
#
# Reads : data/hypoxia_scores_per_sample.csv
#         data/group_comparison_hypoxia_scores.csv  (for the bracket p-value)
# Writes: figures/fig_{tardon_iqm,buffa_mean,ragnum_iqm}_by_group.{pdf,png}

.root <- local({
  hit <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(hit)) {
    dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/"))
  } else {
    getwd()
  }
  while (!file.exists(file.path(d, "common.R")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(.root, "common.R"))

need_pkg("ggplot2", "ggsignif")

scores <- utils::read.csv(file.path(DATA_DIR, "hypoxia_scores_per_sample.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
cmp_tab <- utils::read.csv(file.path(DATA_DIR, "group_comparison_hypoxia_scores.csv"),
                           stringsAsFactors = FALSE)

miss <- setdiff(c(GROUP_COL, names(SCORE_COLS)), names(scores))
if (length(miss)) {
  stop("Scores table missing columns: ", paste(miss, collapse = ", "))
}

df <- scores[!is.na(scores[[GROUP_COL]]) & nzchar(scores[[GROUP_COL]]), , drop = FALSE]
lvls <- sort(unique(df[[GROUP_COL]]))
if (length(lvls) != 2L) {
  stop("Expected two group levels; found: ", paste(lvls, collapse = ", "))
}
df[[GROUP_COL]] <- factor(df[[GROUP_COL]], levels = lvls)

#' Box plot across the two groups with the Wilcoxon p-value as a bracket.
#' outlier.shape = NA because the jitter layer already draws every sample; with
#' the default the same sample would appear twice, once in the box and once
#' offset by the jitter.
plot_score_by_group <- function(df, score_col, title, p_value) {
  finite_mask <- is.finite(df[[score_col]])
  n_per_level <- table(factor(df[[GROUP_COL]][finite_mask], levels = levels(df[[GROUP_COL]])))
  x_tick_labels <- sprintf("%s (n=%d)", names(n_per_level), as.integer(n_per_level))
  names(x_tick_labels) <- names(n_per_level)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[GROUP_COL]], y = .data[[score_col]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[GROUP_COL]]),
                          width = 0.55, outlier.shape = NA,
                          size = 0.5, color = "black", show.legend = FALSE) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.15, height = 0),
                         size = 1.6, alpha = 0.8, color = "black") +
    ggplot2::scale_fill_manual(values = GROUP_FILL_COLORS) +
    ggplot2::scale_x_discrete(labels = x_tick_labels) +
    ggplot2::labs(x = NULL, y = "Hypoxia Score") +
    base_theme()

  if (is.finite(p_value)) {
    y_vals <- df[[score_col]][is.finite(df[[score_col]])]
    y_rng  <- range(y_vals)
    y_top  <- y_rng[2] + 0.08 * diff(y_rng)
    p <- p +
      ggsignif::geom_signif(
        comparisons = list(levels(df[[GROUP_COL]])),
        annotations = signif_annotation(p_value),
        y_position  = y_top,
        tip_length  = 0.015,
        textsize    = 3,
        vjust       = -0.2,
        size        = 0.4
      ) +
      ggplot2::expand_limits(y = y_top + 0.12 * diff(y_rng))
  }
  p
}

pval_for <- function(score) {
  hit <- cmp_tab[cmp_tab$score == score, , drop = FALSE]
  if (!nrow(hit) || !"wilcox_p" %in% names(hit)) NA_real_ else hit$wilcox_p[1]
}

# Panel sizes are per signature: the Tardon panel is the small one used in the
# main figure, the other two are drawn at the larger supplementary size.
fig_size <- list(
  tardon_iqm = list(pdf = c(2.5, 2.0), png = c(3.0, 3.0), dpi = 300),
  buffa_mean = list(pdf = c(3.5, 3.0), png = c(3.5, 3.0), dpi = 150),
  ragnum_iqm = list(pdf = c(3.5, 3.0), png = c(3.5, 3.0), dpi = 150)
)

for (sc in names(SCORE_COLS)) {
  p <- plot_score_by_group(df, sc, SCORE_COLS[[sc]], pval_for(sc))
  sz <- fig_size[[sc]]
  ggplot2::ggsave(file.path(FIG_DIR, paste0("fig_", sc, "_by_group.pdf")),
                  p, width = sz$pdf[1], height = sz$pdf[2])
  ggplot2::ggsave(file.path(FIG_DIR, paste0("fig_", sc, "_by_group.png")),
                  p, width = sz$png[1], height = sz$png[2], dpi = sz$dpi)
}

message("Figure 01 done. Figures in: ", FIG_DIR)
