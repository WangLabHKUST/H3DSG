# Figure 03 - immune-rich cell ratio against hypoxia score, Spearman.
#
# Reads : data/immune_hypoxia/merged_immune_hypoxia_per_sample.csv
#         data/immune_hypoxia/immune_hypoxia_spearman.csv
# Writes: figures/scatter_<score>_vs_<immune>.{pdf,png}
#         figures/scatter_facet_all.{pdf,png}

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

need_pkg("ggplot2")

POINT_SIZE  <- 1.8
POINT_ALPHA <- 0.85
SMOOTH_COL  <- "grey30"

IN_DIR <- file.path(DATA_DIR, "immune_hypoxia")
merged <- utils::read.csv(file.path(IN_DIR, "merged_immune_hypoxia_per_sample.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)
corr_tab <- utils::read.csv(file.path(IN_DIR, "immune_hypoxia_spearman.csv"),
                            stringsAsFactors = FALSE)

score_cols    <- unique(corr_tab$score)
immune_cols   <- unique(corr_tab$immune)
score_labels  <- stats::setNames(corr_tab$score_label[match(score_cols, corr_tab$score)],
                                 score_cols)
immune_labels <- stats::setNames(corr_tab$immune_label[match(immune_cols, corr_tab$immune)],
                                 immune_cols)

merged[[GROUP_COL]] <- factor(merged[[GROUP_COL]])
group_levels <- levels(merged[[GROUP_COL]])
point_color_map <- stats::setNames(GROUP_POINT_COLORS[seq_along(group_levels)], group_levels)

corr_row <- function(score, immune) {
  hit <- corr_tab[corr_tab$score == score & corr_tab$immune == immune, , drop = FALSE]
  if (!nrow(hit)) NULL else hit[1, ]
}

annotation_text <- function(r, multiline_n = TRUE) {
  if (is.null(r) || is.na(r$spearman_rho)) {
    return(sprintf("n=%s", if (is.null(r)) NA else r$n))
  }
  if (multiline_n) {
    sprintf("Spearman rho = %.2f\n%s\nn = %s", r$spearman_rho, format_pval(r$p_value), r$n)
  } else {
    sprintf("rho=%.2f\n%s", r$spearman_rho, format_pval(r$p_value))
  }
}

make_scatter <- function(df, x_col, y_col, x_label, y_label, ann_lab) {
  x_pos <- max(df[[x_col]][is.finite(df[[x_col]])], na.rm = TRUE)
  y_pos <- max(df[[y_col]][is.finite(df[[y_col]])], na.rm = TRUE)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                         color = SMOOTH_COL, fill = "grey80",
                         linewidth = 0.5, alpha = 0.4, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[GROUP_COL]]),
                        size = POINT_SIZE, alpha = POINT_ALPHA, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = point_color_map, name = GROUP_COL) +
    ggplot2::annotate("text", x = x_pos, y = y_pos, label = ann_lab,
                      hjust = 1, vjust = 1, size = 3, lineheight = 1.0) +
    ggplot2::labs(title = "", x = x_label, y = y_label) +
    base_theme(angle_x = 0, legend = "right")
}

for (sc in score_cols) {
  for (ic in immune_cols) {
    p <- make_scatter(merged, x_col = ic, y_col = sc,
                      x_label = immune_labels[[ic]], y_label = "Hypoxia score",
                      ann_lab = annotation_text(corr_row(sc, ic)))
    fname <- sprintf("scatter_%s_vs_%s", sc, ic)
    ggplot2::ggsave(file.path(FIG_DIR, paste0(fname, ".pdf")), p, width = 4, height = 3)
    ggplot2::ggsave(file.path(FIG_DIR, paste0(fname, ".png")), p, width = 4, height = 3, dpi = 300)
  }
}

## ------------ facet grid: score (rows) x immune (cols) ------------
long <- do.call(rbind, lapply(score_cols, function(sc) {
  do.call(rbind, lapply(immune_cols, function(ic) {
    data.frame(
      group        = merged[[GROUP_COL]],
      score_label  = unname(score_labels[[sc]]),
      immune_label = unname(immune_labels[[ic]]),
      score_val    = merged[[sc]],
      immune_val   = merged[[ic]],
      stringsAsFactors = FALSE
    )
  }))
}))
long$score_label  <- factor(long$score_label,  levels = unname(score_labels))
long$immune_label <- factor(long$immune_label, levels = unname(immune_labels))
long$group        <- factor(long$group, levels = group_levels)

ann_df <- do.call(rbind, lapply(score_cols, function(sc) {
  do.call(rbind, lapply(immune_cols, function(ic) {
    sub_x <- merged[[ic]][is.finite(merged[[ic]])]
    sub_y <- merged[[sc]][is.finite(merged[[sc]])]
    data.frame(
      score_label  = unname(score_labels[[sc]]),
      immune_label = unname(immune_labels[[ic]]),
      ann = annotation_text(corr_row(sc, ic), multiline_n = FALSE),
      x = if (!length(sub_x)) NA_real_ else max(sub_x),
      y = if (!length(sub_y)) NA_real_ else max(sub_y),
      stringsAsFactors = FALSE
    )
  }))
}))
ann_df$score_label  <- factor(ann_df$score_label,  levels = unname(score_labels))
ann_df$immune_label <- factor(ann_df$immune_label, levels = unname(immune_labels))

p_facet <- ggplot2::ggplot(long, ggplot2::aes(x = immune_val, y = score_val)) +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x,
                       color = SMOOTH_COL, fill = "grey80",
                       linewidth = 0.5, alpha = 0.4, na.rm = TRUE) +
  ggplot2::geom_point(ggplot2::aes(color = group),
                      size = POINT_SIZE, alpha = POINT_ALPHA, na.rm = TRUE) +
  ggplot2::scale_color_manual(values = point_color_map, name = GROUP_COL) +
  ggplot2::geom_text(data = ann_df, ggplot2::aes(x = x, y = y, label = ann),
                     hjust = 1, vjust = 1, size = 3, lineheight = 1.0,
                     inherit.aes = FALSE) +
  ggplot2::facet_grid(score_label ~ immune_label, scales = "free") +
  ggplot2::labs(x = NULL, y = NULL, title = "Immune ratio vs hypoxia score (Spearman)") +
  base_theme(angle_x = 0, legend = "right")

facet_w <- 2.8 * length(immune_cols) + 1.2
facet_h <- 2.6 * length(score_cols) + 0.8
ggplot2::ggsave(file.path(FIG_DIR, "scatter_facet_all.pdf"), p_facet,
                width = facet_w, height = facet_h)
ggplot2::ggsave(file.path(FIG_DIR, "scatter_facet_all.png"), p_facet,
                width = facet_w, height = facet_h, dpi = 300)

message("Figure 03 done. Figures in: ", FIG_DIR)
