# =============================================================================
# PL vs TL: share of each malignant state per patient, boxplots and Wilcoxon tests
#
# Input : ../../results/malignant_states/05_cellstates/tumor_meta_<DATE_TAG>.csv
#             metadata of the 81,304 tumour cells with cell.type (07_score_states_and_plot.R);
#             the latest tumor_meta_<date>.csv is used unless the environment variable
#             DATE_TAG names one
#         ../../sample_metadata.tsv   group and patient_id of the samples
# Output: ../../results/malignant_states/05_cellstates/
#             state_composition_patient{,_AC_MES}_boxplot{,_p}_<DATE_TAG>.{png,pdf}
#                 "_p" = p value written on the plot; "_AC_MES" = AC-like and MES-like merged
#             state_composition_patient_AC_MES_boxplot_min{,_p}_<DATE_TAG>.{png,pdf}
#                 the same plots as the two AC_MES ones above (same function, same style);
#                 the file names are kept because figures refer to them
#             state_composition_patient_proportions_<DATE_TAG>.tsv   share of every state per
#                 patient (column scheme tells the two state sets apart)
#             state_composition_patient_tests_<DATE_TAG>.tsv         n, mean, median, p, BH
#
# Patient level only
#   - the percentage of each state among the tumour cells of a sample is computed first
#     (a missing state counts as 0), then the samples of one patient (patient_id in
#     sample_metadata.tsv) are averaged: one point per patient. TL_6_1..5 belong to one patient
#     (they share one patient_id) and give one point; there are 8 PL and 7 TL patients.
#   - two state sets:
#       separate  Cycling / OPC-like / AC-like / MES-like / OC-like
#       AC_MES    Cycling / OPC-like / AC/MES-like / OC-like (AC-like + MES-like merged before
#                 the shares are computed)
#     Cycling = S + G2M (the cell.type column)
#   - included_in_prop is not used: it downsamples the cells of TL_6_1..5 and does not change
#     the per-sample shares
#   - test: two-group wilcox.test; the table also gives the BH correction within a state set
#
# Plot size: the plotting area of every facet is fixed at 103 x 146 px at 96 dpi
# (about 2.73 x 3.86 cm), the size used in the figure. The PDF has the real physical size, the
# PNG is 300 dpi. Box plots draw no outlier points (the individual patients are the points).
#
# Runtime a few seconds. R 4.5.1.
# Paths are relative to this directory: run the script from here.
# =============================================================================

SAMPLE_META_FILE <- "../../sample_metadata.tsv"
RES_DIR      <- "../../results/malignant_states"
DIR_05_STATE <- file.path(RES_DIR, "05_cellstates")

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

SAMPLE_META <- read.delim(SAMPLE_META_FILE, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
SIDS <- SAMPLE_META$sid

#' sid -> renamed_id (the name used outside the pipeline)
to_renamed <- function(sid) {
  out <- SAMPLE_META$renamed_id[match(as.character(sid), SAMPLE_META$sid)]
  ifelse(is.na(out), as.character(sid), out)
}

say <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(...))

# =============================================================================
# Plot functions
#
# plot_multi_box(): faceted boxplot (one feature per facet, x = group): per-facet Wilcoxon p value
#   and geom_blank to make room for it under free_y.
#   base_size / line_width scale fonts and line widths together for small figures (defaults
#   13 / 0.8); the jitter has a fixed seed so that the points do not move between runs.
# save_fixed_panel(): fixes the plotting area of every facet to a physical size and then saves
#   the figure (total width and height are computed).
# =============================================================================
GROUP_LEVELS <- c("PL", "TL")
GROUP_COLORS <- c(PL = "#990000", TL = "#2b8cbe")

plot_multi_box <- function(data_long, group_var = "group", facet_var = "celltype", y_var = "proportion",
                           group_levels = GROUP_LEVELS,
                           colors = unname(GROUP_COLORS[GROUP_LEVELS]),
                           show_p = FALSE,
                           dot_size = 2.1,
                           boxplot_width = 0.55,
                           ncol = 3,
                           jitter_seed = 1,
                           base_size = 13,
                           line_width = 0.8) {
  k <- base_size / 13

  p <- ggplot(data_long,
              aes(x = factor(.data[[group_var]], levels = group_levels),
                  y = .data[[y_var]],
                  fill = .data[[group_var]])) +
    geom_boxplot(outlier.shape = NA, width = boxplot_width, linewidth = line_width, alpha = 0.9, color = "black") +
    geom_point(position = position_jitter(width = 0.18, height = 0, seed = jitter_seed),
               size = dot_size, alpha = 0.95, color = "black") +
    scale_fill_manual(values = colors) +
    facet_wrap(~ .data[[facet_var]], scales = "free_y", ncol = ncol) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10 * k),
      strip.text = element_text(size = 11 * k, face = "bold"),
      strip.background = element_rect(fill = "grey92"),
      panel.spacing = unit(k, "lines"),
      axis.title.y = element_text(angle = 90, vjust = 0.5, size = 12 * k)
    ) +
    labs(title = "Comparison across features") +
    ylab(NULL)

  if (show_p) {
    pval_df <- data_long %>%
      group_by(.data[[facet_var]]) %>%
      summarise(
        p_val = tryCatch(
          suppressWarnings(wilcox.test(
            .data[[y_var]][.data[[group_var]] == group_levels[1]],
            .data[[y_var]][.data[[group_var]] == group_levels[2]]
          )$p.value),
          error = function(e) NA_real_
        ),
        y_max = max(.data[[y_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        label   = ifelse(is.na(p_val), "n.s.", paste0("p = ", formatC(p_val, digits = 2, format = "g"))),
        y_label = y_max * 1.12,
        y_blank = y_max * 1.20,
        x_mid   = 1.5
      )

    p <- p +
      geom_blank(data = pval_df, aes(x = x_mid, y = y_blank), inherit.aes = FALSE) +
      geom_text(data = pval_df, aes(x = x_mid, y = y_label, label = label),
                inherit.aes = FALSE, size = 3.5 * k, color = "black")
  }
  p
}

# The plotting area of every facet is panel_w x panel_h inches; saves png / pdf by extension and
# returns the width and height of the whole figure (inches)
save_fixed_panel <- function(p, files, panel_w, panel_h, dpi = 300) {
  g <- ggplotGrob(p)
  pn <- g$layout[grepl("^panel", g$layout$name), ]
  g$widths[unique(pn$l)]  <- unit(panel_w, "in")
  g$heights[unique(pn$t)] <- unit(panel_h, "in")
  grDevices::pdf(NULL)                                   # a device is needed to convert units such as grobwidth
  w <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  h <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  grDevices::dev.off()
  for (f in files) ggsave(f, g, width = w, height = h, units = "in", dpi = dpi)
  invisible(c(width = w, height = h))
}

# =============================================================================
# Input
# =============================================================================
DATE_TAG <- Sys.getenv("DATE_TAG")
if (!nzchar(DATE_TAG)) {
  f <- sort(list.files(DIR_05_STATE, pattern = "^tumor_meta_[0-9]{8}\\.csv$"), decreasing = TRUE)
  if (!length(f)) stop("No tumor_meta_<date>.csv found, run 07_score_states_and_plot.R first")
  DATE_TAG <- sub("^tumor_meta_([0-9]{8})\\.csv$", "\\1", f[1])
}
meta_file <- file.path(DIR_05_STATE, sprintf("tumor_meta_%s.csv", DATE_TAG))
say("input:", meta_file)

SCHEMES <- list(
  separate = list(order = c("Cycling", "OPC-like", "AC-like", "MES-like", "OC-like"),
                  map = identity, file_tag = ""),
  AC_MES   = list(order = c("Cycling", "OPC-like", "AC/MES-like", "OC-like"),
                  map = function(x) ifelse(x %in% c("AC-like", "MES-like"), "AC/MES-like", x), file_tag = "_AC_MES")
)

meta <- read.csv(meta_file, row.names = 1, check.names = FALSE)
stopifnot(setequal(unique(meta$orig.ident), SIDS), !anyNA(meta$cell.type),
          all(meta$cell.type %in% SCHEMES$separate$order))
# The group comes from sample_metadata.tsv; check it against the Group column of the metadata
grp_meta <- unique(data.frame(sid = meta$orig.ident, Group = meta$Group))
stopifnot(all(grp_meta$Group == SAMPLE_META$group[match(grp_meta$sid, SAMPLE_META$sid)]))
say("tumor cells:", nrow(meta), "; samples:", length(unique(meta$orig.ident)),
    "; patients:", length(unique(SAMPLE_META$patient_id)))

## ---- patient-level shares: percentage within the sample (missing = 0) -> mean over the samples of a patient
patient_prop <- function(scheme) {
  s <- SCHEMES[[scheme]]
  data.frame(sid = meta$orig.ident, state = factor(s$map(meta$cell.type), levels = s$order)) %>%
    count(sid, state, name = "n_cells", .drop = FALSE) %>%
    group_by(sid) %>% mutate(proportion = n_cells / sum(n_cells) * 100) %>% ungroup() %>%
    left_join(SAMPLE_META[, c("sid", "group", "patient_id")], by = "sid") %>%
    group_by(patient_id, group, state) %>%
    summarise(proportion = mean(proportion), n_samples = n(),
              samples = paste(to_renamed(sid), collapse = ","), .groups = "drop") %>%
    mutate(scheme = scheme)
}
props <- lapply(setNames(names(SCHEMES), names(SCHEMES)), patient_prop)
stopifnot(all(vapply(props, function(d) all(abs(tapply(d$proportion, d$patient_id, sum) - 100) < 1e-8), logical(1))))

## ---- tests ------------------------------------------------------------------------------
tests <- bind_rows(lapply(props, function(d) d %>%
  group_by(scheme, state) %>%
  summarise(n_PL = sum(group == "PL"), n_TL = sum(group == "TL"),
            mean_PL = mean(proportion[group == "PL"]), mean_TL = mean(proportion[group == "TL"]),
            median_PL = median(proportion[group == "PL"]), median_TL = median(proportion[group == "TL"]),
            wilcox_p = suppressWarnings(wilcox.test(proportion[group == "PL"], proportion[group == "TL"])$p.value),
            .groups = "drop") %>%
  mutate(p_adj_BH = p.adjust(wilcox_p, method = "BH"))))
print(as.data.frame(tests), digits = 3)

write.table(bind_rows(props) %>% select(scheme, patient_id, group, samples, n_samples, state, proportion),
            file.path(DIR_05_STATE, sprintf("state_composition_patient_proportions_%s.tsv", DATE_TAG)),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tests, file.path(DIR_05_STATE, sprintf("state_composition_patient_tests_%s.tsv", DATE_TAG)),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ---- plots (all facets in one row) ---------------------------------------------------------
## Plotting area of a facet = 103 x 146 px at 96 dpi = 1.07 x 1.52 in = 2.73 x 3.86 cm; font
## sizes, points and lines are matched to that size. The PDF has the physical size, the PNG is
## 300 dpi and is scaled to the size above when placed in a figure.
PANEL_W_IN <- 103 / 96
PANEL_H_IN <- 146 / 96
## point size / box line width / axis line width, the same for all plots
DOT_SIZE <- 0.7
BOX_LW   <- 0.4
AXIS_LW  <- 0.4
std_box <- function(df, n_facets, show_p)
  plot_multi_box(df, group_var = "group", y_var = "proportion", facet_var = "state",
                 ncol = n_facets, show_p = show_p,
                 base_size = 9, line_width = BOX_LW, dot_size = DOT_SIZE) +
    theme(axis.line = element_line(linewidth = AXIS_LW), axis.ticks = element_line(linewidth = AXIS_LW),
          strip.background = element_rect(fill = "grey92", linewidth = AXIS_LW)) +
    labs(title = "Tumor state composition (per patient, % of malignant cells)")
for (scheme in names(SCHEMES)) {
  s <- SCHEMES[[scheme]]
  df <- props[[scheme]] %>% mutate(state = factor(state, levels = s$order))
  for (show_p in c(FALSE, TRUE)) {
    stem <- sprintf("state_composition_patient%s_boxplot%s_%s", s$file_tag, if (show_p) "_p" else "", DATE_TAG)
    wh <- save_fixed_panel(std_box(df, length(s$order), show_p),
                           file.path(DIR_05_STATE, paste0(stem, c(".png", ".pdf"))), PANEL_W_IN, PANEL_H_IN)
    say("saved", stem, "{.png,.pdf}", sprintf("%.2f x %.2f in", wh[1], wh[2]))
  }
}

## ---- final-style files: only the AC_MES set -------------------------------------------------
## Same std_box() and the same sizes as the standard AC_MES plots above, so the _min{,_p} files
## have the same content as the two standard AC_MES files; the _min names are kept for the
## figures that refer to them.
df_min <- props$AC_MES %>% mutate(state = factor(state, levels = SCHEMES$AC_MES$order))
for (show_p in c(FALSE, TRUE)) {
  stem <- sprintf("state_composition_patient_AC_MES_boxplot_min%s_%s", if (show_p) "_p" else "", DATE_TAG)
  wh <- save_fixed_panel(std_box(df_min, length(SCHEMES$AC_MES$order), show_p),
                         file.path(DIR_05_STATE, paste0(stem, c(".png", ".pdf"))), PANEL_W_IN, PANEL_H_IN)
  say("saved", stem, "{.png,.pdf}", sprintf("%.2f x %.2f in", wh[1], wh[2]))
}
