# Figure 02 - hypoxia score by MRI tumour enhancement, one figure per signature.
#
# Reads : data/mri/scores_merged_mri_enhancement.csv
#         data/mri/pairwise_wilcox_by_enhancement.csv
# Writes: figures/fig_{tardon,buffa,ragnum}_by_tumor_enhancement.{pdf,png}

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

IN_DIR <- file.path(DATA_DIR, "mri")
merged <- utils::read.csv(file.path(IN_DIR, "scores_merged_mri_enhancement.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)

pw_path <- file.path(IN_DIR, "pairwise_wilcox_by_enhancement.csv")
pw_tab <- if (file.exists(pw_path)) {
  utils::read.csv(pw_path, stringsAsFactors = FALSE)
} else {
  data.frame()
}

out_prefix_map <- c(
  tardon_iqm = "fig_tardon_by_tumor_enhancement",
  buffa_mean = "fig_buffa_by_tumor_enhancement",
  ragnum_iqm = "fig_ragnum_by_tumor_enhancement"
)

for (sc in names(SCORE_COLS)) {
  pw_sc <- if (nrow(pw_tab)) pw_tab[pw_tab$score == sc, , drop = FALSE] else NULL
  plot_by_enhancement(
    merged,
    value_col   = sc,
    value_label = "Hypoxia score",
    out_prefix  = file.path(FIG_DIR, out_prefix_map[[sc]]),
    pw_sc       = pw_sc,
    title       = SCORE_COLS[[sc]],
    dpi         = 150
  )
}

message("Figure 02 done. Figures in: ", FIG_DIR)
