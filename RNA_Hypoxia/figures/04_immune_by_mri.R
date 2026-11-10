# Figure 04 - immune-rich cell ratio by MRI tumour enhancement.
#
# Reads : data/immune_mri/merged_immune_mri_per_sample.csv
#         data/immune_mri/pairwise_wilcox_immune_by_enhancement.csv
#         data/immune_mri/kruskal_wallis_immune_by_enhancement.csv (column order)
# Writes: figures/fig_{T_NK,Myeloid,immune_total}_by_tumor_enhancement.{pdf,png}

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

IN_DIR <- file.path(DATA_DIR, "immune_mri")
merged <- utils::read.csv(file.path(IN_DIR, "merged_immune_mri_per_sample.csv"),
                          stringsAsFactors = FALSE, check.names = FALSE)

pw_path <- file.path(IN_DIR, "pairwise_wilcox_immune_by_enhancement.csv")
pw_tab <- if (file.exists(pw_path)) {
  utils::read.csv(pw_path, stringsAsFactors = FALSE)
} else {
  data.frame()
}

im <- immune_columns_and_labels(merged, IMMUNE_RICH_COMPONENTS)
immune_cols   <- im$cols
immune_labels <- im$labels

for (ic in immune_cols) {
  pw_ic <- if (nrow(pw_tab)) pw_tab[pw_tab$immune == ic, , drop = FALSE] else NULL
  plot_by_enhancement(
    merged,
    value_col   = ic,
    value_label = immune_labels[[ic]],
    out_prefix  = file.path(FIG_DIR, paste0("fig_", ic, "_by_tumor_enhancement")),
    pw_sc       = pw_ic,
    title       = immune_labels[[ic]],
    dpi         = 300
  )
}

message("Figure 04 done. Figures in: ", FIG_DIR)
