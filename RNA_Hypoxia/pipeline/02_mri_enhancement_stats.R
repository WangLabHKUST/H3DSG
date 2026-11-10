# Step 02 - hypoxia scores against MRI tumour enhancement.
#
# Reads : data/hypoxia_scores_per_sample.csv (data_full/ when present, so that the
#         archive keeps the patient ids), raw_data/clinical_table_s1_20260514.xlsx
# Writes: data/mri/scores_merged_mri_enhancement.csv
#         data/mri/enhancement_group_counts.csv
#         data/mri/kruskal_wallis_by_enhancement.csv
#         data/mri/pairwise_wilcox_by_enhancement.csv
#         data/mri/enhancement_group_comparison.txt

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

need_pkg("readxl")

SUBDIR <- "mri"

scores_csv <- file.path(DATA_FULL_DIR, "hypoxia_scores_per_sample.csv")
if (!file.exists(scores_csv)) {
  scores_csv <- file.path(DATA_DIR, "hypoxia_scores_per_sample.csv")
}
if (!file.exists(scores_csv)) {
  stop("Run pipeline/01_hypoxia_scores.R first.")
}

scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE)
scores$SP <- trimws(as.character(scores$SP))

mri <- read_mri_enhancement()

merged <- merge(scores, mri, by = "SP", all.x = FALSE)
merged$tumor_enhancement_short <- short_enhancement_label(merged$tumor_enhancement)
merged <- merged[
  !is.na(merged$tumor_enhancement) &
    nzchar(merged$tumor_enhancement) &
    !is.na(merged$tumor_enhancement_short) &
    nzchar(merged$tumor_enhancement_short),
  ,
  drop = FALSE
]
if (!nrow(merged)) {
  stop("No samples overlap between the scores table and the MRI table.")
}
write_tiered(merged, "scores_merged_mri_enhancement.csv", SUBDIR)

tab_n <- as.data.frame(table(merged$tumor_enhancement_short), stringsAsFactors = FALSE)
names(tab_n) <- c("tumor_enhancement_short", "n")
write_public(tab_n, "enhancement_group_counts.csv", SUBDIR)

kw_tab <- do.call(rbind, lapply(names(SCORE_COLS), function(sc) {
  row <- compare_kruskal(merged[[sc]], merged$tumor_enhancement_short)
  row$score <- sc
  row
}))
rownames(kw_tab) <- NULL
write_public(kw_tab, "kruskal_wallis_by_enhancement.csv", SUBDIR)

pw_tab <- do.call(rbind, lapply(names(SCORE_COLS), function(sc) {
  pw <- pairwise_wilcox_enhancement(merged[[sc]], merged$tumor_enhancement_short)
  if (nrow(pw)) {
    pw$score <- sc
    pw$significant_bh <- pw$p_adj_bh < PAIRWISE_ALPHA
  }
  pw
}))
if (is.null(pw_tab)) {
  pw_tab <- data.frame()
}
if (nrow(pw_tab)) {
  write_public(pw_tab, "pairwise_wilcox_by_enhancement.csv", SUBDIR)
}

## ------------ human-readable summary ------------
stat_lines <- c(
  paste0("scores: data/", basename(scores_csv)),
  paste0("MRI: raw_data/", basename(CLINICAL_XLSX), " (rows 1-", N_COHORT, " only)"),
  paste0("n_samples: ", nrow(merged)),
  paste0("pairwise significance: BH-adjusted p < ", PAIRWISE_ALPHA),
  "",
  "Tumor Enhancement counts (short labels):"
)
for (i in seq_len(nrow(tab_n))) {
  stat_lines <- c(stat_lines,
                  sprintf("  %s: %s", tab_n$tumor_enhancement_short[i], tab_n$n[i]))
}
stat_lines <- c(stat_lines, "", "Kruskal-Wallis (score ~ Tumor Enhancement):")
for (i in seq_len(nrow(kw_tab))) {
  r <- kw_tab[i, ]
  if (!is.na(r$note) && nzchar(r$note)) {
    stat_lines <- c(stat_lines, sprintf("  %s: %s", r$score, r$note))
  } else {
    stat_lines <- c(stat_lines,
                    sprintf("  %s: chi-squared=%s, p=%s (n=%s, k=%s)",
                            r$score, signif(r$statistic, 4), signif(r$p_value, 4),
                            r$n, r$n_groups))
  }
}
if (nrow(pw_tab)) {
  stat_lines <- c(stat_lines, "", "Pairwise Wilcoxon (BH-adjusted):")
  for (i in seq_len(nrow(pw_tab))) {
    r <- pw_tab[i, ]
    sig_mark <- if (isTRUE(r$significant_bh)) " *" else ""
    stat_lines <- c(stat_lines,
                    sprintf("  %s | %s vs %s: p_adj=%s%s",
                            r$score, r$group1, r$group2, signif(r$p_adj_bh, 4), sig_mark))
  }
}
write_lines_public(stat_lines, "enhancement_group_comparison.txt", SUBDIR)

message("Step 02 done. Tables in: ", file.path(DATA_DIR, SUBDIR))
