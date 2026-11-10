# Step 04 - mIHC immune-rich cell ratios against MRI tumour enhancement.
#
# Reads : raw_data/ct_prop_adjusted_31_samples.csv,
#         raw_data/ID_relationship_36cases_NF1_TP53.csv,
#         raw_data/clinical_table_s1_20260514.xlsx
# Writes: data/immune_mri/merged_immune_mri_per_sample.csv
#         data/immune_mri/enhancement_group_counts.csv
#         data/immune_mri/kruskal_wallis_immune_by_enhancement.csv
#         data/immune_mri/pairwise_wilcox_immune_by_enhancement.csv
#         data/immune_mri/enhancement_immune_group_comparison.txt

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

SUBDIR <- "immune_mri"

ct <- read_ct_proportions()
present_components <- attr(ct, "present_components")

id_map <- utils::read.csv(ID_MAP_CSV, stringsAsFactors = FALSE)
miss_m <- setdiff(c("Original_ID", "ID_in_network"), names(id_map))
if (length(miss_m)) {
  stop("ID map missing columns: ", paste(miss_m, collapse = ", "))
}
id_map$Original_ID <- trimws(as.character(id_map$Original_ID))
id_map$SP <- paste0("S", trimws(as.character(id_map$ID_in_network)))

mri <- read_mri_enhancement()

merged <- merge(ct, id_map[, c("Original_ID", "SP"), drop = FALSE],
                by = "Original_ID", all.x = FALSE)
merged <- merge(merged, mri, by = "SP", all.x = FALSE)
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
  stop("No samples overlap between the cell-type table, the ID map and the MRI table.")
}

im <- immune_columns_and_labels(merged, present_components)
immune_cols   <- im$cols
immune_labels <- im$labels

write_tiered(
  merged[, c("SP", "Original_ID", "tumor_enhancement", "tumor_enhancement_short",
             immune_cols), drop = FALSE],
  "merged_immune_mri_per_sample.csv",
  SUBDIR
)

tab_n <- as.data.frame(table(merged$tumor_enhancement_short), stringsAsFactors = FALSE)
names(tab_n) <- c("tumor_enhancement_short", "n")
write_public(tab_n, "enhancement_group_counts.csv", SUBDIR)

kw_tab <- do.call(rbind, lapply(immune_cols, function(ic) {
  row <- compare_kruskal(merged[[ic]], merged$tumor_enhancement_short)
  row$immune <- ic
  row
}))
rownames(kw_tab) <- NULL
write_public(kw_tab, "kruskal_wallis_immune_by_enhancement.csv", SUBDIR)

pw_tab <- do.call(rbind, lapply(immune_cols, function(ic) {
  pw <- pairwise_wilcox_enhancement(merged[[ic]], merged$tumor_enhancement_short)
  if (nrow(pw)) {
    pw$immune <- ic
    pw$significant_bh <- pw$p_adj_bh < PAIRWISE_ALPHA
  }
  pw
}))
if (is.null(pw_tab)) {
  pw_tab <- data.frame()
}
if (nrow(pw_tab)) {
  write_public(pw_tab, "pairwise_wilcox_immune_by_enhancement.csv", SUBDIR)
}

## ------------ human-readable summary ------------
stat_lines <- c(
  paste0("cell-type proportions: raw_data/", basename(CT_PROP_CSV)),
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
stat_lines <- c(stat_lines, "", "Kruskal-Wallis (immune ratio ~ Tumor Enhancement):")
for (i in seq_len(nrow(kw_tab))) {
  r <- kw_tab[i, ]
  if (!is.na(r$note) && nzchar(r$note)) {
    stat_lines <- c(stat_lines, sprintf("  %s: %s", r$immune, r$note))
  } else {
    stat_lines <- c(stat_lines,
                    sprintf("  %s: chi-squared=%s, p=%s (n=%s, k=%s)",
                            r$immune, signif(r$statistic, 4), signif(r$p_value, 4),
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
                            r$immune, r$group1, r$group2, signif(r$p_adj_bh, 4), sig_mark))
  }
}
write_lines_public(stat_lines, "enhancement_immune_group_comparison.txt", SUBDIR)

message("Step 04 done. Tables in: ", file.path(DATA_DIR, SUBDIR))
