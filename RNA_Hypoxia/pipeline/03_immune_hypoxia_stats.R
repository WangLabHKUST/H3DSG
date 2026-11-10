# Step 03 - Spearman correlation between mIHC immune-rich cell ratios
# (T/NK, Myeloid, immune total) and the bulk hypoxia scores.
#
# Reads : data/hypoxia_scores_per_sample.csv (data_full/ when present),
#         raw_data/ct_prop_adjusted_31_samples.csv,
#         raw_data/ID_relationship_36cases_NF1_TP53.csv
# Writes: data/immune_hypoxia/merged_immune_hypoxia_per_sample.csv
#         data/immune_hypoxia/immune_hypoxia_spearman.csv

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

SUBDIR <- "immune_hypoxia"

# Which scores to correlate. The 2026-05 version of this analysis had Buffa and
# Ragnum commented out and ran Tardon alone; all three are enabled here so that
# the table and the facet figure cover the same signatures as every other step.
SCORES_FOR_CORR <- names(SCORE_COLS)

scores_csv <- file.path(DATA_FULL_DIR, "hypoxia_scores_per_sample.csv")
if (!file.exists(scores_csv)) {
  scores_csv <- file.path(DATA_DIR, "hypoxia_scores_per_sample.csv")
}
if (!file.exists(scores_csv)) {
  stop("Run pipeline/01_hypoxia_scores.R first.")
}

scores <- utils::read.csv(scores_csv, stringsAsFactors = FALSE, check.names = FALSE)
need_score_cols <- c("PID", GROUP_COL, SCORES_FOR_CORR)
miss_s <- setdiff(need_score_cols, names(scores))
if (length(miss_s)) {
  stop("Scores table missing columns: ", paste(miss_s, collapse = ", "),
       ". The merge needs PID, so run step 01 with raw_data/ present.")
}
scores$PID <- trimws(as.character(scores$PID))

id_map <- utils::read.csv(ID_MAP_CSV, stringsAsFactors = FALSE)
miss_m <- setdiff(c("Original_ID", "PID"), names(id_map))
if (length(miss_m)) {
  stop("ID map missing columns: ", paste(miss_m, collapse = ", "))
}
id_map$Original_ID <- trimws(as.character(id_map$Original_ID))
id_map$PID         <- trimws(as.character(id_map$PID))

ct <- read_ct_proportions()
present_components <- attr(ct, "present_components")

merged <- merge(ct, id_map[, c("Original_ID", "PID"), drop = FALSE],
                by = "Original_ID", all.x = FALSE)
merged <- merge(merged,
                scores[, c("SP", "PID", GROUP_COL, SCORES_FOR_CORR), drop = FALSE],
                by = "PID", all.x = FALSE)
if (!nrow(merged)) {
  stop("No samples overlap between the cell-type table, the ID map and the scores.")
}

im <- immune_columns_and_labels(merged, present_components)
immune_cols   <- im$cols
immune_labels <- im$labels

write_tiered(
  merged[, c("SP", "PID", "Original_ID", GROUP_COL, immune_cols, SCORES_FOR_CORR),
         drop = FALSE],
  "merged_immune_hypoxia_per_sample.csv",
  SUBDIR
)

## ------------ Spearman table ------------
corr_rows <- list()
for (sc in SCORES_FOR_CORR) {
  for (ic in immune_cols) {
    res <- spearman_safe(merged[[sc]], merged[[ic]])
    corr_rows[[length(corr_rows) + 1L]] <- data.frame(
      score = sc,
      score_label = unname(SCORE_COLS[[sc]]),
      immune = ic,
      immune_label = unname(immune_labels[[ic]]),
      n = res$n,
      spearman_rho = res$rho,
      p_value = res$p_value,
      note = ifelse(is.null(res$note), NA_character_, res$note),
      stringsAsFactors = FALSE
    )
  }
}
corr_tab <- do.call(rbind, corr_rows)
corr_tab$p_adj_bh <- stats::p.adjust(corr_tab$p_value, method = "BH")
write_public(corr_tab, "immune_hypoxia_spearman.csv", SUBDIR)

message("Step 03 done. Tables in: ", file.path(DATA_DIR, SUBDIR))
