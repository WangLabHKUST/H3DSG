# Step 01 - hypoxia signature scores on bulk RNA, and the PL vs TL comparison.
#
# Reads : raw_data/TPM.csv, raw_data/ID_relationship_36cases_NF1_TP53.csv
#         data/{tardon,buffa,ragnum}_genes.csv
# Writes: data/hypoxia_scores_per_sample.csv      (one row per sample)
#         data/score_spearman_correlation.csv     (3x3, between the signatures)
#         data/group_comparison_hypoxia_scores.csv/.txt
#
# Scores: Tardon and Ragnum are summarised with the interquartile mean, Buffa
# with the mean, following Di Giovannantonio et al., Cell Genomics 2025, which
# found those score/signature combinations to perform best.

# Locate the module root (the directory holding common.R) and load shared code.
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

TARDON_GENES_CSV <- file.path(DATA_DIR, "tardon_genes.csv")
BUFFA_GENES_CSV  <- file.path(DATA_DIR, "buffa_genes.csv")
RAGNUM_GENES_CSV <- file.path(DATA_DIR, "ragnum_genes.csv")

transform_expr <- function(mat, mode) {
  switch(
    mode,
    "log2(x+1)" = log2(as.matrix(mat) + 1),
    "log1p" = log1p(as.matrix(mat)),
    "none" = as.matrix(mat),
    stop("Unknown EXPR_TRANSFORM: ", mode)
  )
}

score_samples <- function(mat, genes, score_fun) {
  present <- intersect(rownames(mat), unique(genes))
  if (!length(present)) {
    stop("No signature genes found in bulk rownames.")
  }
  sub <- mat[present, , drop = FALSE]
  scores <- apply(sub, 2L, function(v) score_fun(v))
  n_used <- apply(sub, 2L, function(v) sum(is.finite(v)))
  list(scores = scores, n_genes_used = n_used)
}

## ------------ cohort and expression matrix ------------
meta <- prepare_cohort_metadata(ID_MAP_CSV, "Group")

bulk <- utils::read.csv(BULK_TPM_CSV, check.names = FALSE, stringsAsFactors = FALSE)
rownames(bulk) <- bulk[[1]]
bulk <- bulk[, -1, drop = FALSE]
# Some delivered column headers carry typos. The corrections live in
# raw_data/column_fixes.csv rather than here, so that no patient id is written
# into code that goes to the public repository.
fixes_csv <- file.path(RAW_DIR, "column_fixes.csv")
if (file.exists(fixes_csv)) {
  fixes <- utils::read.csv(fixes_csv, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(fixes))) {
    hit <- colnames(bulk) == trimws(fixes$from[i])
    if (any(hit)) {
      colnames(bulk)[hit] <- trimws(fixes$to[i])
    }
  }
}
for (j in seq_len(ncol(bulk))) {
  bulk[[j]] <- suppressWarnings(as.numeric(bulk[[j]]))
}

keep_pid <- intersect(meta$PID, colnames(bulk))
if (!length(keep_pid)) {
  stop("No overlap between the ID map and the bulk column names.")
}
miss_bulk <- setdiff(meta$PID, colnames(bulk))
if (length(miss_bulk)) {
  warning("Dropped from cohort (in ID map, not in bulk): ",
          paste(miss_bulk, collapse = ", "), call. = FALSE)
}
meta <- meta[meta$PID %in% keep_pid, , drop = FALSE]
meta <- meta[match(keep_pid, meta$PID), , drop = FALSE]
bulk <- bulk[, meta$PID, drop = FALSE]

mat <- transform_expr(bulk, EXPR_TRANSFORM)

## ------------ scores ------------
st_t <- score_samples(mat, read_gene_list(TARDON_GENES_CSV), iqm)
st_b <- score_samples(mat, read_gene_list(BUFFA_GENES_CSV), function(v) mean(v, na.rm = TRUE))
st_r <- score_samples(mat, read_gene_list(RAGNUM_GENES_CSV), iqm)

out <- data.frame(
  SP = meta$SP,
  PID = meta$PID,
  group = meta$group,
  tardon_iqm = as.numeric(st_t$scores),
  tardon_n_genes_used = as.integer(st_t$n_genes_used),
  buffa_mean = as.numeric(st_b$scores),
  buffa_n_genes_used = as.integer(st_b$n_genes_used),
  ragnum_iqm = as.numeric(st_r$scores),
  ragnum_n_genes_used = as.integer(st_r$n_genes_used),
  expr_transform = EXPR_TRANSFORM,
  stringsAsFactors = FALSE
)
write_tiered(out, "hypoxia_scores_per_sample.csv")

cor_s <- stats::cor(
  cbind(out$tardon_iqm, out$buffa_mean, out$ragnum_iqm),
  use = "pairwise.complete.obs",
  method = "spearman"
)
colnames(cor_s) <- rownames(cor_s) <- names(SCORE_COLS)
utils::write.csv(cor_s, file.path(DATA_DIR, "score_spearman_correlation.csv"))

## ------------ PL vs TL ------------
out_cmp <- out[!is.na(out$group) & nzchar(out$group), , drop = FALSE]
tg <- sort(unique(out_cmp$group))
if (length(tg) != 2L) {
  stop("Expected two group levels; found: ", paste(tg, collapse = ", "))
}
g1 <- tg[1]
g2 <- tg[2]

cmps <- lapply(names(SCORE_COLS), function(sc) {
  compare_two_groups(out_cmp[[sc]], out_cmp$group, g1, g2)
})
names(cmps) <- names(SCORE_COLS)

cmp_tab <- do.call(rbind, lapply(names(cmps), function(sc) {
  cmp <- cmps[[sc]]
  if (!is.null(cmp$note)) {
    return(data.frame(score = sc, note = cmp$note, stringsAsFactors = FALSE))
  }
  data.frame(
    score = sc, group1 = g1, group2 = g2,
    n1 = cmp$n1, n2 = cmp$n2,
    median1 = cmp$median1, median2 = cmp$median2,
    median_diff = cmp$median_diff_g1_minus_g2,
    wilcox_W = cmp$wilcox_W, wilcox_p = cmp$wilcox_p,
    conf_int_low = cmp$conf_int_low, conf_int_high = cmp$conf_int_high,
    cliffs_delta = cmp$cliffs_delta_g1_vs_g2,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}))
write_public(cmp_tab, "group_comparison_hypoxia_scores.csv")

fmt_cmp <- function(name, cmp) {
  if (!is.null(cmp$note)) {
    return(c(name, cmp$note, ""))
  }
  c(
    name,
    sprintf("n_%s=%s n_%s=%s median_%s=%s median_%s=%s median_diff=%s",
            g1, cmp$n1, g2, cmp$n2, g1, signif(cmp$median1, 5), g2,
            signif(cmp$median2, 5), signif(cmp$median_diff_g1_minus_g2, 5)),
    sprintf("Wilcoxon p=%s (95%% CI %s, %s) W=%s",
            signif(cmp$wilcox_p, 5), signif(cmp$conf_int_low, 5),
            signif(cmp$conf_int_high, 5), signif(cmp$wilcox_W, 5)),
    sprintf("Cliff delta (%s vs %s)=%s", g1, g2, signif(cmp$cliffs_delta_g1_vs_g2, 5)),
    ""
  )
}

stat_lines <- c(
  paste0("expr_transform: ", EXPR_TRANSFORM),
  paste0("groups: ", g1, " vs ", g2),
  paste0("n_samples: ", nrow(out)),
  "",
  fmt_cmp("Tardon_IQM", cmps$tardon_iqm),
  fmt_cmp("Buffa_mean", cmps$buffa_mean),
  fmt_cmp("Ragnum_IQM", cmps$ragnum_iqm),
  "Spearman correlation:",
  utils::capture.output(print(round(cor_s, 3))),
  ""
)
write_lines_public(stat_lines, "group_comparison_hypoxia_scores.txt")

message("Step 01 done. Tables in: ", DATA_DIR)
