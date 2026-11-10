# Shared configuration and helpers for the RNA hypoxia module.
#
# Every script under pipeline/ and figures/ sources this file, so paths, cohort
# constants, statistics and plot styling are defined exactly once.
#
# Paths are resolved from the location of this file, so the module can be moved
# or renamed without editing a single script.

## ------------ R 4.5.1 workaround ------------
# `+<-` is part of base R only from 4.5.2 onward. ggplot2 4.0+ is compiled under
# 4.5.2 and its bytecode calls `+<-` at runtime, so on R 4.5.1 every `p + layer`
# fails with "could not find function '+<-'". Remove once R is upgraded.
if (!exists("+<-", mode = "function")) `+<-` <- function(e1, e2) e1 + e2

## ------------ paths ------------
# ROOT is the directory holding this file.
ROOT <- local({
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  start <- if (length(hit)) {
    dirname(normalizePath(sub("^--file=", "", hit[1]), winslash = "/", mustWork = FALSE))
  } else {
    normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  }
  d <- start
  repeat {
    if (file.exists(file.path(d, "common.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) {
      stop("common.R not found above ", start, "; run the script from inside the module.")
    }
    d <- parent
  }
})

RAW_DIR       <- file.path(ROOT, "raw_data")   # archive only, never committed
DATA_DIR      <- file.path(ROOT, "data")       # small SP-keyed tables, committed
DATA_FULL_DIR <- file.path(ROOT, "data_full")  # same tables plus patient ids, archive only
FIG_DIR       <- file.path(ROOT, "figures")

# Raw inputs. See raw_data/README.md for where each one came from.
BULK_TPM_CSV  <- file.path(RAW_DIR, "TPM.csv")
ID_MAP_CSV    <- file.path(RAW_DIR, "ID_relationship_36cases_NF1_TP53.csv")
CLINICAL_XLSX <- file.path(RAW_DIR, "clinical_table_s1_20260514.xlsx")
CT_PROP_CSV   <- file.path(RAW_DIR, "ct_prop_adjusted_31_samples.csv")

for (d in c(DATA_DIR, DATA_FULL_DIR, FIG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ------------ cohort constants ------------
N_COHORT        <- 36L     # the MRI table is read down to this row only
GROUP_COL       <- "group"
EXPR_TRANSFORM  <- "log2(x+1)"
PAIRWISE_ALPHA  <- 0.05
SIGNIF_MODE     <- "pval"  # "pval", "stars" or "both"

SAMPLE_ID_COL   <- "Sample ID"        # column names inside CLINICAL_XLSX
ENHANCEMENT_COL <- "Tumor Enhancement"

SCORE_COLS <- c(
  tardon_iqm = "Tardon hypoxia signature (IQM)",
  buffa_mean = "Buffa hypoxia signature (mean)",
  ragnum_iqm = "Ragnum hypoxia signature (IQM)"
)

GROUP_FILL_COLORS  <- c("#990000", "#2b8cbe")  # applied in factor-level order
GROUP_POINT_COLORS <- c("#990000", "#2b8cbe")

ENHANCEMENT_LEVELS <- c("Focal", "Lept", "Without")
ENH_FILL_COLORS <- c(
  Focal   = "#1b9e77",  # teal
  Lept    = "#d95f02",  # orange
  Without = "#7570b3"   # purple
)

# Cell-type columns, case must match ct_prop_adjusted_31_samples.csv exactly.
IMMUNE_RICH_COMPONENTS <- c("T_NK", "Myeloid", "Neutrophil", "B_cell")
IMMUNE_PRIMARY <- c(
  T_NK    = "T/NK ratio",
  Myeloid = "Myeloid ratio"
)

# Patient-level identifiers. Kept out of every table written to data/ so that the
# GitHub version carries the SP01..SP36 sample ids only.
PATIENT_ID_COLS <- c("PID", "Original_ID", "Multiomics_analysis_ID", "sample_X")

## ------------ small utilities ------------
need_pkg <- function(...) {
  for (p in c(...)) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Please install.packages('", p, "')", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Write a table twice: the full one to data_full/, an SP-keyed one to data/.
#'
#' The two differ only in the patient-identifier columns, so a reader of the
#' GitHub version sees the same numbers under their published sample ids.
write_tiered <- function(df, filename, subdir = NULL) {
  full_dir <- if (is.null(subdir)) DATA_FULL_DIR else file.path(DATA_FULL_DIR, subdir)
  pub_dir  <- if (is.null(subdir)) DATA_DIR else file.path(DATA_DIR, subdir)
  dir.create(full_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)

  utils::write.csv(df, file.path(full_dir, filename), row.names = FALSE)
  pub <- df[, setdiff(names(df), PATIENT_ID_COLS), drop = FALSE]
  utils::write.csv(pub, file.path(pub_dir, filename), row.names = FALSE)
  invisible(file.path(pub_dir, filename))
}

#' Write a table that holds no patient identifier to data/ only.
write_public <- function(df, filename, subdir = NULL) {
  pub_dir <- if (is.null(subdir)) DATA_DIR else file.path(DATA_DIR, subdir)
  dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(df, file.path(pub_dir, filename), row.names = FALSE)
  invisible(file.path(pub_dir, filename))
}

write_lines_public <- function(lines, filename, subdir = NULL) {
  pub_dir <- if (is.null(subdir)) DATA_DIR else file.path(DATA_DIR, subdir)
  dir.create(pub_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(lines, file.path(pub_dir, filename))
  invisible(file.path(pub_dir, filename))
}

read_gene_list <- function(path) {
  g <- utils::read.csv(path, stringsAsFactors = FALSE)[, 1]
  g <- unique(trimws(as.character(g)))
  g[nzchar(g) & !is.na(g)]
}

#' Cohort metadata: Original_ID / PID / SP / group, read from the ID map.
prepare_cohort_metadata <- function(meta_csv = ID_MAP_CSV, group_col = "Group") {
  id_tab <- utils::read.csv(meta_csv, stringsAsFactors = FALSE)
  req <- c("Original_ID", "ID_in_network", "PID", group_col)
  if (!all(req %in% names(id_tab))) {
    stop("ID map must contain: ", paste(req, collapse = ", "))
  }
  id_tab$orig_chr <- trimws(as.character(id_tab$Original_ID))
  id_tab$PID      <- trimws(as.character(id_tab$PID))
  id_tab$SP       <- paste0("S", trimws(as.character(id_tab$ID_in_network)))
  id_tab$group    <- trimws(as.character(id_tab[[group_col]]))
  id_tab
}

#' MRI tumour enhancement for the first N_COHORT rows of the clinical table.
read_mri_enhancement <- function(xlsx = CLINICAL_XLSX, sheet = 1L, n_cohort = N_COHORT) {
  need_pkg("readxl")
  mri <- as.data.frame(readxl::read_excel(xlsx, sheet = sheet))
  if (!all(c(SAMPLE_ID_COL, ENHANCEMENT_COL) %in% names(mri))) {
    stop("MRI table must contain: ",
         paste(c(SAMPLE_ID_COL, ENHANCEMENT_COL), collapse = ", "))
  }
  if (nrow(mri) < n_cohort) {
    stop("MRI table has fewer than ", n_cohort, " rows.")
  }
  mri <- mri[seq_len(n_cohort), , drop = FALSE]
  mri$SP <- trimws(as.character(mri[[SAMPLE_ID_COL]]))
  mri$tumor_enhancement <- trimws(as.character(mri[[ENHANCEMENT_COL]]))
  mri[, c("SP", "tumor_enhancement"), drop = FALSE]
}

#' Cell-type proportions plus the immune_total column, keyed by Original_ID.
read_ct_proportions <- function(csv = CT_PROP_CSV) {
  ct <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  # The first column holds the X-prefixed sample id but may have an empty
  # header, so index by position rather than by name.
  ct$sample_X    <- trimws(as.character(ct[[1L]]))
  ct$Original_ID <- sub("^X", "", ct$sample_X)

  miss <- setdiff(IMMUNE_RICH_COMPONENTS, names(ct))
  if (length(miss)) {
    warning("Cell-type CSV missing immune component columns: ",
            paste(miss, collapse = ", "), call. = FALSE)
  }
  present <- intersect(IMMUNE_RICH_COMPONENTS, names(ct))
  if (!length(present)) {
    stop("None of IMMUNE_RICH_COMPONENTS found in cell-type CSV.")
  }
  ct$immune_total <- rowSums(ct[, present, drop = FALSE], na.rm = TRUE)
  attr(ct, "present_components") <- present
  ct
}

#' Labels for the immune columns actually present in `df`.
immune_columns_and_labels <- function(df, present_components) {
  cols <- c(names(IMMUNE_PRIMARY),
            if ("immune_total" %in% names(df)) "immune_total" else NULL)
  cols <- cols[cols %in% names(df)]
  labels <- c(
    IMMUNE_PRIMARY,
    immune_total = sprintf("Immune total (%s)",
                           paste(present_components, collapse = " + "))
  )
  list(cols = cols, labels = labels[cols])
}

short_enhancement_label <- function(x) {
  x <- trimws(as.character(x))
  out <- rep(NA_character_, length(x))
  out[grepl("Focal", x, ignore.case = TRUE)] <- "Focal"
  out[grepl("Leptomeningeal", x, ignore.case = TRUE)] <- "Lept"
  out[grepl("^Without$", x, ignore.case = TRUE)] <- "Without"
  out
}

## ------------ statistics ------------
#' Interquartile mean: the mean of the values between the 25th and 75th
#' percentile, which limits the influence of individual outlying genes.
iqm <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) {
    return(NA_real_)
  }
  q1 <- unname(stats::quantile(x, 0.25, type = 7, names = FALSE))
  q3 <- unname(stats::quantile(x, 0.75, type = 7, names = FALSE))
  mid <- x[x >= q1 & x <= q3]
  if (!length(mid)) {
    return(NA_real_)
  }
  mean(mid)
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (!length(x) || !length(y)) {
    return(NA_real_)
  }
  dom <- 0
  for (xi in x) {
    dom <- dom + sum(xi > y) - sum(xi < y)
  }
  dom / (length(x) * length(y))
}

compare_two_groups <- function(vals, grp, g1, g2) {
  a <- vals[grp == g1]
  b <- vals[grp == g2]
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2L || length(b) < 2L) {
    return(list(note = "Insufficient finite values per group for testing."))
  }
  wt <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE, conf.int = TRUE))
  list(
    n1 = length(a), n2 = length(b),
    median1 = stats::median(a, na.rm = TRUE),
    median2 = stats::median(b, na.rm = TRUE),
    median_diff_g1_minus_g2 = stats::median(a, na.rm = TRUE) - stats::median(b, na.rm = TRUE),
    wilcox_W = unname(wt$statistic),
    wilcox_p = wt$p.value,
    conf_int_low = wt$conf.int[1],
    conf_int_high = wt$conf.int[2],
    cliffs_delta_g1_vs_g2 = cliffs_delta(a, b)
  )
}

compare_kruskal <- function(vals, grp) {
  ok <- is.finite(vals) & !is.na(grp) & nzchar(trimws(as.character(grp)))
  v <- vals[ok]
  g <- factor(trimws(as.character(grp[ok])))
  if (length(unique(g)) < 2L) {
    return(data.frame(
      n = length(v), n_groups = length(unique(g)),
      statistic = NA_real_, p_value = NA_real_,
      note = "Fewer than two enhancement groups with data.",
      stringsAsFactors = FALSE
    ))
  }
  kt <- stats::kruskal.test(v ~ g)
  data.frame(
    n = length(v),
    n_groups = length(levels(g)),
    statistic = unname(kt$statistic),
    p_value = kt$p.value,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

pairwise_wilcox_enhancement <- function(vals, grp) {
  ok <- is.finite(vals) & !is.na(grp) & nzchar(trimws(as.character(grp)))
  v <- vals[ok]
  g <- factor(trimws(as.character(grp[ok])))
  if (length(unique(g)) < 2L) {
    return(data.frame())
  }
  pw <- stats::pairwise.wilcox.test(v, g, p.adjust.method = "BH", exact = FALSE)
  pmat <- pw$p.value
  if (is.null(pmat)) {
    return(data.frame())
  }
  if (is.matrix(pmat)) {
    idx <- which(!is.na(pmat), arr.ind = TRUE)
    data.frame(
      group1 = rownames(pmat)[idx[, 1L]],
      group2 = colnames(pmat)[idx[, 2L]],
      p_adj_bh = as.numeric(pmat[idx]),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      group1 = names(pmat),
      group2 = NA_character_,
      p_adj_bh = as.numeric(pmat),
      stringsAsFactors = FALSE
    )
  }
}

lookup_pair_pval <- function(pw_sc, a, b) {
  if (is.null(pw_sc) || !nrow(pw_sc)) return(NA_real_)
  hit <- pw_sc[
    (pw_sc$group1 == a & pw_sc$group2 == b) |
      (pw_sc$group1 == b & pw_sc$group2 == a),
  ]
  if (!nrow(hit)) NA_real_ else hit$p_adj_bh[1]
}

spearman_safe <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) {
    return(list(n = sum(ok), rho = NA_real_, p_value = NA_real_,
                note = "n<3 finite pairs"))
  }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(n = sum(ok), rho = unname(ct$estimate), p_value = ct$p.value, note = NA_character_)
}

## ------------ formatting and plot styling ------------
format_pval <- function(p) {
  if (!is.finite(p)) return("p = NA")
  if (p < 1e-4) return(sprintf("p = %.1e", p))
  sprintf("p = %.3g", p)
}

pval_stars <- function(p) {
  if (!is.finite(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  "ns"
}

signif_annotation <- function(p, mode = SIGNIF_MODE) {
  switch(
    mode,
    pval  = format_pval(p),
    stars = pval_stars(p),
    both  = paste0(format_pval(p), " ", pval_stars(p))
  )
}

#' Shared plot theme. `angle_x = 45` for the grouped box plots, `0` for scatters;
#' `legend` is "none" for the box plots and "right" for the scatters.
base_theme <- function(angle_x = 45, legend = "none") {
  x_text <- if (angle_x == 0) {
    ggplot2::element_text(size = 9)
  } else {
    ggplot2::element_text(angle = angle_x, hjust = 1, size = 9)
  }
  ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.4),
      strip.text       = ggplot2::element_text(size = 9),
      axis.text.x      = x_text,
      axis.text.y      = ggplot2::element_text(size = 9),
      axis.title       = ggplot2::element_text(size = 10),
      plot.title       = ggplot2::element_text(size = 10, hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(size = 8, hjust = 0.5, color = "grey30"),
      legend.position  = legend,
      legend.title     = ggplot2::element_text(size = 9),
      legend.text      = ggplot2::element_text(size = 8)
    )
}

#' Box plot of one numeric column across the enhancement groups, with the
#' BH-adjusted pairwise p-values drawn as brackets.
#' Outliers are never drawn by the box layer: the jitter layer already shows
#' every sample, so outlier.shape = NA avoids plotting those points twice.
plot_by_enhancement <- function(df, value_col, value_label, out_prefix,
                                pw_sc = NULL, title = value_label,
                                dpi = 150) {
  need_pkg("ggplot2", "ggsignif")
  plot_df <- df[
    is.finite(df[[value_col]]) &
      !is.na(df$tumor_enhancement_short) &
      nzchar(df$tumor_enhancement_short),
    ,
    drop = FALSE
  ]
  if (!nrow(plot_df)) {
    warning("No finite values to plot for ", value_col, call. = FALSE)
    return(invisible(NULL))
  }

  lv <- ENHANCEMENT_LEVELS[ENHANCEMENT_LEVELS %in% unique(plot_df$tumor_enhancement_short)]
  plot_df$tumor_enhancement_short <- factor(plot_df$tumor_enhancement_short, levels = lv)

  n_per_level <- table(plot_df$tumor_enhancement_short)
  x_tick_labels <- sprintf("%s (n=%d)", names(n_per_level), as.integer(n_per_level))
  names(x_tick_labels) <- names(n_per_level)

  my_comparisons <- utils::combn(lv, 2L, simplify = FALSE)
  ann_txt <- vapply(
    my_comparisons,
    function(pair) signif_annotation(lookup_pair_pval(pw_sc, pair[1], pair[2])),
    character(1)
  )

  y_vals  <- plot_df[[value_col]]
  y_max   <- max(y_vals, na.rm = TRUE)
  y_min   <- min(y_vals, na.rm = TRUE)
  y_range <- y_max - y_min
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1
  step        <- y_range * 0.12
  y_positions <- y_max + step * seq_along(my_comparisons)
  y_top       <- y_max + step * (length(my_comparisons) + 0.8)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data$tumor_enhancement_short,
    y = .data[[value_col]]
  )) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data$tumor_enhancement_short),
                          width = 0.55, outlier.shape = NA,
                          size = 0.5, color = "black", show.legend = FALSE) +
    ggplot2::geom_jitter(position = ggplot2::position_jitter(width = 0.15, height = 0),
                         size = 1.6, alpha = 0.8, color = "black") +
    ggplot2::scale_fill_manual(values = ENH_FILL_COLORS) +
    ggplot2::scale_x_discrete(labels = x_tick_labels) +
    ggsignif::geom_signif(
      comparisons = my_comparisons,
      annotations = ann_txt,
      y_position  = y_positions,
      tip_length  = 0.015,
      textsize    = 3,
      vjust       = -0.2,
      size        = 0.4
    ) +
    ggplot2::expand_limits(y = y_top) +
    ggplot2::labs(
      title    = title,
      subtitle = sprintf("MRI Tumor Enhancement (n = %s, first %s cases)",
                         nrow(plot_df), N_COHORT),
      x = NULL,
      y = value_label
    ) +
    base_theme()

  w <- 3.5 + 0.3 * max(0, length(lv) - 2)
  h <- 3.2 + 0.2 * length(my_comparisons)
  ggplot2::ggsave(paste0(out_prefix, ".pdf"), p, width = w, height = h)
  ggplot2::ggsave(paste0(out_prefix, ".png"), p, width = w, height = h, dpi = dpi)
  invisible(p)
}
