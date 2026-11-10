# Shared helpers for the two scripts in this directory.

# ---- p value formatting ----
# Rendered by ggtext::element_markdown, so *p* shows as an italic p.
# Small p values use scientific notation written as "2.2e-4" (no leading zero in
# the exponent).
fmt_p_sci <- function(p, digits = 2) {
  if (is.na(p)) return("NA")
  if (p < 0.001) sub("e-0", "e-", sprintf(paste0("%.", digits - 1, "e"), p), fixed = TRUE)
  else formatC(signif(p, digits), format = "g")
}

fmt_p <- function(p, italic = TRUE) {
  pl <- if (italic) "*p*" else "p"
  if (is.na(p)) paste0(pl, " = NA")
  else if (p < 0.0001) paste0(pl, " < 0.0001")
  else paste0(pl, " = ", formatC(signif(p, 2), format = "g"))
}

# ---- Mut / WT normalisation ----
# NF1 and TP53 are character columns in Table S1 and missing values are the
# literal string "NA".
mut_status <- function(x) {
  x <- trimws(as.character(x))
  ifelse(x %in% c("Mut", "Mutant", "MUT"), "Mut",
  ifelse(x %in% c("WT", "Wildtype", "wt"), "WT", NA_character_))
}

# ---- sentinel check for the validation cohort ----
# Run this before any PL vs TL comparison on the validation cohort. SP41, SP80
# and SP82 are the only three samples where the current table and the legacy
# files disagree: SP41 is a stale PL in the legacy survival files, and SP80 and
# SP82 were swapped in the 2026-09-14 update (SP80 TL -> PL, SP82 PL -> TL). A
# mismatch means the input is out of date, or the sample IDs use the old
# numbering (old SP70-SP82 = current SP74-SP86). Run it on the unfiltered table:
# SP80 has no OS and drops out of any survival subset.
VALIDATION_SENTINELS <- c(SP41 = "TL", SP80 = "PL", SP82 = "TL")

check_validation_sentinels <- function(data,
                                       id_col = "Sample ID",
                                       group_col = "Group/Predict Group",
                                       stop_on_fail = TRUE) {
  ids <- trimws(as.character(data[[id_col]]))
  grp <- trimws(as.character(data[[group_col]]))
  norm <- ifelse(grepl("^(PL|Group1|subtype.?1)", grp, ignore.case = TRUE), "PL",
          ifelse(grepl("^(TL|TML|Group2|subtype.?2)", grp, ignore.case = TRUE), "TL", NA))

  res <- data.frame(sample = names(VALIDATION_SENTINELS),
                    expected = unname(VALIDATION_SENTINELS),
                    observed = NA_character_,
                    status = NA_character_,
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(res))) {
    hit <- which(ids == res$sample[i])
    if (length(hit) == 0) {
      res$status[i] <- "ABSENT"
    } else {
      res$observed[i] <- paste(unique(norm[hit]), collapse = "/")
      res$status[i] <- if (identical(res$observed[i], res$expected[i])) "OK" else "MISMATCH"
    }
  }

  cat("validation sentinel check:\n")
  for (i in seq_len(nrow(res))) {
    cat(sprintf("  %-5s expected %-2s observed %-6s [%s]\n",
                res$sample[i], res$expected[i],
                ifelse(is.na(res$observed[i]), "-", res$observed[i]), res$status[i]))
  }

  bad <- res[res$status == "MISMATCH", ]
  if (nrow(bad) > 0) {
    msg <- paste0("validation sentinels do not match; the input is wrong ",
                  "(stale table, or old SP_ID numbering): ",
                  paste0(bad$sample, " expected ", bad$expected,
                         " observed ", bad$observed, collapse = "; "))
    if (stop_on_fail) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  if (any(res$status == "ABSENT")) {
    cat(sprintf("  note: %s not present in this data, not checked\n",
                paste(res$sample[res$status == "ABSENT"], collapse = ", ")))
  }
  invisible(res)
}

# ---- stacked mutation proportion plot with a Fisher exact test ----
# tab: 2x2 contingency table, rows = c("Mut","WT"), columns = the two groups in
# display order.
# style:
#   "legend_title"  -- With/Without legend on the right, group counts and the
#                      Fisher p in the title block, no text inside the plot
#   "compact_inbar" -- no legend, only the P value on top, y axis reads
#                      "<gene> mutation Proportion (%)"
plot_mut_prop_fisher <- function(tab,
                                 gene,
                                 group_labels  = colnames(tab),
                                 style         = c("legend_title", "compact_inbar"),
                                 legend_labels = c("With", "Without"),   # mutant, wildtype
                                 fill_cols     = c("#fec44f", "#2c7fb8"),
                                 title         = NULL,
                                 bar_width     = 0.9,
                                 x_text_angle  = 60,
                                 base_size     = 9,
                                 title_size    = 9.5,
                                 subtitle_size = 7) {
  style <- match.arg(style)
  stopifnot(identical(dim(tab), c(2L, 2L)))
  ft   <- fisher.test(tab)
  n_by <- colSums(tab)
  mut  <- tab["Mut", ]

  # as.table() lets both a matrix and a table object be melted
  plot_dat <- as.data.frame(as.table(as.matrix(tab)), stringsAsFactors = FALSE)
  names(plot_dat) <- c("Status", "Group", "n")
  plot_dat <- plot_dat %>%
    group_by(Group) %>%
    mutate(Proportion = n / sum(n)) %>%
    ungroup() %>%
    mutate(
      Status_lab = factor(ifelse(Status == "Mut", legend_labels[1], legend_labels[2]),
                          levels = legend_labels),
      Group_lab  = factor(unname(group_labels[match(Group, colnames(tab))]),
                          levels = unname(group_labels))
    )

  # "<group> <mut>/<group total> (<mutant percent>)", e.g. "SP DMG 15/36 (42%)"
  one_grp <- function(i) sprintf("%s %d/%d (%.0f%%)", group_labels[i], mut[i], n_by[i],
                                 100 * mut[i] / n_by[i])

  if (style == "compact_inbar") {
    plot_dat$Pct <- 100 * plot_dat$Proportion
    return(
      ggplot(plot_dat, aes(x = Group_lab, y = Pct, fill = Status_lab)) +
        geom_col(width = bar_width) +
        scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                           expand = expansion(mult = c(0, 0.02))) +
        scale_fill_manual(values = unname(fill_cols)) +
        labs(x = '', y = sprintf("*%s* mutation Proportion (%%)", gene),
             title    = paste0("*P* = ", fmt_p_sci(ft$p.value)),
             # too narrow for one line, so this is split in two
             subtitle = paste0(one_grp(1), "<br>vs ", one_grp(2))) +
        theme_classic(base_size = base_size) +
        theme(
          legend.position = "none",
          axis.text.x  = element_text(color = "black", size = base_size - 0.5,
                                      angle = x_text_angle,
                                      hjust = if (x_text_angle == 0) 0.5 else 1),
          axis.text.y  = element_text(color = "black", size = base_size - 0.5),
          axis.title.y = ggtext::element_markdown(color = "black", size = base_size),
          plot.title    = ggtext::element_markdown(size = base_size, colour = "black",
                                                   hjust = 0.5, margin = margin(b = 1)),
          plot.subtitle = ggtext::element_markdown(size = subtitle_size, colour = "black",
                                                   hjust = 0.5, margin = margin(b = 3)),
          plot.margin   = margin(t = 4, r = 4, b = 2, l = 2)
        )
    )
  }

  ggplot(plot_dat, aes(x = Group_lab, y = Proportion, fill = Status_lab)) +
    geom_col(width = bar_width) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = unname(fill_cols)) +
    labs(
      x = '', y = 'Proportion', fill = '',
      title    = if (is.null(title)) paste0(gene, ' mutation status') else title,
      # <br> is rendered as a line break by element_markdown; a narrow panel
      # clips a single line, so this is always split in three
      subtitle = paste0(one_grp(1), "<br>vs ", one_grp(2),
                        "<br>Fisher's exact ", fmt_p(ft$p.value))
    ) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.32, "cm"),
      legend.text      = element_text(size = base_size - 1.5),
      legend.margin    = margin(l = -4),
      axis.text.x  = element_text(color = "black", size = base_size,
                                  angle = x_text_angle, hjust = 1),
      axis.text.y  = element_text(color = "black", size = base_size - 1),
      axis.title.y = element_text(color = "black", size = base_size + 1),
      plot.title    = ggtext::element_markdown(size = title_size, face = "bold",
                                               colour = "black", hjust = 0.5,
                                               margin = margin(b = 1)),
      plot.subtitle = ggtext::element_markdown(size = subtitle_size, colour = "black",
                                               hjust = 0.5, margin = margin(b = 4)),
      plot.margin = margin(t = 4, r = 2, b = 2, l = 2)
    )
}

# ---- save with the panel, not the canvas, at a fixed physical size ----
# What has to match across panels is the physical size of the bars, not the size
# of the canvas. With a fixed canvas the panel is squeezed by the number of title
# lines, the height of rotated axis labels and the width of the legend, so the
# bars come out a different size in every figure. This fixes the panel instead
# and lets the canvas grow around it.
#   with two categories and bar_width = 0.9, one bar is panel_width / 2 * 0.9
save_fixed_panel <- function(p, filename,
                             panel_width  = 1.7,   # inch
                             panel_height = 2.05,  # inch
                             margin = 0.06,        # canvas margin, inch
                             verbose = TRUE, ...) {
  g  <- ggplot2::ggplotGrob(p)
  wi <- g$layout$l[g$layout$name == "panel"]
  hi <- g$layout$t[g$layout$name == "panel"]
  g$widths[wi]  <- grid::unit(panel_width,  "in")
  g$heights[hi] <- grid::unit(panel_height, "in")

  total_w <- sum(grid::convertWidth(g$widths,   "in", valueOnly = TRUE)) + 2 * margin
  total_h <- sum(grid::convertHeight(g$heights, "in", valueOnly = TRUE)) + 2 * margin
  if (verbose) {
    cat(sprintf('  panel %.2f x %.2f in -> canvas %.2f x %.2f in (one bar %.3f in): %s\n',
                panel_width, panel_height, total_w, total_h,
                panel_width / 2 * 0.9, basename(filename)))
  }
  ggplot2::ggsave(filename, g, width = total_w, height = total_h, ...)
  invisible(c(total_w = total_w, total_h = total_h))
}
