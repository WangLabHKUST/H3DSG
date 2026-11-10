# Shared helpers for the script in this directory.

# ---- p value formatting ----
# Rendered by ggtext::element_markdown, so *p* shows as an italic p.
fmt_p <- function(p, italic = TRUE) {
  pl <- if (italic) "*p*" else "p"
  if (is.na(p)) paste0(pl, " = NA")
  else if (p < 0.0001) paste0(pl, " < 0.0001")
  else paste0(pl, " = ", formatC(signif(p, 2), format = "g"))
}

# ---- Yes / No normalisation ----
# The treatment columns of Table S1 are character and missing values are the
# literal string "NA".
yesno <- function(x) {
  x <- trimws(as.character(x))
  ifelse(x %in% c("Yes", "yes", "YES"), "Yes",
  ifelse(x %in% c("No", "no", "NO"), "No", NA_character_))
}

# ---- sentinel check for the validation cohort ----
# Run this before any PL vs TL comparison that includes the validation cohort.
# SP41, SP80 and SP82 are the only three samples where the current table and the
# legacy files disagree: SP41 is a stale PL in the legacy survival files, and
# SP80 and SP82 were swapped in the 2026-09-14 update (SP80 TL -> PL,
# SP82 PL -> TL). A mismatch means the input is out of date, or the sample IDs
# use the old numbering (old SP70-SP82 = current SP74-SP86). Run it on the
# unfiltered table: SP80 has no OS and drops out of every survival subset.
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

# ---- survminer font customisation ----
# ggsurvplot returns a list of ggplots; this walks it and applies the fonts so
# that titles and axis labels are rendered through element_markdown, which is
# what lets the *p* in a title show as an italic p.
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL,
                              font.xtickslab = NULL, font.ytickslab = NULL){
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}

# ---- append the group size to a legend label ----
# ggsurvplot legend.labs accepts newlines, so (n=xx) goes on its own line.
#   legend_labs_with_n(c("H3-DSG-PL","H3-DSG-TL"), c(20, 28))
#   -> "H3-DSG-PL\n(n=20)"  "H3-DSG-TL\n(n=28)"
# The order of n must match the strata order of the survfit, which is the order
# of the levels of the grouping factor.
legend_labs_with_n <- function(labels, n, sep = "\n") {
  stopifnot(length(labels) == length(n))
  sprintf("%s%s(n=%d)", labels, sep, as.integer(n))
}
