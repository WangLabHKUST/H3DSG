# =============================================================================
# Shared by every script in figures/: palettes, orders, labels, theme, saving.
#
# Scripts load it with  source("panel_helpers.R")  and are run from this directory.
# Nothing here reads data or opens an object.
#
# Samples are labelled by `renamed_id` (never `sid`). Patients are grouped by `patient_id`:
# TL_6 .. TL_10 are five sites of one tumour and count as one patient in any per-patient
# statistic.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# ---- major cell types --------------------------------------------------------------
# `cell_type` in data/cell_level/cell_annotation.tsv.gz carries the object's own labels
# (tumor, oligo, ...); figures use the labels below.
CELLTYPE_LABELS <- c(tumor = "Malignant", myeloid = "Myeloid", neutrophil = "Neutrophil",
                     "T/NK" = "T/NK", B = "B", vascular = "Vascular", oligo = "Oligo")

# Display order: legend order, and the order of a stack from the bottom (Malignant) up.
CELLTYPE_ORDER <- unname(CELLTYPE_LABELS)

# Stacked bars and dot plots.
CELLTYPE_COLORS <- c(Malignant = "#c63932", Myeloid = "#8e69b7", Neutrophil = "#babe2b",
                     "T/NK" = "#228b22", B = "#216096", Vascular = "#ff8c00", Oligo = "#da80be")

# The UMAP of the major cell types uses a slightly different set of the same hues.
CELLTYPE_COLORS_UMAP <- c(Malignant = "#C1272D", Myeloid = "#8E66A2", Neutrophil = "#B5B35C",
                          "T/NK" = "#4E9A4E", B = "#3274A1", Vascular = "#E67E22", Oligo = "#D980B9")

# Raw `cell_type` value(s) -> ordered factor of display labels. A value that is not in
# CELLTYPE_LABELS stops the script rather than turning into NA.
celltype_label <- function(x) {
  x <- as.character(x)
  bad <- setdiff(unique(x), names(CELLTYPE_LABELS))
  if (length(bad)) stop("unknown cell_type: ", paste(bad, collapse = ", "))
  factor(unname(CELLTYPE_LABELS[x]), levels = CELLTYPE_ORDER)
}

# ---- malignant states (colours of the six states, as in the original analysis) -------
MALIGNANT_STATE_ORDER  <- c("G2M", "S", "OPC-like", "OC-like", "MES-like", "AC-like")
MALIGNANT_STATE_COLORS <- c("G2M" = "#c63932", "S" = "#e36b58", "OPC-like" = "#1f9a8e",
                            "OC-like" = "#3bc2b6", "MES-like" = "#ef2c8d", "AC-like" = "#ff5ca8")

# ---- groups --------------------------------------------------------------------------
GROUP_LEVELS <- c("PL", "TL")
GROUP_LABELS <- c(PL = "H3-DSG-PL", TL = "H3-DSG-TL")
GROUP_COLORS <- c(PL = "#990000", TL = "#2b8cbe")

# Samples in the order of the figures: PL before TL, then by the number in `renamed_id`
# (PL_1, PL_2, ..., PL_11; numeric, not alphabetical).
order_samples <- function(renamed_id) {
  ids <- unique(as.character(renamed_id))
  grp <- sub("_.*$", "", ids)
  num <- as.integer(sub("^[A-Z]+_", "", ids))
  ids[order(match(grp, GROUP_LEVELS), num)]
}

# ---- theme and saving ----------------------------------------------------------------
theme_panel <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(axis.text = element_text(colour = "black"), axis.title = element_text(colour = "black"),
          plot.title = element_text(colour = "black", hjust = 0.5), strip.background = element_blank())
}

# PDF and PNG next to the calling script; the PNG has an opaque white background.
save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(stem)
}

# ---- PL vs TL box plot ---------------------------------------------------------------
# One box and the points of one measurement. `d` has columns `group` (PL / TL) and `percent`.
# Box plots never draw outliers: every point is drawn once, by the jittered layer.
group_box <- function(d, title, p_value, ylab = NULL, breaks = waiver()) {
  d$group <- factor(d$group, levels = GROUP_LEVELS)
  ptxt <- sprintf("%.1e", p_value)
  ggplot(d, aes(group, percent, fill = group)) +
    geom_boxplot(width = 0.55, linewidth = 0.8, alpha = 0.9, colour = "black", outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.18, height = 0, seed = 1),
               size = 2.1, alpha = 0.95, colour = "black") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(breaks = breaks) +
    labs(x = NULL, y = ylab, title = bquote(atop(.(title), italic(P) ~ "=" ~ .(ptxt)))) +
    theme_panel() +
    theme(legend.position = "none", axis.text = element_text(colour = "black", size = 12),
          axis.title.y = element_text(size = 14))
}
