# Extended Data Fig. 2f: pseudobulk expression of the 39 HOX genes in single-cell RNA-seq, log2(CPM + 1).
# Rows: thalamus DMG (brain), pons DMG (brain), H3-DSG-TL and H3-DSG-PL (spinal cord). Columns: HOXA, HOXB, HOXC, HOXD.
# One dot is one patient (brain DMG: one sample of a patient); box colour follows the paralog number (1 = anterior, 13 = posterior).
#
# Input (data/):
#   ThalamusDMG_scRNA.HOX.txt          pseudobulk counts and CPM of the HOX genes, 4 thalamus DMG samples (public single-cell data
#                                      GSE210568; sample IDs are those of the public data set)
#   H3PonsDMG_scRNA.HOX.txt            the same for 13 H3 K27M pons DMG samples of that data set
#   hox_pseudobulk_by_patient.tsv      the same for the 15 spinal cord patients, made by hox_pseudobulk_build.R
# Output: hox_pseudobulk_boxplot.pdf
#
# The brain tables were made with the same procedure as hox_pseudobulk_build.R (sum of counts of the 39 HOX genes over the cells of a
# sample, divided by the total counts of the sample) from the public 10x data; the raw matrices are not part of this repository.
# Spinal cord groups: the group column of the project metadata (PL / TL), as in hox_pseudobulk_build.R. The Source Data of the
# figure group the patients by Table S1.
# Run from the Cell_of_Origin/ directory:  Rscript hox_pseudobulk_boxplot.R

library(ggplot2)
library(cowplot)

clusters <- c("HOXA", "HOXB", "HOXC", "HOXD")

# ---- data ----
read_brain <- function(f) {
  d <- read.delim(f)
  d <- d[order(d$HOX, d$int), ]
  d$gene_name <- factor(d$gene, levels = d$gene[!duplicated(d$gene)])
  d
}
thalamus <- read_brain("data/ThalamusDMG_scRNA.HOX.txt")
pons     <- read_brain("data/H3PonsDMG_scRNA.HOX.txt")
stopifnot(length(unique(thalamus$SID)) == 4, length(unique(pons$SID)) == 13)

spinal <- read.delim("data/hox_pseudobulk_by_patient.tsv", na.strings = "NA")
spinal <- spinal[!is.na(spinal$CPM), ]                                    # genes absent from the reference of one sample are missing, not 0
spinal <- spinal[order(spinal$HOX, spinal$int), ]
spinal$gene_name <- factor(spinal$gene, levels = unique(spinal$gene))
n_patients <- tapply(spinal$patient_id, spinal$group, function(x) length(unique(x)))
stopifnot(identical(as.vector(n_patients[c("PL", "TL")]), c(8L, 7L)))
spinal_tl <- spinal[spinal$group == "TL", ]
spinal_pl <- spinal[spinal$group == "PL", ]

# ---- one panel: one HOX cluster of one row ----
box_panel <- function(d, h, ylim = NULL, scale_limits = FALSE) {
  p <- ggplot(d[d$HOX == h, ], aes(x = gene_name, y = log2(CPM + 1))) +
    geom_boxplot(width = 0.5, aes(color = int, fill = int), outlier.shape = NA, show.legend = FALSE, alpha = 0.5) +
    geom_jitter(aes(color = int), pch = 20, size = 0.1, position = position_jitter(width = 0.1, seed = 42), show.legend = FALSE) +
    scale_color_gradient2(high = "#f7fcfd", mid = "#8c96c6", low = "#4d004b", midpoint = 7, limits = c(1, 13)) +
    scale_fill_gradient2(high = "#f7fcfd", mid = "#8c96c6", low = "#4d004b", midpoint = 7, limits = c(1, 13)) +
    theme_bw() +
    theme(axis.text = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
          panel.grid = element_blank(), plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}
ylim_spinal <- c(0, max(log2(spinal$CPM + 1)) * 1.05)                     # one y range for the spinal cord rows, so that they compare

# thalamus: y axis 0-10 (values above 10 are not drawn); pons: y axis free per panel; spinal cord: common y axis
thalamus_panels <- lapply(clusters, function(h) box_panel(thalamus, h) + scale_y_continuous(limits = c(0, 10)))
pons_panels     <- lapply(clusters, function(h) box_panel(pons, h))
tl_panels       <- lapply(clusters, function(h) box_panel(spinal_tl, h, ylim = ylim_spinal))
pl_panels       <- lapply(clusters, function(h) box_panel(spinal_pl, h, ylim = ylim_spinal))

# ---- labels: column titles, row names, gene names below the last row ----
txt <- function(label, size = 7, colour = "black", hjust = 0.5, x = 0.5) ggdraw() + draw_label(label, size = size, colour = colour, hjust = hjust, x = x)
gene_labels <- lapply(clusters, function(h)
  ggplot(data.frame(gene = factor(levels(spinal$gene_name)[grepl(h, levels(spinal$gene_name))], levels = levels(spinal$gene_name)), y = 1),
         aes(gene, y)) + geom_blank() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "grey30", size = 6), axis.text.y = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), rect = element_blank(),
          plot.margin = margin(r = 8)))

plots <- c(list(NULL), lapply(clusters, txt, size = 8),
           list(txt("H3 K27M\nThalamus DMG\n(Brain)", hjust = 1, x = 0.95)), thalamus_panels,
           list(txt("H3 K27M\nPons DMG\n(Brain)", hjust = 1, x = 0.95)), pons_panels,
           list(txt("H3-DSG-TL\n(Spinal cord)", colour = "#2b8cbe", hjust = 1, x = 0.95)), tl_panels,
           list(txt("H3-DSG-PL\n(Spinal cord)", colour = "#990000", hjust = 1, x = 0.95)), pl_panels,
           list(NULL), gene_labels)
g <- plot_grid(plotlist = plots, ncol = 5, rel_widths = c(0.55, 1, 1, 1, 1), rel_heights = c(0.12, 1, 1, 1, 1, 0.45))
ggsave("hox_pseudobulk_boxplot.pdf", g, width = 437 / 72, height = 335 / 72)    # panels 96 pt wide, as in the original figure
