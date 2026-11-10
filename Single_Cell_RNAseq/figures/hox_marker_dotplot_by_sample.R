# =============================================================================
# Expression of spinal-cord and brain position marker genes in each sample
#
# Input : ../data/panel_data/hox_dotplot_stats_by_sample_celltype.tsv
#         ../data/panel_data/hox_gene_annotation.tsv
#         ../data/panel_data/sample_resection_order.tsv
# Output: hox_marker_dotplot_by_sample.pdf / .png in this directory
#
# One panel per sample (19), in the order of the resection site from cervical to thoracic
# (caudal end of the tumour, then cranial end); columns are immune, stromal and malignant
# cells. Rows are the HOX genes of the brainstem, cervical and thoracic levels and the
# thalamic markers ZIC3, ZIC1, ZIC4 and OTX1; the gene name is coloured by the level it marks.
# Colour = mean log-normalised expression of the gene over all cells of the group; dot area
# = percentage of cells in which the gene is detected.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

CELLTYPE_LEVELS <- c("immune", "stromal", "malignant")
CELLTYPE_LABELS <- c(immune = "Immune", stromal = "Stromal", malignant = "Malignant")
LEVEL_COLORS <- c(Thalamus = "#35c1ee", Brainstem = "#e8000b", Cervical = "#e259c2", Thoracic = "#8b90c9")

hox <- read.delim("../data/panel_data/hox_dotplot_stats_by_sample_celltype.tsv")
ann <- read.delim("../data/panel_data/hox_gene_annotation.tsv")
smp <- read.delim("../data/panel_data/sample_resection_order.tsv", colClasses = c(patient_id = "character"))

genes <- ann[!is.na(ann$order_all_samples), ]
genes <- genes[order(genes$order_all_samples), ]                 # first level at the bottom of the plot
d <- hox[hox$gene %in% genes$gene, ]
d$gene <- factor(d$gene, levels = genes$gene)
d$celltype <- factor(d$celltype, levels = CELLTYPE_LEVELS, labels = CELLTYPE_LABELS[CELLTYPE_LEVELS])
# panel = sample, in resection order; the title is the patient ID (a multi-site patient has several panels)
smp <- smp[order(smp$resection_order), ]
d$panel <- factor(d$renamed_id, levels = smp$renamed_id)
panel_title <- setNames(smp$patient_id, smp$renamed_id)

p <- ggplot(d, aes(celltype, gene)) +
  geom_point(aes(size = pct_exp, colour = mean_log_expr)) +
  scale_colour_gradient(low = "white", high = "blue", limits = c(0, NA), name = "Avg. Exp") +
  scale_size(range = c(0, 5), name = "Pct. Exp", breaks = c(20, 40, 60)) +
  facet_grid(. ~ panel, labeller = labeller(panel = panel_title)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), panel.spacing = unit(0.15, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = LEVEL_COLORS[genes$location_group]))

save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
}
save_figure(p, "hox_marker_dotplot_by_sample", width = 14, height = 6.5)
