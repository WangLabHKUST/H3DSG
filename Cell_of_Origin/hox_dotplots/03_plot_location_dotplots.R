# =============================================================================
# hox_dotplots / 03_plot_location_dotplots.R -- dot plots of the four position gene sets
#
# Input : ../data/hox_dotplots/location_dotplot_by_group.tsv
#         ../data/hox_dotplots/location_dotplot_by_sample.tsv
# Output: hox_location_dotplot_by_group.pdf     7.5 x 3.5 in
#         hox_location_dotplot_by_sample.pdf    13.65 x 2.16 cm, no text (see below)
#
# Rows are the four gene sets (Thalamus, Brainstem, Cervical, Thoracic: the genes that mark
# these levels), columns the immune, stromal and malignant cells. Colour = average
# expression of the gene set, dot area = percentage of cells expressing at least one gene
# of the set.
#
# by group : the public thalamic and pontine DMG cells and the H3-DSG-PL and H3-DSG-TL
#            samples of this cohort split into Cervical and Thoracic (the group of TL_4 is
#            set in 00_config.R). Average expression is clipped at SPLIT_CMAX.
# by sample: each of the 19 samples (one panel each) in the order of the resection site,
#            from cervical to thoracic. Colour bar 0 to COLOR_MAX. The figure is drawn
#            without axis text, titles or legends, and with small dots (dot.scale 2.2),
#            to fit the small page; the sample labels and the legends are added when the
#            figure is assembled. Sample order, left to right: column panel_order of the table.
#
# This script reads the two tables only and opens no object.
# =============================================================================
source("00_config.R")
source("dotplot_df.R")
suppressPackageStartupMessages(library(ggplot2))

# ---- by group ------------------------------------------------------------------------
g <- read_table("location_dotplot_by_group.tsv")
group_order <- c("Thalamus", "Brainstem", "Cervical-PL", "Cervical-TL", "Thoracic-PL", "Thoracic-TL")
stopifnot(setequal(g$group, group_order), setequal(g$celltype, CELLTYPES), setequal(g$gene_set, GENE_SET_LEVELS))
group_df <- data.frame(features.plot = factor(g$gene_set, levels = rev(GENE_SET_LEVELS)),   # y axis, bottom to top
                       id = factor(g$celltype, levels = CELLTYPES),
                       avg.exp = g$avg_exp, pct.exp = g$pct_exp,
                       orig.ident = factor(g$group, levels = group_order))
p_group <- DotPlotDF(
  group_df,
  group.by = "id", facet.by = "orig.ident", facet.nrow = 1, flip = TRUE,
  scale = FALSE, log.transform = FALSE, scale.by = "size", dot.scale = 8,
  col.max = SPLIT_CMAX,                      # average expression clipped at SPLIT_CMAX
  size.breaks = seq(0, 40, 10), size.labels = NULL, size.limits = c(0, NA), scale.max = NA,
  color.scale = avg_color_scale(cmax = NA_real_))
ggsave("hox_location_dotplot_by_group.pdf", p_group, width = 7.5, height = 3.5, dpi = 300, limitsize = FALSE)
message("saved hox_location_dotplot_by_group.pdf")

# ---- by sample -------------------------------------------------------------------------
s <- read_table("location_dotplot_by_sample.tsv")
stopifnot(length(unique(s$sid)) == 19, setequal(s$celltype, CELLTYPES), setequal(s$gene_set, GENE_SET_LEVELS))
panel_levels <- unique(s$sid[order(s$panel_order)])                    # left to right
sample_df <- data.frame(features.plot = factor(s$gene_set, levels = rev(GENE_SET_LEVELS)),
                        id = factor(s$celltype, levels = CELLTYPES),
                        avg.exp = s$avg_exp, pct.exp = s$pct_exp,
                        orig.ident = factor(s$sid, levels = panel_levels))
sample_df <- sample_df[order(sample_df$orig.ident, sample_df$id), , drop = FALSE]
# every panel has its own three x ticks: the group column is "sample <US> cell type", which is
# unique over the figure, with facet.scales = "free_x"; the sample prefix is cut off again in the labels
ID_N_SEP <- "\u001f"
sample_df$id_n <- paste0(as.character(sample_df$orig.ident), ID_N_SEP, as.character(sample_df$id))
sample_df$id_n <- factor(sample_df$id_n, levels = unique(sample_df$id_n))
p_sample <- DotPlotDF(
  sample_df,
  group.by = "id_n", facet.by = "orig.ident", facet.nrow = 1,
  facet.drop = FALSE, facet.scales = "free_x", flip = TRUE, sort.levels = TRUE,
  scale = FALSE, log.transform = FALSE, scale.by = "size", dot.scale = 2.2,
  size.limits = c(0, NA), scale.max = NA, size.breaks = NULL, size.labels = NULL,
  color.scale = avg_color_scale()
) + scale_x_discrete(labels = function(x) sub(paste0("^.*", ID_N_SEP), "", x))
ggsave("hox_location_dotplot_by_sample.pdf", p_sample + nude_theme(),
       width = 13.65, height = 2.16, units = "cm", dpi = 300, limitsize = FALSE)
message("saved hox_location_dotplot_by_sample.pdf")
