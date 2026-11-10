# =============================================================================
# UMAP of the cells that carry the H3 K27 mutant transcript
#
# Input : ../data/cell_level/cell_annotation.tsv.gz   Alt_freq, depth_Alt, depth_Total
#         ../data/cell_level/umap_coordinates.tsv.gz  umap_corr_1, umap_corr_2
# Output: umap_k27m_cells.pdf / .png in this directory
#
# Alt_freq = depth_Alt / depth_Total, the fraction of the reads over the H3 K27 position of a
# cell that carry the mutant allele. A cell is coloured when Alt_freq >= 0.4 (THR), whatever
# the number of reads; it is then coloured by Alt_freq. Every other cell is grey: cells with
# no read over the position (NA) and cells with Alt_freq < 0.4. All 132,343 cells are drawn,
# malignant or not. Same embedding as the UMAP of the major cell types (umap_corr_1/2).
# =============================================================================

source("panel_helpers.R")

THR <- 0.4

ann  <- read.delim(gzfile("../data/cell_level/cell_annotation.tsv.gz"), colClasses = c(patient_id = "character"))
umap <- read.delim(gzfile("../data/cell_level/umap_coordinates.tsv.gz"))
stopifnot(identical(ann$cell_id, umap$cell_id))

plot_df <- data.frame(cell_id = ann$cell_id, umap_1 = umap$umap_corr_1, umap_2 = umap$umap_corr_2,
                      cell_type = celltype_label(ann$cell_type), depth_Alt = ann$depth_Alt,
                      depth_Total = ann$depth_Total, Alt_freq = ann$Alt_freq)
plot_df$coloured <- !is.na(plot_df$Alt_freq) & plot_df$Alt_freq >= THR

bg <- plot_df[!plot_df$coloured, ]
fg <- plot_df[plot_df$coloured, ]
fg <- fg[order(fg$Alt_freq), ]                      # high frequency on top

# axis arrows in panel units, lower left corner
arrows <- list(
  annotation_custom(grid::segmentsGrob(x0 = unit(c(0.02, 0.02), "npc"), x1 = unit(c(0.12, 0.02), "npc"),
                                       y0 = unit(c(0.04, 0.04), "npc"), y1 = unit(c(0.04, 0.30), "npc"),
                                       arrow = arrow(length = unit(0.06, "in"), type = "closed"),
                                       gp = grid::gpar(lwd = 1.5, fill = "black")),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
  annotation_custom(grid::textGrob("UMAP1", x = unit(0.07, "npc"), y = unit(0.0, "npc"), vjust = 1,
                                   gp = grid::gpar(fontsize = 8)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
  annotation_custom(grid::textGrob("UMAP2", x = unit(0.0, "npc"), y = unit(0.17, "npc"), rot = 90, vjust = 0,
                                   gp = grid::gpar(fontsize = 8)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf))

p <- ggplot() +
  geom_point(data = bg, aes(umap_1, umap_2), colour = "grey85", size = 0.4, stroke = 0) +
  geom_point(data = fg, aes(umap_1, umap_2, colour = Alt_freq), size = 0.4, stroke = 0) +
  scale_colour_viridis_c(option = "plasma", direction = -1, limits = c(THR, 1), breaks = seq(THR, 1, 0.2)) +
  labs(colour = "Alteration frequency of H3 mutation") +
  guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5, direction = "horizontal",
                                  theme = theme(legend.key.width = unit(1.6, "in"), legend.key.height = unit(0.1, "in")))) +
  arrows +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "top", legend.title = element_text(size = 8), legend.text = element_text(size = 8),
        plot.margin = margin(4, 4, 14, 14))

save_figure(p, "umap_k27m_cells", width = 4, height = 3)
