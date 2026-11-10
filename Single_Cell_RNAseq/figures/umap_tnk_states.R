# UMAP of the T/NK cells, coloured by state.
#
# Reads  ../data/panel_data/tnk_cells.tsv.gz   (written by ../pipeline/99b_export_tnk_panels.R)
# Writes umap_tnk_states.pdf and umap_tnk_states.png next to this script.
#
# The embedding is the one computed on the 4,201 T/NK cells alone (30 PCs of the CCA-integrated
# T/NK object, UMAP on PCs 1-10), not the UMAP of all cells.

suppressPackageStartupMessages({ library(ggplot2); library(ggrepel) })

# ---- T/NK states: object name, name shown in figures, colour (RColorBrewer "Paired") ----------
TNK_STATES  <- c("GZMK+ CD8 T", "CD8 Tex", "CD8 T_m", "CD4 T_m", "CD56-CD16+ NK", "Stress T",
                 "CD56+CD16- NK", "Prolif", "Treg", "Naive T", "IFN")
TNK_LABELS  <- c("GZMK+ CD8 T", "exhausted CD8 T", "CD8 memory T", "CD4 memory T", "CD56-CD16+ NK", "Stress T",
                 "CD56+CD16- NK", "proliferative T", "regulatory T", "Naive T", "Interferon T")
TNK_COLOURS <- setNames(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                          "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"), TNK_LABELS)

cells <- read.delim(gzfile("../data/panel_data/tnk_cells.tsv.gz"), colClasses = c(cell_id = "character"))
stopifnot(all(cells$tnk_state %in% TNK_STATES))

# the plotted data: one row per T/NK cell
umap_df <- data.frame(cell_id = cells$cell_id, renamed_id = cells$renamed_id, group = cells$group,
                      tnk_state = cells$tnk_state,
                      state = factor(TNK_LABELS[match(cells$tnk_state, TNK_STATES)], levels = TNK_LABELS),
                      umap_1 = cells$umap_1, umap_2 = cells$umap_2)

# state names sit at the median of their cells, as DimPlot(label = TRUE) places them
centres <- aggregate(cbind(umap_1, umap_2) ~ state, data = umap_df, FUN = median)

# axis arrows in the lower-left corner instead of axes
x0 <- min(umap_df$umap_1); y0 <- min(umap_df$umap_2)
len_x <- diff(range(umap_df$umap_1)) * 0.12; len_y <- diff(range(umap_df$umap_2)) * 0.12

p <- ggplot(umap_df, aes(umap_1, umap_2, colour = state)) +
  geom_point(size = 0.4, shape = 16) +
  geom_text_repel(data = centres, aes(label = state), colour = "black", size = 3.5, seed = 1,
                  max.overlaps = Inf, box.padding = 0.2, segment.color = NA) +
  annotate("segment", x = x0, y = y0, xend = x0 + len_x, yend = y0, arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + len_y, arrow = arrow(length = unit(2, "mm"), type = "closed")) +
  annotate("text", x = x0 + len_x / 2, y = y0 - len_y * 0.18, label = "UMAP1", size = 3.2) +
  annotate("text", x = x0 - len_x * 0.14, y = y0 + len_y / 2, label = "UMAP2", size = 3.2, angle = 90) +
  scale_colour_manual(values = TNK_COLOURS) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(5, 5, 15, 15))

ggsave("umap_tnk_states.pdf", p, width = 6, height = 5)
ggsave("umap_tnk_states.png", p, width = 6, height = 5, dpi = 300, bg = "white")
