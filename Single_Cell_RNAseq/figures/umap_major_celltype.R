# =============================================================================
# UMAP of all 132,343 cells coloured by major cell type
#
# Input : ../data/cell_level/cell_annotation.tsv.gz   cell_type
#         ../data/cell_level/umap_coordinates.tsv.gz  umap_corr_1, umap_corr_2
# Output: umap_major_celltype.pdf / .png in this directory
#
# The embedding is `umap_pc20_correlation` (20 PCs, correlation metric), the columns
# umap_corr_1 / umap_corr_2; the panel is not drawn on umap_1 / umap_2. Cells are drawn in
# alphabetical order of the label (Seurat's DimPlot(order = TRUE)), so that the small
# populations lie on top of the malignant cells. Points are rasterised in a 512 x 512 grid
# with point size 3, as DimPlot(raster = TRUE, pt.size = 3) does.
# =============================================================================

source("panel_helpers.R")
suppressPackageStartupMessages(library(scattermore))

ann  <- read.delim(gzfile("../data/cell_level/cell_annotation.tsv.gz"),
                   colClasses = c(patient_id = "character"))[, c("cell_id", "sid", "cell_type")]
umap <- read.delim(gzfile("../data/cell_level/umap_coordinates.tsv.gz"))
stopifnot(identical(ann$cell_id, umap$cell_id))

umap_df <- data.frame(cell_id = ann$cell_id, umap_1 = umap$umap_corr_1, umap_2 = umap$umap_corr_2,
                      cell_type = celltype_label(ann$cell_type))
n_cells   <- nrow(umap_df)
n_samples <- length(unique(ann$sid))

# drawing order: alphabetical, last on top
DRAW_ORDER <- c("B", "Malignant", "Myeloid", "Neutrophil", "Oligo", "T/NK", "Vascular")
stopifnot(setequal(DRAW_ORDER, CELLTYPE_ORDER))
draw_df <- umap_df[order(match(as.character(umap_df$cell_type), DRAW_ORDER)), ]

# labels at the median of each cell type, as DimPlot(label = TRUE)
label_df <- do.call(rbind, lapply(split(umap_df, umap_df$cell_type), function(d)
  data.frame(cell_type = d$cell_type[1], x = median(d$umap_1), y = median(d$umap_2))))

# the labels of the small populations are moved off the points
NUDGE <- list(Oligo = c(1.2, -0.6), Vascular = c(2.2, -0.4), Neutrophil = c(1.6, -1.0),
              "T/NK" = c(0, 1.7), B = c(0.9, 0.9))
for (ct in names(NUDGE)) {
  i <- label_df$cell_type == ct
  label_df$x[i] <- label_df$x[i] + NUDGE[[ct]][1]; label_df$y[i] <- label_df$y[i] + NUDGE[[ct]][2]
}

# axis arrows in the empty lower left corner
ax0 <- c(min(umap_df$umap_1), min(umap_df$umap_2) + 1.2); AX_LEN <- 3.2

p <- ggplot(draw_df, aes(umap_1, umap_2, colour = cell_type)) +
  geom_scattermore(pointsize = 3, pixels = c(512, 512)) +
  scale_colour_manual(values = CELLTYPE_COLORS_UMAP) +
  geom_text(data = label_df, aes(x, y, label = cell_type), colour = "black", size = 4, inherit.aes = FALSE) +
  annotate("text", x = max(umap_df$umap_1), y = max(umap_df$umap_2), hjust = 1, vjust = 1, size = 4.2,
           label = sprintf("N = %d samples\nn = %s cells", n_samples, format(n_cells, big.mark = ","))) +
  annotate("segment", x = ax0[1], xend = ax0[1] + AX_LEN, y = ax0[2], yend = ax0[2],
           arrow = arrow(length = unit(0.1, "in"), type = "closed"), linewidth = 0.7) +
  annotate("segment", x = ax0[1], xend = ax0[1], y = ax0[2], yend = ax0[2] + AX_LEN,
           arrow = arrow(length = unit(0.1, "in"), type = "closed"), linewidth = 0.7) +
  annotate("text", x = ax0[1] + AX_LEN / 2, y = ax0[2] - 0.5, label = "UMAP1", size = 4, vjust = 1) +
  annotate("text", x = ax0[1] - 0.4, y = ax0[2] + AX_LEN / 2, label = "UMAP2", size = 4, angle = 90, vjust = 0) +
  coord_fixed() +
  theme_void() + theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5))

save_figure(p, "umap_major_celltype", width = 6.5, height = 5)
