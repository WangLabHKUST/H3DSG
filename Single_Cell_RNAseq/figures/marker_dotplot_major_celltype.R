# =============================================================================
# Marker genes of the major cell types
#
# Input : ../data/panel_data/marker_gene_expression_nonzero.tsv.gz   (pipeline/99b_export_major_panels.R)
#         ../data/cell_level/cell_annotation.tsv.gz                  cell_type
# Output: marker_dotplot_major_celltype.pdf / .png in this directory
#
# The dot plot is computed here the way Seurat's DotPlot(scale = TRUE) does, from the
# per-cell expression of the marker genes (log-normalised RNA values; a gene absent from a
# cell is 0) over all 132,343 cells:
#   avg.exp         mean over the cells of the cell type of expm1(expression)
#   pct.exp         percentage of the cells of the cell type with expression > 0
#   avg.exp.scaled  log1p(avg.exp), z-scored across the seven cell types for each gene and
#                   clipped to [-2.5, 2.5] (colour)
# Colour = avg.exp.scaled, size = pct.exp (radius proportional, 0 to 100).
# =============================================================================

source("panel_helpers.R")

MARKER_GENES <- c("PTPRZ1", "SOX2", "PTPRC", "AIF1", "CD14", "FCGR3B", "S100A9", "CSF3R", "CD3D",
                  "IL32", "NKG7", "CD79A", "MS4A1", "COL1A1", "VWF", "PLP1", "MAG")   # top to bottom

ann <- read.delim(gzfile("../data/cell_level/cell_annotation.tsv.gz"),
                  colClasses = c(patient_id = "character"))[, c("cell_id", "renamed_id", "cell_type")]
ex  <- read.delim(gzfile("../data/panel_data/marker_gene_expression_nonzero.tsv.gz"),
                  colClasses = c("character", "character", "numeric"))
# the expression table is keyed by renamed_id + the 16-bp barcode (see pipeline/99b_export_major_panels.R)
ann$cell_key <- paste0(ann$renamed_id, "_", sub("^.*_([ACGT]{16})(-1)?$", "\\1", ann$cell_id))
stopifnot(setequal(ex$gene, MARKER_GENES), !anyDuplicated(ann$cell_key), all(ex$cell_key %in% ann$cell_key))

ann$celltype <- celltype_label(ann$cell_type)
n_cells <- table(ann$celltype)

ex$celltype <- ann$celltype[match(ex$cell_key, ann$cell_key)]
key <- paste(ex$gene, ex$celltype, sep = "\r")
sum_expm1 <- tapply(expm1(ex$expr), key, sum)
n_pos     <- tapply(ex$expr > 0, key, sum)

dot_df <- expand.grid(gene = MARKER_GENES, celltype = CELLTYPE_ORDER, stringsAsFactors = FALSE)
k <- paste(dot_df$gene, dot_df$celltype, sep = "\r")
dot_df$n_cells <- as.integer(n_cells[dot_df$celltype])
dot_df$avg.exp <- ifelse(is.na(sum_expm1[k]), 0, sum_expm1[k]) / dot_df$n_cells
dot_df$pct.exp <- 100 * ifelse(is.na(n_pos[k]), 0, n_pos[k]) / dot_df$n_cells
dot_df$avg.exp.scaled <- ave(log1p(dot_df$avg.exp), dot_df$gene, FUN = function(v) pmin(pmax(as.vector(scale(v)), -2.5), 2.5))

plot_df <- dot_df
plot_df$gene     <- factor(plot_df$gene, levels = rev(MARKER_GENES))
plot_df$celltype <- factor(plot_df$celltype, levels = CELLTYPE_ORDER)

p <- ggplot(plot_df, aes(celltype, gene)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled)) +
  scale_colour_gradientn(colours = c("lightgrey", "blue"), limits = range(plot_df$avg.exp.scaled)) +
  scale_radius(range = c(0, 5), limits = c(0, 100)) +
  labs(x = NULL, y = NULL, colour = "Avg.exp", size = "Pct.exp") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_figure(p, "marker_dotplot_major_celltype", width = 5, height = max(4, 0.25 * length(MARKER_GENES)))
