# =============================================================================
# Overlap of the six metaprograms with the cell states of Liu et al. 2022 (PMID 36471067)
#
# Input : ../data/metaprograms/metaprogram_overlap_liu2022.tsv
# Output: metaprogram_overlap_liu2022.pdf / .png in this directory
#
# Colour = number of genes shared by the 50 signature genes of a metaprogram of this study
# (rows, mpg1 .. mpg6) and the gene set of a state of Liu et al. (columns). The three OPC-like
# gene sets of Liu et al. are merged into one. Rows are ordered so that the best match of each
# metaprogram lies on the diagonal. The colour scale is that of pheatmap (blue = 0, red = highest).
# =============================================================================

suppressPackageStartupMessages(library(pheatmap))

ovl <- read.delim("../data/metaprograms/metaprogram_overlap_liu2022.tsv")
rows <- unique(ovl[order(ovl$row_order), c("metaprogram", "row_order")])
cols <- unique(ovl[order(ovl$col_order), c("liu_state", "col_order")])
mat <- matrix(NA_integer_, nrow(rows), nrow(cols),
              dimnames = list(sub("^M", "mpg", rows$metaprogram), cols$liu_state))
mat[cbind(match(ovl$metaprogram, rows$metaprogram), match(ovl$liu_state, cols$liu_state))] <- ovl$shared_genes
stopifnot(!anyNA(mat), nrow(mat) == 6, ncol(mat) == 6)

ph <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, border_color = "white",
               display_numbers = FALSE, angle_col = 45, cellwidth = 36, cellheight = 36,
               main = "Metaprograms in Liu et al", silent = TRUE)

draw_panel <- function() {
  grid::grid.newpage(); grid::grid.draw(ph$gtable)
  grid::grid.text("Metaprograms defined in this study", x = grid::unit(0.015, "npc"), rot = 90, gp = grid::gpar(fontsize = 10))
}
pdf("metaprogram_overlap_liu2022.pdf", width = 4.6, height = 4.2); draw_panel(); invisible(dev.off())
png("metaprogram_overlap_liu2022.png", width = 4.6, height = 4.2, units = "in", res = 300, bg = "white")
draw_panel(); invisible(dev.off())
