# =============================================================================
# Similarity of the NMF programs of all samples, grouped into metaprograms
#
# Input : ../data/metaprograms/program_similarity_matrix.tsv
#         ../data/metaprograms/program_metaprogram_assignment.tsv
# Output: program_similarity_matrix.pdf / .png in this directory
#
# Similarity of two programs = number of genes their top 100 genes share (0 to 100). It is not
# a correlation and not a Jaccard index. Rows and columns are the 84 programs that belong to
# one of the six metaprograms (mpg1 .. mpg6; the 23 programs left out of every metaprogram
# are not drawn), ordered by metaprogram and, within a metaprogram, by hierarchical
# clustering of the 107 programs on this similarity. Annotation bars: sample of the program,
# metaprogram, and the Pearson correlation between the usage of the program in the cells of its
# sample and the S and the G2M phase score of the cells.
# =============================================================================

suppressPackageStartupMessages({ library(pheatmap); library(viridis) })

# white-to-dark colour scale used for the similarity (0 = white, 100 = dark purple)
custom_magma <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10),
                  rev(viridis::magma(323, begin = 0.18)))

sim <- as.matrix(read.delim("../data/metaprograms/program_similarity_matrix.tsv", row.names = 1, check.names = FALSE))
asg <- read.delim("../data/metaprograms/program_metaprogram_assignment.tsv", check.names = FALSE)
stopifnot(identical(rownames(sim), colnames(sim)), setequal(asg$program, rownames(sim)))

drawn <- asg[!is.na(asg$block_order), ]
drawn <- drawn[order(drawn$block_order), ]                      # metaprogram blocks, left to right
stopifnot(nrow(drawn) == 84)
sim_drawn <- sim[drawn$program, drawn$program]

ann <- data.frame(Sample = drawn$sample, Metaprogram = sub("^M", "mpg", drawn$metaprogram),
                  S.Score = drawn$cor_S_score, G2M.Score = drawn$cor_G2M_score, row.names = drawn$program)
# one colour per sample, in the order the samples were processed (hue palette)
sample_order <- unique(asg$sample[order(asg$extraction_order)])
ann_colors <- list(Sample = setNames(scales::hue_pal()(length(sample_order)), sample_order))
gaps <- cumsum(table(factor(ann$Metaprogram, levels = paste0("mpg", 1:6))))

ph <- pheatmap(sim_drawn, cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps, gaps_col = gaps,
               annotation_col = ann, annotation_colors = ann_colors, color = custom_magma,
               show_colnames = FALSE, show_rownames = FALSE, border_color = NA, silent = TRUE)

save_figure <- function(ph, stem, width, height, dpi = 300) {
  pdf(paste0(stem, ".pdf"), width = width, height = height); grid::grid.draw(ph$gtable); invisible(dev.off())
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi, bg = "white")
  grid::grid.draw(ph$gtable); invisible(dev.off())
}
save_figure(ph, "program_similarity_matrix", width = 9, height = 9.5)
