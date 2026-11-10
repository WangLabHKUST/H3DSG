# Extended Data Fig. 2g: chromatin accessibility (bulk ATAC-seq) and H3K27me3 (bulk CUT&Tag) at the promoters of 39 HOX genes,
# H3-DSG-TL (n = 2) and H3-DSG-PL (n = 5). Each strip shows the mean over the tumours of a group; darker = higher.
#
# Input: data/hox_atac_k27me3_promoter_signal.tsv   one row per gene, assay and tumour (columns gene, value, HOX, int, type, ID);
#        ID is the SP ID of the tumour, type is ATAC or H3K27Me3, value the mean signal in the HOX promoter.
# Output: hox_atac_k27me3_strips.pdf
#
# How the table was made: the ATAC-seq and CUT&Tag reads of each tumour were converted to CPM-normalised coverage (bigWig, then
# bedGraph; these files are not deposited). For each of the 39 HOX genes (GRCh37.75 gene coordinates) the signal was averaged over
# the promoter, +/- 2500 bp around the transcription start site, giving one value per gene, assay and tumour.
# Eight tumours were profiled; the eighth, an H3 wild-type tumour used as control, is not part of the figure and is not in the table.
#
# Groups: H3-DSG-TL = SP32 and SP36; H3-DSG-PL = SP02, SP03, SP05, SP09 and SP10. These are the groups of Table S1 and of the
# sample sheet of the sequencing data (group A = PL, group B = TL). Both are set here by SP ID.
# Run from the Cell_of_Origin/ directory:  Rscript hox_atac_k27me3_strips.R

library(ggplot2)
library(cowplot)

df <- read.delim("data/hox_atac_k27me3_promoter_signal.tsv")
thalamus_like <- c("SP32", "SP36")
stopifnot(length(unique(df$ID)) == 7, all(thalamus_like %in% df$ID), nrow(df) == 7 * 2 * 39)
df$tumorgroup <- ifelse(df$ID %in% thalamus_like, "TL", "PL")
df$gene <- factor(df$gene, levels = unique(df$gene[order(df$HOX, df$int)]))        # HOXA1 ... HOXD13

# group means per gene and assay
dfacmean     <- aggregate(value ~ gene + tumorgroup, FUN = mean, data = df[df$type == "ATAC", ])
dfk27me3mean <- aggregate(value ~ gene + tumorgroup, FUN = mean, data = df[df$type == "H3K27Me3", ])

# ---- one strip: one HOX cluster, one group, one assay ----
# midpoint of the colour scale is chosen per assay and cluster (as in the original figure); the scale spans the values of the strip
strip <- function(d, cluster, grp, midpoint) {
  ggplot(d[grepl(cluster, d$gene) & d$tumorgroup == grp, ], aes(x = gene, y = 1)) +
    geom_tile(aes(fill = value), color = "black", show.legend = F) +
    scale_fill_gradient2(low = "#f7fcfd", mid = "#8c96c6", high = "#4d004b", midpoint = midpoint) +
    theme_bw() + theme(axis.ticks = element_blank(), axis.text = element_blank(), rect = element_blank(),
                       panel.grid = element_blank(), axis.title = element_blank())
}
clusters <- c("HOXA", "HOXB", "HOXC", "HOXD")
mid_k27me3 <- c(HOXA = 0.3, HOXB = 0.3, HOXC = 3, HOXD = 0.3)
mid_atac   <- c(HOXA = 0.1, HOXB = 0.1, HOXC = 0.1, HOXD = 0.1)
block <- function(grp) c(lapply(clusters, function(h) strip(dfk27me3mean, h, grp, mid_k27me3[h])),
                         rep(list(NULL), 4),                                       # negative-height row below: pulls the two strips together
                         lapply(clusters, function(h) strip(dfacmean, h, grp, mid_atac[h])))

# gene names under the last strips
gene_labels <- lapply(clusters, function(h)
  ggplot(data.frame(gene = factor(levels(df$gene)[grepl(h, levels(df$gene))], levels = levels(df$gene)), y = 1), aes(gene, y)) +
    geom_blank() + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "grey40"), axis.text.y = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), rect = element_blank()))

plots <- c(block("TL"), rep(list(NULL), 4), block("PL"), gene_labels)
# rows: TL H3K27me3, spacer, TL ATAC | gap | PL H3K27me3, spacer, PL ATAC | gene names
g <- plot_grid(plotlist = plots, ncol = 4,
               rel_heights = c(0.1, -0.05, 0.1, 0.06, 0.1, -0.05, 0.1, 0.18),
               rel_widths = rep(1 / 4, 4))
# Group names and assay names were added to the submitted panel afterwards (left: H3-DSG-TL above H3-DSG-PL;
# right: H3K27me3 above ATAC-seq within each group).
ggsave("hox_atac_k27me3_strips.pdf", g, width = 416 / 72, height = 146 / 72)       # strip size of the original figure
