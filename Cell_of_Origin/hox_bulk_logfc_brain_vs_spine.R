# Fig. 2a: log2 fold change of gene expression in bulk RNA-seq, spinal cord tumours against brain DMGs.
#   x axis: H3-DSG-PL (pons-like) against pons DMG
#   y axis: H3-DSG-TL (thalamus-like) against thalamus DMG
# One dot per protein-coding gene; HOX genes in red and labelled, all other genes as a density-coloured cloud.
#
# Input (data/):
#   Spinalcord.H3K27M.expression.RPKM.txt     RPKM of 40 spinal cord tumours (columns are SP IDs), hg19
#   H3K27altered_project_metadata_updated.tsv project metadata; TSNE.cluster A = pons-like (PL), B = thalamus-like (TL)
#   CGGA.H3K27M.hg38.expression.RPKM.txt      RPKM of H3 K27M gliomas of the CGGA cohort (public, hg38); five thalamic ones are used
#   StJude.DIPG.RPKM.txt                      RPKM of five pons DMGs of the St Jude PeCan cohort (public)
#   Homo_sapiens.GRCh37.75.chr.genes.txt      gene biotypes (Ensembl GRCh37.75); only protein-coding genes are kept
# Output: hox_bulk_logfc_brain_vs_spine.pdf
#
# Grouping: this script groups the tumours by the project metadata (TSNE.cluster A / B), as in the submitted figure.
# The Source Data of the figure group them by Table S1, which moves SP02 and SP42 from TL to PL.
# Run from the Cell_of_Origin/ directory:  Rscript hox_bulk_logfc_brain_vs_spine.R

library(ggplot2)
library(ggrepel)

# ---- data ----
spine <- read.delim("data/Spinalcord.H3K27M.expression.RPKM.txt")
meta  <- read.delim("data/H3K27altered_project_metadata_updated.tsv")
spine_pl <- spine[, names(spine) %in% meta$ID2[meta$TSNE.cluster == "A"]]   # pons-like
spine_tl <- spine[, names(spine) %in% meta$ID2[meta$TSNE.cluster == "B"]]   # thalamus-like
stopifnot(ncol(spine_pl) == 16, ncol(spine_tl) == 24)                        # all 40 tumours matched

cgga <- read.delim("data/CGGA.H3K27M.hg38.expression.RPKM.txt", row.names = 1)
thalamus_dmg <- cgga[, names(cgga) %in% c("CGGA_791", "CGGA_P266", "CGGA_1185", "CGGA_1320", "CGGA_1735")]
stopifnot(ncol(thalamus_dmg) == 5)

pons_all <- read.delim("data/StJude.DIPG.RPKM.txt")
pons_all <- pons_all[!is.na(pons_all$gene), ]
pons_all <- pons_all[!duplicated(pons_all$gene), ]                           # first row of a duplicated gene name
stopifnot(all(c("SJHGG002", "SJHGG003", "SJHGG004", "SJHGG102", "SJHGG145") %in% names(pons_all)))
pons_dmg <- as.matrix(pons_all[, 1:5]); rownames(pons_dmg) <- pons_all$gene

biotype <- read.delim("data/Homo_sapiens.GRCh37.75.chr.genes.txt", header = FALSE)
protein_coding <- biotype$V7[biotype$V2 == "protein_coding"]

# ---- log2 fold change of the group medians (RPKM) against brain DMG ----
# lfc = log2(median RPKM in spinal tumours + 1) - log2(median RPKM in brain DMGs + 1)
lfc_vs_brain <- function(spine_grp, brain, genes) {
  genes <- genes[genes %in% protein_coding]
  med_spine <- apply(as.matrix(spine_grp[genes, ]), 1, median)
  med_brain <- apply(as.matrix(brain[genes, ]), 1, median)
  setNames(log2(med_spine + 1) - log2(med_brain + 1), genes)
}
lfc_pons     <- lfc_vs_brain(spine_pl, pons_dmg,     unique(intersect(rownames(spine_pl), rownames(pons_dmg))))
lfc_thalamus <- lfc_vs_brain(spine_tl, thalamus_dmg, unique(intersect(rownames(spine_pl), rownames(thalamus_dmg))))

mg <- data.frame(gene = intersect(names(lfc_pons), names(lfc_thalamus)))
mg$logFC_pons     <- lfc_pons[mg$gene]
mg$logFC_thalamus <- lfc_thalamus[mg$gene]
mg$lab <- ifelse(startsWith(mg$gene, "HOX"), mg$gene, NA)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  dens$z[cbind(ix, iy)]
}
mg$density <- get_density(mg$logFC_pons, mg$logFC_thalamus, n = 100)

# ---- plot ----
# Axis titles were added to the submitted panel afterwards:
#   x: "Log2 (H3-DSG-PL / pons DMG)", y: "Log2 (H3-DSG-TL / thalamus DMG)"
p <- ggplot() +
  geom_point(data = mg[!(mg$logFC_thalamus == 0 & mg$logFC_pons < -5), ], show.legend = FALSE,
             mapping = aes(color = density, x = logFC_pons, y = logFC_thalamus), alpha = 0.5, size = 0.1) +
  geom_point(mg[!is.na(mg$lab), ], mapping = aes(x = logFC_pons, y = logFC_thalamus), size = 0.5, col = "red") +
  geom_vline(xintercept = 0, lty = 2, col = "#dddddd") +
  geom_hline(yintercept = 0, lty = 2, col = "#dddddd") +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = "#dddddd") +
  scale_x_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  scale_color_viridis_c() +
  geom_text_repel(data = mg[(mg$logFC_pons > 1 & mg$logFC_pons < 3.5) | (mg$logFC_thalamus > 1 & mg$logFC_thalamus < 3.5), ],
                  mapping = aes(x = logFC_pons, y = logFC_thalamus, label = lab),
                  bg.color = "white", bg.r = 0.25, col = "red", min.segment.length = 0,
                  segment.size = 0.25, size = 1.5, max.overlaps = 20, seed = 1) +
  theme_classic() + theme(axis.text = element_text(color = "black")) +
  labs(x = " ", y = " ")
ggsave("hox_bulk_logfc_brain_vs_spine.pdf", p, width = 154 / 72, height = 129 / 72)   # page size of the original figure
