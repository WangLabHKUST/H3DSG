# Volcano plot of the HALLMARK_E2F_TARGETS genes, H3-DSG-PL against H3-DSG-TL (bulk RNA-seq, 36 discovery tumours).
#
# The original plotting script of this panel was not kept. This script was rebuilt from the submitted panel and
# reproduces it: 189 points, 162 red and 27 grey (see README.md).
#
# Axes
#   x   log2(mean TPM in H3-DSG-PL - mean TPM in H3-DSG-TL): the log2 of the DIFFERENCE between the two group means.
#       It is not a fold change (ratio). A gene whose mean is not higher in H3-DSG-PL has no x value and is not drawn.
#   y   -log10 P of a two-sided Student t test (equal variances) on TPM.
#   Red: x > 1 and P < 0.05. Dashed lines: x = -1, x = 1 and P = 0.05.
#
# Input
#   ../Cell_of_Origin/data/Raw_SP_AB_RNA_count_data_with_gene_length.csv   gene-level read counts (hg19), gene length
#   HALLMARK_E2F_TARGETS_genes.tsv    the 200 genes of the gene set (MSigDB 2025.1.Hs)
#   sample_groups_discovery36.tsv     subtype of the 36 tumours
# Output
#   E2F_targets_volcano_table.tsv     one row per gene of the gene set
#   E2F_targets_volcano.pdf
#
# Run from this folder: Rscript E2F_targets_volcano_plot.R

library(ggplot2)
library(ggrepel)

p_thre <- 0.05
x_thre <- 1
# genes labelled on the submitted panel
show <- c("ORC6", "GINS1", "PLK1", "LMNB1", "SPC25", "SMC1A", "UBE2T", "STMN1", "ESPL1", "RAN", "EZH2", "MSH2",
          "TOP2A", "PLK4", "MKI67", "CDC20", "BRCA1", "CHEK1", "CDK4", "CHEK2", "AURKA", "MYC", "CDKN2A")

# TPM from the read counts: counts / gene length (kb), scaled to 1e6 per tumour over all genes of the table
d0 <- read.csv("../Cell_of_Origin/data/Raw_SP_AB_RNA_count_data_with_gene_length.csv", strip.white = TRUE,
               check.names = FALSE)
counts <- as.matrix(d0[, -(1:2)])
rownames(counts) <- d0$Geneid
rpk <- counts / (d0$Genelength / 1000)
tpm <- t(t(rpk) / colSums(rpk)) * 1e6

grp <- read.delim("sample_groups_discovery36.tsv", stringsAsFactors = FALSE)
stopifnot(nrow(grp) == 36, !anyDuplicated(grp$Patient_ID), all(grp$Patient_ID %in% colnames(tpm)),
          all(grp$Group %in% c("H3-DSG-PL", "H3-DSG-TL")))
pl <- grp$Patient_ID[grp$Group == "H3-DSG-PL"]
tl <- grp$Patient_ID[grp$Group == "H3-DSG-TL"]

# the gene set, matched to the count table by the current gene symbol
gs <- read.delim("HALLMARK_E2F_TARGETS_genes.tsv", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(gs$gene_symbol), !anyDuplicated(d0$Geneid[d0$Geneid %in% gs$gene_symbol]))
stat <- function(g) {
  if (!g %in% rownames(tpm)) return(c(NA, NA, NA))
  a <- tpm[g, pl]
  b <- tpm[g, tl]
  c(mean(a), mean(b), t.test(a, b, var.equal = TRUE)$p.value)
}
s <- t(sapply(gs$gene_symbol, stat))
info <- data.frame(gene = gs$gene_symbol, mean_TPM_PL = s[, 1], mean_TPM_TL = s[, 2], stringsAsFactors = FALSE)
info$difference_of_means <- info$mean_TPM_PL - info$mean_TPM_TL
info$log2_difference <- NA_real_
pos <- which(info$difference_of_means > 0)
info$log2_difference[pos] <- log2(info$difference_of_means[pos])
info$t_test_P <- s[, 3]
info$plotted <- !is.na(info$log2_difference)
info$lab <- ifelse(!info$plotted, NA, ifelse(info$log2_difference > x_thre & info$t_test_P < p_thre, "up", "ns"))
info$show <- ifelse(info$plotted & info$gene %in% show, "yes", "no")
info$note <- ""
info$note[is.na(info$mean_TPM_PL)] <- "not in the count table under this symbol; not drawn"
info$note[!is.na(info$mean_TPM_PL) & !info$plotted] <- "mean TPM not higher in H3-DSG-PL, log2 of the difference undefined; not drawn"
write.table(info, "E2F_targets_volcano_table.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("%d genes in the set, %d in the count table, %d plotted (%d red, %d grey)\n", nrow(info),
            sum(!is.na(info$mean_TPM_PL)), sum(info$plotted), sum(info$lab == "up", na.rm = TRUE),
            sum(info$lab == "ns", na.rm = TRUE)))

d <- info[info$plotted, ]
p1 <- ggplot(d, aes(log2_difference, -log10(t_test_P))) +
  geom_point(aes(color = lab), size = 1.8) +
  geom_hline(yintercept = -log10(p_thre), linetype = "dashed") +
  geom_vline(xintercept = c(-x_thre, x_thre), linetype = "dashed") +
  geom_text_repel(aes(label = ifelse(show == "yes", gene, ""), color = lab), max.overlaps = 1000, size = 3,
                  box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"),
                  segment.color = "black", show.legend = FALSE, seed = 1) +
  scale_color_manual(values = c(ns = "grey", up = "#FF3300")) +
  scale_x_continuous(limits = c(-4.5, 8.5), breaks = seq(-4, 8, 2)) +
  labs(x = "Log2 (difference in mean TPM, H3-DSG-PL - H3-DSG-TL)", y = "-log10(P value)") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")
pdf("E2F_targets_volcano.pdf", width = 5, height = 6)
print(p1)
dev.off()
