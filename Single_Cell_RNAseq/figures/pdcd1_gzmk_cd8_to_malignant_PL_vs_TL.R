# PDCD1+ GZMK+ CD8 T cells per malignant cell, PL vs TL.
#
# Reads  ../data/panel_data/tnk_cells.tsv.gz and tnk_samples.tsv
#        (written by ../pipeline/99b_export_tnk_panels.R)
# Writes pdcd1_gzmk_cd8_to_malignant_PL_vs_TL.pdf and .png next to this script.
#
# Numerator: cells of the T/NK state GZMK+ CD8 T whose log-normalised PDCD1 value (RNA assay) is
# above PDCD1_CUT. Denominator: malignant cells of the sample. The ratio is computed per sample
# and expressed in units of 10^-2 (a percentage); a patient with several samples (TL_6 to TL_10)
# is the mean of its samples, so the comparison has one point per patient.
# PL_9 has no T/NK cells in the object (16 T/NK cells, below the 30-cell minimum): its ratio is
# counted as 0, so PL has 8 patients and TL 7. Wilcoxon rank-sum test, two-sided, normal
# approximation.

suppressPackageStartupMessages(library(ggplot2))

GROUP_COLOURS <- c(PL = "#990000", TL = "#2b8cbe")
PDCD1_CUT <- 1      # log-normalised value; expm1(1) = 1.7 normalised counts per 10,000

cells   <- read.delim(gzfile("../data/panel_data/tnk_cells.tsv.gz"), colClasses = c(cell_id = "character"))
samples <- read.delim("../data/panel_data/tnk_samples.tsv", colClasses = c(n_tnk = "integer", n_malignant = "integer"))

# per sample, all 19
sample_df <- samples
pos <- cells$tnk_state == "GZMK+ CD8 T" & cells$PDCD1_lognorm > PDCD1_CUT
sample_df$n_pdcd1_gzmk_cd8 <- as.vector(table(factor(cells$sid[pos], levels = samples$sid)))
sample_df$ratio_e2 <- 100 * sample_df$n_pdcd1_gzmk_cd8 / sample_df$n_malignant
sample_df$ratio_e2[sample_df$n_tnk == 0] <- 0

# the plotted data: one row per patient
box_df <- aggregate(ratio_e2 ~ group + patient_id, data = sample_df, FUN = mean)
box_df$group <- factor(box_df$group, levels = c("PL", "TL"))
test <- wilcox.test(ratio_e2 ~ group, data = box_df, exact = FALSE)

set.seed(1)
p <- ggplot(box_df, aes(group, ratio_e2, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.2, size = 4,
           label = sprintf("italic(P)==%s", shQuote(sprintf("%.1e", test$p.value))), parse = TRUE) +
  scale_fill_manual(values = GROUP_COLOURS) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = expression("PDCD1+ GZMK+ CD8 T / malignant cells (" * 10^-2 * ")")) +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(colour = "black"))

ggsave("pdcd1_gzmk_cd8_to_malignant_PL_vs_TL.pdf", p, width = 2.6, height = 4)
ggsave("pdcd1_gzmk_cd8_to_malignant_PL_vs_TL.png", p, width = 2.6, height = 4, dpi = 300, bg = "white")
