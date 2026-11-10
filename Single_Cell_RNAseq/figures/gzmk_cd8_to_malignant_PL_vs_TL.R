# GZMK+ CD8 T cells per 1,000 malignant cells, PL vs TL.
#
# Reads  ../data/panel_data/tnk_cells.tsv.gz and tnk_samples.tsv
#        (written by ../pipeline/99b_export_tnk_panels.R)
# Writes gzmk_cd8_to_malignant_PL_vs_TL.pdf and gzmk_cd8_to_malignant_PL_vs_TL.png next to this script.
#
# GZMK+ CD8 T is the T/NK state of that name (see umap_tnk_states.R), not a marker threshold.
# The counts of the samples of a patient are added before dividing (TL_6 to TL_10 are one patient),
# so the comparison has one point per patient. PL_9 has no T/NK cells in the object (16 T/NK cells,
# below the 30-cell minimum): it enters with 0 GZMK+ CD8 T cells, so PL has 8 patients and TL 7.
# Wilcoxon rank-sum test, two-sided, normal approximation.

suppressPackageStartupMessages(library(ggplot2))

GROUP_COLOURS <- c(PL = "#990000", TL = "#2b8cbe")

cells   <- read.delim(gzfile("../data/panel_data/tnk_cells.tsv.gz"), colClasses = c(cell_id = "character"))
samples <- read.delim("../data/panel_data/tnk_samples.tsv", colClasses = c(n_tnk = "integer", n_malignant = "integer"))

# per sample, all 19
sample_df <- samples
sample_df$n_gzmk_cd8 <- as.vector(table(factor(cells$sid[cells$tnk_state == "GZMK+ CD8 T"], levels = samples$sid)))

# the plotted data: one row per patient
box_df <- aggregate(cbind(n_gzmk_cd8, n_malignant) ~ group + patient_id, data = sample_df, FUN = sum)
box_df$per_1000 <- 1000 * box_df$n_gzmk_cd8 / box_df$n_malignant
box_df$group <- factor(box_df$group, levels = c("PL", "TL"))
test <- wilcox.test(per_1000 ~ group, data = box_df, exact = FALSE)

set.seed(1)
p <- ggplot(box_df, aes(group, per_1000, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.2, size = 4,
           label = sprintf("italic(P)==%s", shQuote(sprintf("%.1e", test$p.value))), parse = TRUE) +
  scale_fill_manual(values = GROUP_COLOURS) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "GZMK+ CD8 T per 1,000 malignant cells") +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(colour = "black"))

ggsave("gzmk_cd8_to_malignant_PL_vs_TL.pdf", p, width = 2.6, height = 4)
ggsave("gzmk_cd8_to_malignant_PL_vs_TL.png", p, width = 2.6, height = 4, dpi = 300, bg = "white")
