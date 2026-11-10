# Fraction of GZMK+ CD8 T cells among the T cells, PL vs TL.
#
# Reads  ../data/panel_data/tnk_cells.tsv.gz   (written by ../pipeline/99b_export_tnk_panels.R)
# Writes gzmk_cd8_fraction_PL_vs_TL.pdf and gzmk_cd8_fraction_PL_vs_TL.png next to this script.
#
# GZMK+ CD8 T is the T/NK state of that name (see umap_tnk_states.R), not a marker threshold.
# T cells are all T/NK cells of the object except the two NK states. The fraction is computed per
# sample; a patient with several samples (TL_6 to TL_10) is the mean of its samples, so the
# comparison has one point per patient. A sample with no T/NK cells in the object (PL_9) has no
# fraction and its patient is left out. Wilcoxon rank-sum test, two-sided, normal approximation.

suppressPackageStartupMessages(library(ggplot2))

GROUP_COLOURS <- c(PL = "#990000", TL = "#2b8cbe")

cells <- read.delim(gzfile("../data/panel_data/tnk_cells.tsv.gz"), colClasses = c(cell_id = "character"))
cells$is_T <- !grepl("NK", cells$tnk_state)

# per sample (only samples with T/NK cells are in the table)
sample_df <- do.call(rbind, lapply(split(cells, cells$sid), function(d)
  data.frame(renamed_id = d$renamed_id[1], group = d$group[1], patient_id = d$patient_id[1],
             n_T = sum(d$is_T), n_gzmk_cd8 = sum(d$tnk_state == "GZMK+ CD8 T"))))
sample_df$pct <- 100 * sample_df$n_gzmk_cd8 / sample_df$n_T

# the plotted data: one row per patient
box_df <- aggregate(pct ~ group + patient_id, data = sample_df, FUN = mean)
box_df$group <- factor(box_df$group, levels = c("PL", "TL"))
test <- wilcox.test(pct ~ group, data = box_df, exact = FALSE)

set.seed(1)
p <- ggplot(box_df, aes(group, pct, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7) +
  geom_jitter(width = 0.1, height = 0, size = 2) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.2, size = 4,
           label = sprintf("italic(P)==%s", shQuote(sprintf("%.1e", test$p.value))), parse = TRUE) +
  scale_fill_manual(values = GROUP_COLOURS) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = NULL, y = "GZMK+ CD8 T in T cells (%)") +
  theme_classic() +
  theme(legend.position = "none", axis.text = element_text(colour = "black"))

ggsave("gzmk_cd8_fraction_PL_vs_TL.pdf", p, width = 2.6, height = 4)
ggsave("gzmk_cd8_fraction_PL_vs_TL.png", p, width = 2.6, height = 4, dpi = 300, bg = "white")
