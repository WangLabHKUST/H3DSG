# Dot plot of marker genes across the T/NK states.
#
# Reads  ../data/panel_data/tnk_marker_dotplot.tsv   (written by ../pipeline/99b_export_tnk_panels.R)
# Writes marker_dotplot_tnk_states.pdf and marker_dotplot_tnk_states.png next to this script.
#
# As Seurat's DotPlot draws it: dot size is the percentage of cells with a non-zero value, colour is
# the mean expression (mean of the de-logged RNA-assay values), as log1p z-scored across the 11
# states and clipped at +-2.5. Both are computed in 99b_export_tnk_panels.R.

suppressPackageStartupMessages(library(ggplot2))

TNK_STATES <- c("GZMK+ CD8 T", "CD8 Tex", "CD8 T_m", "CD4 T_m", "CD56-CD16+ NK", "Stress T",
                "CD56+CD16- NK", "Prolif", "Treg", "Naive T", "IFN")
TNK_LABELS <- c("GZMK+ CD8 T", "exhausted CD8 T", "CD8 memory T", "CD4 memory T", "CD56-CD16+ NK", "Stress T",
                "CD56+CD16- NK", "proliferative T", "regulatory T", "Naive T", "Interferon T")

dot_df <- read.delim("../data/panel_data/tnk_marker_dotplot.tsv")
stopifnot(all(dot_df$tnk_state %in% TNK_STATES))
dot_df$state <- factor(TNK_LABELS[match(dot_df$tnk_state, TNK_STATES)], levels = rev(TNK_LABELS))
dot_df$gene  <- factor(dot_df$gene, levels = unique(dot_df$gene))

p <- ggplot(dot_df, aes(gene, state, size = pct_exp, colour = avg_exp_scaled)) +
  geom_point() +
  scale_radius(range = c(0, 6), name = "Pct.exp") +
  scale_colour_gradient(low = "lightgrey", high = "blue", name = "Avg.exp") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(colour = "black"))

ggsave("marker_dotplot_tnk_states.pdf", p, width = 9.5, height = 3.5)
ggsave("marker_dotplot_tnk_states.png", p, width = 9.5, height = 3.5, dpi = 300, bg = "white")
