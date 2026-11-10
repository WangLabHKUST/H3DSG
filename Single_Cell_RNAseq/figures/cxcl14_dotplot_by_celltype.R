# =============================================================================
# CXCL14 expression by cell type and subtype (PL / TL)
#
# Input : ../data/panel_data/cxcl14_dotplot_data.csv
# Output: cxcl14_dotplot_by_celltype.pdf / .png in this directory
#
# One dot per cell type and group. Colour = avg.exp, the mean of the de-logged CXCL14
# expression; size = pct.exp, the percentage of cells with CXCL14 > 0 (radius proportional
# to pct.exp, capped at 30%). OPC-like and OC-like are merged into OPC/OC-like, AC-like and
# MES-like into AC/MES-like, G2M and S into Cycling.
#
# The table holds two versions of the malignant states. VERSION = "cell.type" is the state
# assignment of the released objects; "cell.type_0811" is an earlier one that the released
# objects do not carry.
# =============================================================================

source("panel_helpers.R")

VERSION <- "cell.type"
PCT_MAX <- 30

dots <- read.csv("../data/panel_data/cxcl14_dotplot_data.csv", stringsAsFactors = FALSE)
stopifnot(VERSION %in% dots$version)

# raw state -> label on the panel, top to bottom
STATE_LABELS <- c("OPC/OC" = "OPC/OC-like", "AC/MES" = "AC/MES-like", "Cycling" = "Cycling",
                  "vascular" = "Vascular", "T/NK" = "T", "B" = "B", "myeloid" = "Myeloid",
                  "neutrophil" = "Neutrophil", "oligo" = "Oligodendrocyte")

plot_df <- dots[dots$version == VERSION, ]
stopifnot(setequal(plot_df$state, names(STATE_LABELS)), !anyDuplicated(plot_df[c("state", "Group")]),
          all(plot_df$pct.exp <= PCT_MAX))
plot_df$state_label <- factor(STATE_LABELS[plot_df$state], levels = rev(STATE_LABELS))
plot_df$group_label <- factor(GROUP_LABELS[plot_df$Group], levels = GROUP_LABELS)

p <- ggplot(plot_df) +
  geom_point(aes(x = group_label, y = state_label, colour = avg.exp, size = pct.exp)) +
  scale_colour_gradient2(low = "white", high = "blue", name = "Avg.exp") +
  scale_radius(range = c(0, 6), limits = c(0, PCT_MAX), breaks = c(0, PCT_MAX / 2, PCT_MAX), name = "Pct.exp") +
  guides(size = guide_legend(order = 1),
         colour = guide_colourbar(order = 2, theme = theme(legend.key.height = unit(0.8, "in"),
                                                            legend.key.width = unit(0.15, "in")))) +
  theme_classic() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1, colour = unname(GROUP_COLORS[names(GROUP_LABELS)])),
        strip.background = element_blank()) +
  labs(x = "", y = "")

save_figure(p, "cxcl14_dotplot_by_celltype", width = 3.4, height = 3.2)
