# Composition of the T/NK states in PL and TL: one stacked bar per subtype.
#
# Reads  ../data/panel_data/tnk_cells.tsv.gz   (written by ../pipeline/99b_export_tnk_panels.R)
# Writes tnk_composition_PL_vs_TL.pdf and tnk_composition_PL_vs_TL.png next to this script.
#
# The cells of each subtype are pooled. The five samples of the patient with five sites (TL_6 to
# TL_10) contribute a random 20% of their cells each (flag in_pooled_bar, drawn with seed 1 in
# 99b_export_tnk_panels.R); every other sample contributes all of its cells.

suppressPackageStartupMessages(library(ggplot2))

# ---- T/NK states: object name, name shown in figures, colour (RColorBrewer "Paired") ----------
TNK_STATES  <- c("GZMK+ CD8 T", "CD8 Tex", "CD8 T_m", "CD4 T_m", "CD56-CD16+ NK", "Stress T",
                 "CD56+CD16- NK", "Prolif", "Treg", "Naive T", "IFN")
TNK_LABELS  <- c("GZMK+ CD8 T", "exhausted CD8 T", "CD8 memory T", "CD4 memory T", "CD56-CD16+ NK", "Stress T",
                 "CD56+CD16- NK", "proliferative T", "regulatory T", "Naive T", "Interferon T")
TNK_COLOURS <- setNames(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
                          "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"), TNK_LABELS)

cells <- read.delim(gzfile("../data/panel_data/tnk_cells.tsv.gz"), colClasses = c(cell_id = "character"))
stopifnot(all(cells$tnk_state %in% TNK_STATES))

# per sample and state: all cells, and the cells entering the bars (PL samples first, then TL)
smp <- unique(cells[, c("renamed_id", "group")])
smp <- smp[order(smp$group, as.integer(sub("^[A-Z]+_", "", smp$renamed_id))), ]
by_sample <- function(keep) table(factor(cells$renamed_id[keep], levels = smp$renamed_id),
                                  factor(cells$tnk_state[keep], levels = TNK_STATES))
per_sample <- data.frame(renamed_id = rep(smp$renamed_id, times = length(TNK_STATES)),
                         group = rep(smp$group, times = length(TNK_STATES)),
                         tnk_state = rep(TNK_STATES, each = nrow(smp)),
                         cells_total = as.vector(by_sample(rep(TRUE, nrow(cells)))),
                         cells_in_bar = as.vector(by_sample(cells$in_pooled_bar)))

# the plotted data: cells and percentage per subtype and state
bar_cells <- cells[cells$in_pooled_bar, ]
bar_df <- as.data.frame(table(group = factor(bar_cells$group, levels = c("PL", "TL")),
                              tnk_state = factor(bar_cells$tnk_state, levels = TNK_STATES)),
                        stringsAsFactors = FALSE)
names(bar_df)[3] <- "n"
bar_df$pct   <- 100 * bar_df$n / ave(bar_df$n, bar_df$group, FUN = sum)
bar_df$group <- factor(bar_df$group, levels = c("PL", "TL"))
bar_df$state <- factor(TNK_LABELS[match(bar_df$tnk_state, TNK_STATES)], levels = TNK_LABELS)

n_samples <- length(unique(cells$renamed_id))
title <- sprintf("T cells (N=%d samples, n=%d cells)", n_samples, nrow(cells))

p <- ggplot(bar_df, aes(group, pct, fill = state)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = TNK_COLOURS, name = NULL) +
  scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = title) +
  theme_classic() +
  theme(plot.title = element_text(size = 10, hjust = 0.5), axis.text = element_text(colour = "black"))

ggsave("tnk_composition_PL_vs_TL.pdf", p, width = 4.5, height = 5)
ggsave("tnk_composition_PL_vs_TL.png", p, width = 4.5, height = 5, dpi = 300, bg = "white")
