# =============================================================================
# Cell type composition of every sample, with the group average
#
# Input : ../data/sample_level/cell_type_composition.tsv   n cells per sample and cell type
# Output: celltype_composition_by_sample.pdf / .png in this directory
#
# One stacked bar per sample, the seven major cell types as a percentage of all cells of the
# sample (the 132,343 cells of the released object; every cell counts, malignant or not).
# "Average" is the unweighted mean of the sample percentages of the group: 8 PL samples and
# 11 TL samples, where TL_6 .. TL_10 are five sites of one tumour and each counts as a sample.
# =============================================================================

source("panel_helpers.R")
suppressPackageStartupMessages(library(patchwork))

comp <- read.delim("../data/sample_level/cell_type_composition.tsv", colClasses = c(patient_id = "character"))

# every sample x cell type; a type that is absent from a sample has 0 cells
samples <- unique(comp[, c("renamed_id", "group", "patient_id")])
stopifnot(!anyDuplicated(samples$renamed_id))
grid <- merge(samples, data.frame(category = names(CELLTYPE_LABELS)), by = NULL)
grid$n <- comp$n[match(paste(grid$renamed_id, grid$category), paste(comp$renamed_id, comp$category))]
grid$n[is.na(grid$n)] <- 0L
stopifnot(all(comp$category %in% names(CELLTYPE_LABELS)))
grid$percent <- 100 * grid$n / ave(grid$n, grid$renamed_id, FUN = sum)
grid$celltype <- celltype_label(grid$category)

# group averages
avg <- aggregate(percent ~ group + celltype, data = grid, FUN = mean)
avg$renamed_id <- "Average"
bars <- rbind(grid[, c("renamed_id", "group", "celltype", "percent")], avg[, c("renamed_id", "group", "celltype", "percent")])

plot_df <- bars
stopifnot(all(abs(tapply(plot_df$percent, plot_df[c("renamed_id", "group")], sum) - 100) < 1e-8, na.rm = TRUE))

draw_group <- function(g, show_y) {
  d <- plot_df[plot_df$group == g, ]
  d$renamed_id <- factor(d$renamed_id, levels = c(order_samples(d$renamed_id[d$renamed_id != "Average"]), "Average"))
  ggplot(d, aes(renamed_id, percent, fill = celltype)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.92, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = CELLTYPE_COLORS, name = NULL) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = if (show_y) "Single cell proportion (%)" else NULL, title = GROUP_LABELS[[g]]) +
    theme_panel() +
    theme(plot.title = element_text(colour = GROUP_COLORS[[g]], size = 13),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
          axis.text.y = if (show_y) element_text(size = 11) else element_blank(),
          axis.ticks.y = if (show_y) element_line() else element_blank(),
          axis.line.y = if (show_y) element_line() else element_blank(),
          axis.title.y = element_text(size = 13), plot.margin = margin(10, 4, 10, 4))
}

n_pl <- length(unique(plot_df$renamed_id[plot_df$group == "PL"]))
n_tl <- length(unique(plot_df$renamed_id[plot_df$group == "TL"]))
p <- (draw_group("PL", TRUE) + draw_group("TL", FALSE)) +
  plot_layout(widths = c(n_pl, n_tl), guides = "collect") & theme(legend.position = "right")

save_figure(p, "celltype_composition_by_sample", width = (n_pl + n_tl) * 0.45 + 1.6, height = 6)
