# =============================================================================
# Proportion of cycling and MES-like cells among the malignant cells, PL vs TL
#
# Input : ../data/sample_level/malignant_state_composition.tsv
# Output: malignant_state_composition_PL_vs_TL.pdf / .png in this directory
#
# One point per patient. For each sample the percentage of its malignant cells in a state is
# taken over all malignant cells of the sample (S and G2M are added up as Cycling); the
# samples of one patient are then averaged. TL_6 .. TL_10 are five sites of one patient and
# give one point (15 patients: 8 PL, 7 TL). Two-sided Wilcoxon rank sum test between the
# patient values of PL and TL (exact: no ties); the P values are not adjusted.
#
# All malignant cells of a sample are used, not only the cells flagged included_in_prop:
# that flag thins the five samples of the multi-site patient for the plot of the cells
# (malignant_state_plane_PL_vs_TL.R) and is not used for a per-sample percentage.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

GROUP_LEVELS <- c("PL", "TL")
GROUP_COLORS <- c(PL = "#990000", TL = "#2b8cbe")
STATES_SHOWN <- c("Cycling", "MES-like")
ALL_STATES   <- c("Cycling", "OPC-like", "AC-like", "MES-like", "OC-like")

comp <- read.delim("../data/sample_level/malignant_state_composition.tsv", colClasses = c(patient_id = "character"))
comp$state <- ifelse(comp$category %in% c("S", "G2M"), "Cycling", comp$category)
stopifnot(all(comp$state %in% ALL_STATES))

# percentage of the sample's malignant cells in each state, 0 where the state is absent
samples <- unique(comp[, c("renamed_id", "group", "patient_id")])
sample_pct <- merge(samples, data.frame(state = ALL_STATES), by = NULL)
n_state <- tapply(comp$n, list(comp$renamed_id, comp$state), sum)
sample_pct$n <- mapply(function(r, s) { v <- n_state[r, s]; if (is.na(v)) 0 else v }, sample_pct$renamed_id, sample_pct$state)
sample_pct$percent <- 100 * sample_pct$n / ave(sample_pct$n, sample_pct$renamed_id, FUN = sum)

# per patient: mean of the percentages of the patient's samples
plot_df <- aggregate(percent ~ patient_id + group + state, data = sample_pct, FUN = mean)
plot_df <- merge(plot_df, aggregate(renamed_id ~ patient_id, data = samples, FUN = length), by = "patient_id")
names(plot_df)[names(plot_df) == "renamed_id"] <- "n_samples"
plot_df$state <- factor(plot_df$state, levels = ALL_STATES)
plot_df$group <- factor(plot_df$group, levels = GROUP_LEVELS)
stopifnot(all(table(plot_df$group[plot_df$state == "Cycling"]) == c(8, 7)))

# test per state
tests <- do.call(rbind, lapply(ALL_STATES, function(s) {
  d <- plot_df[plot_df$state == s, ]
  data.frame(state = s, n_PL = sum(d$group == "PL"), n_TL = sum(d$group == "TL"),
             p = wilcox.test(percent ~ group, data = d)$p.value)
}))

one_panel <- function(s, ylab = NULL, breaks = waiver()) {
  d <- plot_df[plot_df$state == s, ]
  p <- tests$p[tests$state == s]
  ggplot(d, aes(group, percent, fill = group)) +
    geom_boxplot(width = 0.7, linewidth = 0.8, colour = "black", outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.8, colour = "black") +
    scale_fill_manual(values = GROUP_COLORS) +
    scale_y_continuous(breaks = breaks, limits = c(0, NA)) +
    labs(x = NULL, y = ylab, title = bquote(atop(.(s), italic(P) ~ "=" ~ .(sprintf("%.1e", p))))) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 12),
          axis.text = element_text(colour = "black"))
}
p <- one_panel("Cycling", ylab = "Single cell proportion (%)") + one_panel("MES-like")

save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
}
save_figure(p, "malignant_state_composition_PL_vs_TL", width = 4.4, height = 3.6)
