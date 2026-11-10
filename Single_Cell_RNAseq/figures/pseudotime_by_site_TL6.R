# =============================================================================
# Mean pseudotime of the malignant cells at each sampling site of the multi-site patient
#
# Input : ../data/panel_data/pseudotime_TL6_malignant_cells.tsv
# Output: pseudotime_by_site_TL6.pdf / .png in this directory
#
# One point per sampling site: the mean slingshot pseudotime of the malignant cells of the
# site (pt_mean, the average over the lineages a cell belongs to), with the 95% confidence
# interval mean +/- 1.96 x s.e.m. (cells are the unit, n = cells of the site). Sites are
# ordered by mean pseudotime. The bracket between two neighbouring sites gives the P value
# of a two-sided Welch t-test between the cells of the two sites ("ns" when P >= 0.05).
# The sites are five samples of one tumour: P values test cells, not patients, and are not
# adjusted for the four comparisons.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

pt <- read.delim("../data/panel_data/pseudotime_TL6_malignant_cells.tsv", colClasses = c(patient_id = "character"))
stopifnot(length(unique(pt$patient_id)) == 1, !anyNA(pt$pt_mean))

# site label of the figure: patient and site number (1 = most cranial), e.g. "SP73-1"
smp <- read.delim("../data/panel_data/sample_resection_order.tsv", colClasses = c(patient_id = "character"))
pt$site_label <- paste(pt$patient_id, smp$site_number[match(pt$renamed_id, smp$renamed_id)], sep = "-")

# mean, s.e.m. and 95% CI per site; sites ordered by mean (lowest first)
summ <- do.call(rbind, lapply(split(pt, pt$site_label), function(d) data.frame(
  site = d$site_label[1], renamed_id = d$renamed_id[1], tumor_location = d$tumor_location[1],
  n = nrow(d), mean = mean(d$pt_mean), se = sd(d$pt_mean) / sqrt(nrow(d)))))
summ$lo <- summ$mean - 1.96 * summ$se
summ$hi <- summ$mean + 1.96 * summ$se
summ <- summ[order(summ$mean), ]
rownames(summ) <- NULL
summ$site <- factor(summ$site, levels = summ$site)

# Welch t-test between each site and the next one up
tests <- do.call(rbind, lapply(2:nrow(summ), function(i) {
  a <- pt$pt_mean[pt$site_label == summ$site[i]]; b <- pt$pt_mean[pt$site_label == summ$site[i - 1]]
  tt <- t.test(a, b)
  data.frame(upper = as.character(summ$site[i]), lower = as.character(summ$site[i - 1]),
             y_upper = i, y_lower = i - 1, difference = mean(a) - mean(b),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value)
}))
tests$label <- ifelse(tests$p < 0.05, sprintf("%.1e", tests$p), "ns")

x_br <- max(summ$hi) + 0.35        # x of the brackets
p <- ggplot(summ, aes(x = mean, y = site)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y", linewidth = 0.4) +
  geom_point(aes(colour = mean), size = 3.2) +
  geom_text(aes(label = sprintf("%.2f", mean)), vjust = -1, size = 3.5) +
  geom_segment(data = tests, aes(x = x_br, xend = x_br, y = y_lower + 0.1, yend = y_upper - 0.1), inherit.aes = FALSE) +
  geom_segment(data = tests, aes(x = x_br - 0.08, xend = x_br, y = y_lower + 0.1, yend = y_lower + 0.1), inherit.aes = FALSE) +
  geom_segment(data = tests, aes(x = x_br - 0.08, xend = x_br, y = y_upper - 0.1, yend = y_upper - 0.1), inherit.aes = FALSE) +
  geom_text(data = tests, aes(x = x_br + 0.12, y = (y_lower + y_upper) / 2, label = label),
            inherit.aes = FALSE, hjust = 0, size = 3.5) +
  scale_colour_viridis_c() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.02))) +
  coord_cartesian(xlim = c(min(summ$lo) - 0.2, x_br + 1.1), clip = "off") +
  labs(x = "Mean pseudotime (95% CI)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", axis.text = element_text(colour = "black"),
        plot.margin = margin(6, 6, 6, 6))

save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
}
save_figure(p, "pseudotime_by_site_TL6", width = 4.2, height = 3)
