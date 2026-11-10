# =============================================================================
# PAX3, NKX6-1 and HOX gene expression by cell type at the five sampling sites of the
# multi-site patient
#
# Input : ../data/panel_data/pax3_nkx61_dotplot_stats_by_sample_celltype.tsv
#         ../data/panel_data/hox_dotplot_stats_by_sample_celltype.tsv
#         ../data/panel_data/hox_gene_annotation.tsv
#         ../data/panel_data/sample_resection_order.tsv
#         ../data/panel_data/pseudotime_TL6_malignant_cells.tsv   (order of the sites)
# Output: pax3_nkx61_hox_dotplot_by_celltype_site.pdf / .png in this directory
#
# One panel per sampling site, in the order of increasing mean pseudotime of its malignant
# cells; columns are immune, stromal and malignant cells. Upper block: PAX3 and NKX6-1,
# colour = mean log-normalised expression over all cells of the group (red). Lower block:
# HOX genes, same measure (blue); the genes are those of the HOX display order that are detected (> 0 % of cells) in at
# least one group of these five samples. Dot area = percentage of cells expressing the gene.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

CELLTYPE_LEVELS <- c("immune", "stromal", "malignant")
CELLTYPE_LABELS <- c(immune = "Immune", stromal = "Stromal", malignant = "Malignant")
AVG_MAX_PAX <- 0.24     # upper end of the PAX3 / NKX6-1 colour scale (values above are drawn at the end colour)
AVG_MAX_HOX <- 0.12     # same for the HOX genes
PAX_COLORS  <- c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d")   # RColorBrewer "Reds"

pax <- read.delim("../data/panel_data/pax3_nkx61_dotplot_stats_by_sample_celltype.tsv")
hox <- read.delim("../data/panel_data/hox_dotplot_stats_by_sample_celltype.tsv")
ann <- read.delim("../data/panel_data/hox_gene_annotation.tsv")
smp <- read.delim("../data/panel_data/sample_resection_order.tsv", colClasses = c(patient_id = "character"))
pt  <- read.delim("../data/panel_data/pseudotime_TL6_malignant_cells.tsv", colClasses = c(patient_id = "character"))

# the samples of the multi-site patient, labelled "<patient>-<site number>", ordered by pseudotime
sites <- smp[smp$patient_id %in% pt$patient_id, ]
sites$site <- paste(sites$patient_id, sites$site_number, sep = "-")
sites$mean_pseudotime <- tapply(pt$pt_mean, pt$renamed_id, mean)[sites$renamed_id]
sites <- sites[order(sites$mean_pseudotime), ]
site_levels <- sites$site
stopifnot(nrow(sites) == 5)

prep <- function(d) {
  d <- d[d$renamed_id %in% sites$renamed_id, ]
  d$site <- factor(sites$site[match(d$renamed_id, sites$renamed_id)], levels = site_levels)
  d$celltype <- factor(d$celltype, levels = CELLTYPE_LEVELS, labels = CELLTYPE_LABELS[CELLTYPE_LEVELS])
  d
}
pax_df <- prep(pax)
pax_df$gene <- factor(pax_df$gene, levels = c("PAX3", "NKX6-1"))        # NKX6-1 on the upper row

hox_df <- prep(hox)
hox_df <- hox_df[hox_df$gene %in% ann$gene[!is.na(ann$order_multisite)], ]
detected <- tapply(hox_df$pct_exp, hox_df$gene, max) > 0
hox_genes <- ann$gene[!is.na(ann$order_multisite)]
hox_genes <- hox_genes[order(ann$order_multisite[match(hox_genes, ann$gene)])]
hox_genes <- hox_genes[detected[hox_genes]]
hox_df <- hox_df[hox_df$gene %in% hox_genes, ]
hox_df$gene <- factor(hox_df$gene, levels = rev(hox_genes))             # first gene on top

size_scale <- scale_size(range = c(0, 7), limits = c(0, 32), breaks = c(0, 5, 10, 15, 20), name = "% expr.")
base <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(colour = "black"), strip.background = element_blank(),
        panel.spacing = unit(0.6, "lines"))

p_pax <- ggplot(pax_df, aes(celltype, gene)) +
  geom_point(aes(size = pct_exp, colour = mean_log_expr)) +
  scale_colour_gradientn(colours = PAX_COLORS, limits = c(0, AVG_MAX_PAX), oob = scales::squish,
                         name = "NKX6-1/PAX3\navg. expr.") +
  size_scale + facet_grid(. ~ site) + labs(x = NULL, y = NULL) + base +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

p_hox <- ggplot(hox_df, aes(celltype, gene)) +
  geom_point(aes(size = pct_exp, colour = mean_log_expr)) +
  scale_colour_gradient(low = "lightgray", high = "#1a08fe", limits = c(0, AVG_MAX_HOX), oob = scales::squish,
                        name = "HOX genes\navg. expr.") +
  size_scale + facet_grid(. ~ site) + labs(x = NULL, y = NULL) + base +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_blank())

p <- (p_pax / p_hox) + plot_layout(heights = c(2, 23), guides = "collect")

save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
}
save_figure(p, "pax3_nkx61_hox_dotplot_by_celltype_site", width = 7.5, height = 8)
