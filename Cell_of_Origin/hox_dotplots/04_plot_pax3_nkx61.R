# =============================================================================
# hox_dotplots/04_plot_pax3_nkx61.R -- PAX3 and NKX6-1: public brain DMG vs H3-DSG-PL / H3-DSG-TL
#
# Input : ../data/hox_dotplots/pax3_nkx61_dotplot.tsv
# Output: pax3_nkx61_celltype_dotplot.pdf     4 x 3 in; the four groups x immune / stromal / malignant cells
#         pax3_nkx61_malignant_dotplot.pdf    3 x 3 in; malignant cells only, the four groups x two genes
#
# Groups: pontine and thalamic = the public pontine and thalamic DMG cells; PL and TL = the
# H3-DSG-PL and H3-DSG-TL samples of this cohort. Colour = average log-normalised expression,
# dot area = percentage of cells in which the gene is detected.
# Malignant-cell figure: percentages above 30 are drawn as 30; average expression above
# 0.35 takes the top colour.
#
# This script reads the table only and opens no object.
# =============================================================================
source("00_config.R")
source("dotplot_df.R")
suppressPackageStartupMessages(library(ggplot2))

pax <- read_table("pax3_nkx61_dotplot.tsv")
panel_levels <- c("pontine", "thalamic", "PL", "TL")       # first group at the bottom of the plot
stopifnot(setequal(pax$group, panel_levels), setequal(pax$celltype, CELLTYPES), setequal(pax$gene, c("PAX3", "NKX6-1")))
pax_df <- data.frame(features.plot = pax$gene,
                     id = factor(pax$celltype, levels = CELLTYPES),
                     avg.exp = pax$avg_exp, pct.exp = pax$pct_exp,
                     orig.ident = factor(pax$group, levels = panel_levels),
                     stringsAsFactors = FALSE)

# ---- by cell type group ------------------------------------------------------------------
p_celltype <- DotPlotDF(
  pax_df,
  group.by = "id", facet.by = "orig.ident", facet.nrow = 1, flip = TRUE, sort.levels = TRUE,
  cols = c("white", "blue"), scale.by = "size", dot.scale = 6, x.text.angle = 0,
  scale = FALSE, log.transform = FALSE,
  color.name = "Avg. Exp", size.name = "Pct. Exp"
) + theme(axis.text.y = element_text(hjust = 1), axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
ggsave("pax3_nkx61_celltype_dotplot.pdf", p_celltype, width = 4, height = 3, dpi = 300)
message("saved pax3_nkx61_celltype_dotplot.pdf")

# ---- malignant cells only (Seurat layout: x = gene, y = group) -------------------------------
p_malignant <- DotPlotDF(
  pax_df[pax_df$id == "malignant", ],
  group.by = "orig.ident", sort.levels = TRUE,
  features = c("PAX3", "NKX6-1"),     # PAX3 on the left, NKX6-1 on the right
  scale = FALSE, log.transform = FALSE,
  scale.by = "size", dot.scale = 6, scale.max = 30, size.limits = c(0, 30),   # pct > 30 is drawn as 30
  x.text.angle = 0,
  cols = RColorBrewer::brewer.pal(5, "Reds"),
  color.limits = c(0, 0.35),          # larger values take the top colour (squish)
  color.name = "Avg. Exp", size.name = "Pct. Exp"
) + theme(axis.text.y = element_text(hjust = 1))
ggsave("pax3_nkx61_malignant_dotplot.pdf", p_malignant, width = 3, height = 3, dpi = 300)
message("saved pax3_nkx61_malignant_dotplot.pdf")
