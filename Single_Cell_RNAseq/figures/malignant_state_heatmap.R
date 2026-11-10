# =============================================================================
# Heat map of the malignant cell states and their signature genes
#
# Input : ../data/panel_data/malignant_state_heatmap_5000cells.tsv.gz
# Output: malignant_state_heatmap.pdf / .png in this directory
#
# Columns are 5,000 malignant cells drawn at random (seed 123) from the 81,304 malignant
# cells; rows are the 297 genes of the six metaprogram signatures (top 50 genes of each,
# genes shared by two signatures once). Cells are grouped by state (S, G2M, OPC-like,
# AC-like, MES-like, OC-like) and the genes by the signature they come from, in the same
# order. Values are z-scores of the log-normalised expression, scaled over all 132,343
# cells of the cohort and clipped at -2.5 and 2.5 (as Seurat's DoHeatmap draws them). The
# table holds the column position of each cell (draw_position, the gaps between the state
# blocks are empty columns); a cell is named by its sample and 16-nt barcode. Only a few genes per state are named.
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

STATE_ORDER  <- c("S", "G2M", "OPC-like", "AC-like", "MES-like", "OC-like")
STATE_COLORS <- c("S" = "#d53e4f", "G2M" = "#d53e4f", "OPC-like" = "#f46d43",
                  "AC-like" = "#fdae61", "MES-like" = "#fcc5c0", "OC-like" = "#fee08b")
HEAT_COLORS  <- rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7",
                      "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"))
LABEL_GENES  <- c("PCNA", "MCM7", "BRCA1", "TOP2A", "MKI67", "CDK1", "CSPG5", "PDGFRA", "OLIG2",
                  "APOE", "GFAP", "AQP4", "CD44", "CAV1", "CD63", "MOG", "PLP1", "MBP")

heat <- read.delim(gzfile("../data/panel_data/malignant_state_heatmap_5000cells.tsv.gz"),
                   check.names = FALSE, colClasses = c(cell = "character", sid = "character"))
meta_cols <- c("cell", "draw_position", "sid", "renamed_id", "state")
genes <- setdiff(colnames(heat), meta_cols)                      # top to bottom
heat$state <- factor(heat$state, levels = STATE_ORDER)
stopifnot(nrow(heat) == 5000, length(genes) == 297, all(LABEL_GENES %in% genes), !anyNA(heat$state))

long <- data.frame(x = rep(heat$draw_position, times = length(genes)),
                   y = rep(seq(length(genes), 1), each = nrow(heat)),          # first gene on top
                   value = as.vector(as.matrix(heat[, genes])))
n_col <- max(heat$draw_position)

# extent and centre of each state block
blocks <- do.call(rbind, lapply(STATE_ORDER, function(s) {
  x <- heat$draw_position[heat$state == s]
  data.frame(state = s, xmin = min(x) - 0.5, xmax = max(x) + 0.5)
}))
blocks$mid <- (blocks$xmin + blocks$xmax) / 2
# S and G2M are both "Cycling" in the labels
labels <- data.frame(label = c("Cycling", STATE_ORDER[-(1:2)]),
                     x = c((blocks$xmin[1] + blocks$xmax[2]) / 2, blocks$mid[-(1:2)]))

# rows of each signature (a gene in two signatures belongs to the first); the three genes named
# for a state are spread evenly over the rows of its block
sig <- read.delim("../data/metaprograms/metaprogram_signature_genes.tsv")
row_state <- vapply(genes, function(g) STATE_ORDER[min(match(sig$state[sig$gene == g], STATE_ORDER))], "")
stopifnot(!anyNA(row_state), !is.unsorted(match(row_state, STATE_ORDER)))
gene_labels <- do.call(rbind, lapply(STATE_ORDER, function(s) {
  named <- LABEL_GENES[row_state[LABEL_GENES] == s]
  rows <- which(row_state == s)                                        # positions from the top
  centre <- rows[1] - 1 + (seq_along(named) - 0.5) * length(rows) / length(named)
  data.frame(gene = named, y = length(genes) - centre + 1)
}))

p_bar <- ggplot() +
  geom_rect(data = blocks, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1, fill = state)) +
  geom_text(data = labels, aes(x = x, y = 1.4, label = label), angle = 30, hjust = 0, size = 3.2) +
  scale_fill_manual(values = STATE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(0.5, n_col + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3), expand = c(0, 0)) +
  coord_cartesian(clip = "off") + theme_void() + theme(plot.margin = margin(0, 28, 0, 0))

p_heat <- ggplot(long, aes(x, y, fill = value)) +
  geom_raster() +
  geom_vline(xintercept = blocks$xmax[-nrow(blocks)] + 0.5, colour = "grey40", linewidth = 0.3) +
  geom_text(data = gene_labels, aes(x = n_col + 20, y = y, label = gene), inherit.aes = FALSE, hjust = 0, size = 3) +
  scale_fill_gradientn(colours = HEAT_COLORS, limits = c(-2.5, 2.5), na.value = "white", name = "Exp.") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, n_col + 0.5), ylim = c(0.5, length(genes) + 0.5), clip = "off") +
  theme_void() + theme(legend.position = "right", legend.box.spacing = unit(1.5, "cm"), plot.margin = margin(0, 28, 0, 0))

p <- p_bar / p_heat + plot_layout(heights = c(0.12, 1))

save_figure <- function(p, stem, width, height, dpi = 300) {
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
}
save_figure(p, "malignant_state_heatmap", width = 6.4, height = 5)
