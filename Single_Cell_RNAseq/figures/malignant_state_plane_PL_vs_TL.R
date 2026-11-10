# =============================================================================
# Malignant cells of H3-DSG-PL and H3-DSG-TL in the plane of the cell states
#
# Input : ../data/cell_level/malignant_state_scores.tsv.gz   six module scores per cell
#         ../data/cell_level/cell_annotation.tsv.gz          group, sample, included_in_prop
# Output: malignant_state_plane_PL_vs_TL.pdf / .png in this directory
#
# Each point is a malignant cell. x = MES-like score - AC-like score; y = cycling score =
# larger of the G2M and S scores - larger of the MES-like and AC-like scores. Colour = the
# state with the highest score, with S and G2M merged into Cycling and OPC-like and OC-like
# into OPC/OC-like. The percentages are the share of each state among the cells drawn for
# the group.
#
# The cells are the 74,090 malignant cells flagged included_in_prop (all malignant cells,
# except that the five samples of the multi-site patient are down-sampled). The group with
# more cells (PL, 54,764 against 19,326) is then thinned to the size of the smaller one by
# a simple random draw, seed 123, so that both panels show 19,326 cells and their densities
# can be compared.
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

SEED <- 123
GROUP_LEVELS <- c("PL", "TL")
GROUP_LABELS <- c(PL = "H3-DSG-PL", TL = "H3-DSG-TL")
GROUP_COLORS <- c(PL = "#990000", TL = "#2b8cbe")
STATE_LEVELS <- c("Cycling", "OPC/OC-like", "AC-like", "MES-like")
STATE_FILL   <- c("MES-like" = "#c2a5cf", "AC-like" = "#edc9c566", "Cycling" = "#a94c3c66", "OPC/OC-like" = "#d7927066")
STATE_COLOUR <- c("MES-like" = "#c2a5cf", "AC-like" = "#edc9c5", "Cycling" = "#a94c3c", "OPC/OC-like" = "#d79270")

sc  <- read.delim(gzfile("../data/cell_level/malignant_state_scores.tsv.gz"), colClasses = c(cell_id = "character"))
ann <- read.delim(gzfile("../data/cell_level/cell_annotation.tsv.gz"), colClasses = c(cell_id = "character"))
j <- match(sc$cell_id, ann$cell_id)
stopifnot(!anyNA(j))
sc$group <- ann$group[j]
sc$renamed_id <- ann$renamed_id[j]
sc$included_in_prop <- ann$included_in_prop[j] %in% c(TRUE, "TRUE")

m <- sc[sc$included_in_prop, ]
stopifnot(nrow(m) == 74090)
m$state <- ifelse(m$state %in% c("S", "G2M"), "Cycling", ifelse(m$state %in% c("OC-like", "OPC-like"), "OPC/OC-like", m$state))
stopifnot(all(m$state %in% STATE_LEVELS))
m$x <- m$MES_like - m$AC_like
m$y <- pmax(m$G2M, m$S) - pmax(m$MES_like, m$AC_like)

# thin the larger group to the size of the smaller one (cells in the order of the table)
n_min <- min(table(m$group))
set.seed(SEED)
keep <- unlist(lapply(split(m$cell_id, m$group), function(cells) if (length(cells) > n_min) sample(cells, n_min) else cells),
               use.names = FALSE)
plot_df <- m[m$cell_id %in% keep, c("cell_id", "renamed_id", "group", "state", "x", "y")]
stopifnot(all(table(plot_df$group) == n_min))
plot_df$group <- factor(plot_df$group, levels = GROUP_LEVELS)
plot_df$state <- factor(plot_df$state, levels = STATE_LEVELS)

# share of each state among the cells drawn for the group, and where to print it
pct <- as.data.frame(prop.table(table(group = plot_df$group, state = plot_df$state), 1) * 100, responseName = "percent")
pct$label <- sprintf("%.1f%%", pct$percent)
pos <- data.frame(state = STATE_LEVELS, x = c(0, -0.15, -1.15, 0.85), y = c(1.85, 0.3, -0.8, -0.8))
pct <- merge(pct, pos, by = "state")

p <- ggplot(plot_df, aes(x, y)) +
  geom_point(aes(fill = state, colour = state), size = 0.1, shape = 21, alpha = 0.5) +
  geom_text(data = pct, aes(label = label), size = 3.6) +
  facet_grid(. ~ group, labeller = labeller(group = GROUP_LABELS)) +
  scale_fill_manual(values = STATE_FILL, breaks = STATE_LEVELS, name = NULL) +
  scale_colour_manual(values = STATE_COLOUR, breaks = STATE_LEVELS, name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1, fill = unname(STATE_COLOUR[STATE_LEVELS]))), fill = "none") +
  coord_cartesian(xlim = c(-2, 1.5)) +
  labs(x = "AC-like  <-->  MES-like", y = "Cycling score") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), legend.position = "bottom",
        axis.text = element_text(colour = "black"), panel.spacing = unit(0.6, "lines"),
        strip.text = element_text(colour = "black"))
# group names in the colour of the group
g <- ggplot_gtable(ggplot_build(p))
strip_idx <- grep("strip-t", g$layout$name)
for (i in seq_along(strip_idx)) {
  idx <- which(vapply(g$grobs[[strip_idx[i]]]$grobs[[1]]$children, function(x) inherits(x, "titleGrob"), TRUE))
  if (length(idx)) g$grobs[[strip_idx[i]]]$grobs[[1]]$children[[idx[1]]]$children[[1]]$gp$col <- GROUP_COLORS[[i]]
}

save_figure <- function(g, stem, width, height, dpi = 300) {
  pdf(paste0(stem, ".pdf"), width = width, height = height); grid::grid.draw(g); invisible(dev.off())
  png(paste0(stem, ".png"), width = width, height = height, units = "in", res = dpi, bg = "white")
  grid::grid.draw(g); invisible(dev.off())
}
save_figure(g, "malignant_state_plane_PL_vs_TL", width = 5.2, height = 3.6)
