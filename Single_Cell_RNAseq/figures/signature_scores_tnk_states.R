# Mean score of published T cell signatures (Zheng et al., Science 2021) in each T/NK state.
#
# Reads  ../data/panel_data/tnk_signature_scores.tsv   (written by ../pipeline/99b_export_tnk_panels.R)
# Writes signature_scores_tnk_states.pdf and signature_scores_tnk_states.png next to this script.
#
# tnk_signature_scores.tsv holds, for every state and every one of the 40 meta-cluster signatures,
# the mean AddModuleScore of the state's cells. Here each signature is z-scored across the 11
# states, and only the three signatures with the highest z in each state are drawn (30 of 40).
# z is truncated at +-3; rows and columns are ordered by hierarchical clustering of the drawn
# values (Euclidean distance, complete linkage, pheatmap defaults). States are named as in the object.

suppressPackageStartupMessages(library(pheatmap))

TNK_STATES <- c("GZMK+ CD8 T", "CD8 Tex", "CD8 T_m", "CD4 T_m", "CD56-CD16+ NK", "Stress T",
                "CD56+CD16- NK", "Prolif", "Treg", "Naive T", "IFN")
TOP_N <- 3

sc <- read.delim("../data/panel_data/tnk_signature_scores.tsv", check.names = FALSE)
stopifnot(all(sc$tnk_state %in% TNK_STATES))
signatures <- unique(sc$signature)

# state x signature matrix of mean module scores. Rows enter the clustering in alphabetical order of
# the state name, as they did in the analysis (rowsum() sorts them); the order of equal-distance
# leaves in the dendrogram depends on it.
states <- sort(TNK_STATES)
mean_mat <- matrix(sc$mean_score[match(paste(rep(states, times = length(signatures)),
                                              rep(signatures, each = length(states))),
                                        paste(sc$tnk_state, sc$signature))],
                   nrow = length(states), dimnames = list(states, signatures))
stopifnot(!anyNA(mean_mat))

z_all <- scale(mean_mat)                 # each signature across the states
z_all[is.na(z_all)] <- 0
top <- unique(unlist(apply(z_all, 1, function(r) names(sort(r, decreasing = TRUE))[seq_len(TOP_N)])))
z_shown <- z_all[, colnames(z_all) %in% top, drop = FALSE]      # columns keep their original order
z_plot  <- pmin(pmax(z_shown, -3), 3)

# page size grows with the labels, 15 pt cells
ov <- max(nchar(colnames(z_plot))) * 0.05
w  <- ncol(z_plot) * 15 / 72 + max(nchar(rownames(z_plot))) * 0.07 + 2 * ov + 2.6
h  <- nrow(z_plot) * 15 / 72 + ov + 1.3

pdf(NULL)                                # text is measured on a device; nothing is drawn here
ph <- pheatmap(z_plot, cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45,
               cellwidth = 15, cellheight = 15,
               color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
               breaks = seq(-3, 3, length.out = 101),
               main = "Mean signature score: T/NK states vs. Zheng et al. (Science, 2021)",
               silent = TRUE)
invisible(dev.off())

pdf("signature_scores_tnk_states.pdf", width = w, height = h)
grid::grid.draw(ph$gtable)
invisible(dev.off())
png("signature_scores_tnk_states.png", width = w, height = h, units = "in", res = 300, bg = "white")
grid::grid.newpage()
grid::grid.draw(ph$gtable)
invisible(dev.off())
