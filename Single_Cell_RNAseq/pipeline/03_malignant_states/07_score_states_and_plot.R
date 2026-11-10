# =============================================================================
# Malignant cell states: score the tumour cells on the six metaprograms, assign the state,
# plot, and tabulate the PL / TL composition
#
# Input : ../../objects/integrated/integrated_minf400_only_PL_6_4_cca_anno.rds   CCA-integrated
#             object of the 19 samples (132,610 cells), major cell types in integrated_snn_res.1_anno
#         ../../results/malignant_states/04_metaprograms/plateau_v1_all19/
#             cellular_states_metaprograms.txt, metaprogram_vs_reference_overlap.csv
#         ../../inputs/annotations/tumor_mp_state_names_plateau_v1_all19.csv   M1..M6 -> state name
#         ../../inputs/frozen/included_in_prop.csv   frozen TL_6 downsampling flag
# Output: ../../results/malignant_states/05_cellstates/
#             tumor_meta_<DATE_TAG>.csv            metadata of the tumour cells (81,304 rows)
#             tumor_states.csv                     scores and state of all cells (NA for non-tumour
#                                                  cells): the per-cell table the state columns of
#                                                  the final object are taken from
#             PL_TL_state_proportions_<DATE_TAG>.tsv   cells and % per state within PL / TL (long)
#         ../../results/malignant_states/04_metaprograms/plateau_v1_all19/
#             heatmap.jpeg, heatmap_6state_doheatmap.jpeg, triangle_plot_data{,_downsample}.csv
#             fig2g_state_triangle2{,_downsample}.{jpg,pdf} and fig2g_state_triangle2{,_downsample}_prop.csv
#                 (Group x state percentages of each of the two triangle plots)
#
# DATE_TAG is the date in the file names (environment variable DATE_TAG, default today,
# yyyymmdd). PLOT=0 in the environment scores and writes the tables but makes no plots.
#
# Outline
#   1. check that M1..M6 still mean what the annotation table says: the best K27scDSG match of
#      every metaprogram has to agree with the table, otherwise stop
#   2. drop cluster myeloid_mixed_T, scale the signature genes on all cells, subset the tumour
#      cells (integrated_snn_res.1_anno == "tumor")
#   3. AddModuleScore on the top 50 genes of each metaprogram (nbin = 24, ctrl = 100, seed = 1,
#      data layer; inputs in the order S / G2M / OPC-like / AC-like / MES-like / OC-like);
#      which.max = State2025, max = state.score
#   4. PL / TL state proportions
#   5-7. heatmap of the scaled signature genes by state; triangle plot (AC-like <-> MES-like
#      against a cycling score) by Group, also with PL and TL downsampled to the same number
#      of cells; Seurat DoHeatmap
#
# Runtime about 15 min (reading the integrated object dominates). Seurat 5.2.1, R 4.5.1.
# Paths are relative to this directory: run the script from here.
# =============================================================================

OBJ              <- "../../objects/integrated/integrated_minf400_only_PL_6_4_cca_anno.rds"
ANNOT_DIR        <- "../../inputs/annotations"
FROZEN_INCL_PROP <- "../../inputs/frozen/included_in_prop.csv"
RES_DIR    <- "../../results/malignant_states"
MP_DIR     <- file.path(RES_DIR, "04_metaprograms", "plateau_v1_all19")
OUT_DIR    <- file.path(RES_DIR, "05_cellstates")

DATE_TAG <- Sys.getenv("DATE_TAG", unset = format(Sys.Date(), "%Y%m%d"))
DO_PLOT  <- Sys.getenv("PLOT") != "0"

suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste(...))

P <- list(
  state_name_file  = "tumor_mp_state_names_plateau_v1_all19.csv",
  genes_per_module = 50,
  state_order      = c("S", "G2M", "OPC-like", "AC-like", "MES-like", "OC-like"),
  am_nbin = 24, am_ctrl = 100, am_seed = 1, am_slot = "data",
  tumor_anno_col = "integrated_snn_res.1_anno", tumor_label = "tumor",
  exclude_label  = "myeloid_mixed_T",
  incl_prop_csv  = FROZEN_INCL_PROP,          # frozen TL_6 downsampling flag
  heatmap_idents = "integrated_snn_res.2",   # subset(downsample = ) draws per ident; the active.ident of the original object was resolution 2
  heatmap_downsample = 5000, heatmap_shuffle_seed = 123,
  doheatmap_n = 5000, doheatmap_seed = 123,
  tri_downsample_seed = 123   # equal-size triangle plot: the larger group is drawn down to the size of the smaller
)

## ---- 1. metaprogram -> state name; check first that the numbering has not changed -------
state_map   <- read.csv(file.path(ANNOT_DIR, P$state_name_file), stringsAsFactors = FALSE,
                        check.names = FALSE, colClasses = "character")
state_names <- setNames(state_map$state, state_map$metaprogram)
ovl <- read.csv(file.path(MP_DIR, "metaprogram_vs_reference_overlap.csv"), row.names = 1, check.names = FALSE)
best <- setNames(colnames(ovl)[apply(ovl, 1, which.max)], rownames(ovl))
say("metaprogram M1..M6 best K27scDSG match:", paste(sprintf("%s=%s", names(best), best), collapse = " "))
if (!setequal(names(best), names(state_names)) || any(best[names(state_names)] != state_names))
  stop("The metaprogram numbering / meaning does not agree with ", P$state_name_file, ", rename them by hand first:\n",
       paste(sprintf("%s: annotation=%s reference_best=%s", names(state_names), state_names,
                     best[names(state_names)]), collapse = "\n"))

mods <- read.delim(file.path(MP_DIR, "cellular_states_metaprograms.txt"), check.names = FALSE, colClasses = "character")
colnames(mods) <- unname(state_names[colnames(mods)])
features <- lapply(mods, function(col) { v <- col[!is.na(col) & nzchar(trimws(col))]; unique(v)[seq_len(P$genes_per_module)] })

## ---- 2. object ---------------------------------------------------------------------------
say("loading annotated object")
obj <- readRDS(OBJ)
base_meta <- obj@meta.data                  # basic metadata of the 132,610 cells, used for the all-cell table in 3b
DefaultAssay(obj) <- "RNA"
obj <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"))   # only RNA is used; scale.data is computed below
obj <- subset(obj, cells = colnames(obj)[obj@meta.data[[P$tumor_anno_col]] != P$exclude_label])
Idents(obj) <- P$heatmap_idents
say("cells after dropping", P$exclude_label, ":", ncol(obj))

features_present <- lapply(features, function(g) intersect(g, rownames(obj)))[P$state_order]
say("genes present per state:", paste(sprintf("%s=%d", names(features_present), lengths(features_present)), collapse = " "))
obj <- ScaleData(obj, features = unique(unlist(features_present)), verbose = FALSE)   # scaled on all cells

tumor_cells <- colnames(obj)[obj@meta.data[[P$tumor_anno_col]] == P$tumor_label]
tumor_obj <- subset(obj, cells = tumor_cells)
rm(obj); invisible(gc())
say("tumor cells:", ncol(tumor_obj))

## ---- 3. scoring --------------------------------------------------------------------------
score_cols <- gsub("[^A-Za-z0-9]", "_", names(features_present))
tumor_obj <- AddModuleScore(tumor_obj, features = features_present, name = "MPSCORE",
                            nbin = P$am_nbin, ctrl = P$am_ctrl, seed = P$am_seed, slot = P$am_slot)
raw_cols <- paste0("MPSCORE", seq_along(features_present))
names(tumor_obj@meta.data)[match(raw_cols, names(tumor_obj@meta.data))] <- score_cols
md <- tumor_obj@meta.data
md$State2025   <- factor(names(features_present)[apply(md[, score_cols], 1, which.max)], levels = names(features_present))
md$state.score <- apply(md[, score_cols], 1, max)
md$cell.type.detailed <- as.character(md$State2025)
md$cell.type   <- ifelse(md$cell.type.detailed %in% c("S", "G2M"), "Cycling", md$cell.type.detailed)
md$State2025_tri <- ifelse(md$cell.type.detailed %in% c("S", "G2M"), "Cycling",
                    ifelse(md$cell.type.detailed %in% c("OC-like", "OPC-like"), "OC/OPC-like", md$cell.type.detailed))
md$patient_id  <- sub("^([^_]+_[^_]+)(?:_.*)?$", "\\1", md$orig.ident)

## TL_6 downsampling: TL_6_1..5 are five sites of one tumour. The draw that keeps 1,803 of their
## 9,017 tumour cells was made with sample() without a seed and cannot be redone, so the flag is
## read from the frozen included_in_prop.csv. It only affects the PL / TL proportion table and the
## triangle plots, not the scores and not the patient-level boxplots of
## 07b_state_composition_boxplot.R. The random sampling inside the plots uses fixed seeds.
ref <- read.csv(P$incl_prop_csv, row.names = "cell", check.names = FALSE)[, "included_in_prop", drop = FALSE]
md$included_in_prop <- as.logical(ref[rownames(md), "included_in_prop"])
stopifnot(!anyNA(md$included_in_prop))
tl6 <- md$orig.ident %in% paste0("TL_6_", 1:5)
say(sprintf("included_in_prop: TL_6_1..5 kept %d / %d; other samples all TRUE = %s",
            sum(md$included_in_prop[tl6]), sum(tl6), all(md$included_in_prop[!tl6])))
tumor_obj@meta.data <- md

meta_file <- file.path(OUT_DIR, sprintf("tumor_meta_%s.csv", DATE_TAG))
write.csv(md, meta_file)
say("wrote", meta_file, "(", nrow(md), "cells x", ncol(md), "cols )")

## ---- 3b. table of all cells ---------------------------------------------------------------
## All 132,610 cells with their basic metadata, plus the six scores / State2025 / state.score of
## the tumour cells (NA for the others). cell.type = State2025 of the tumour cells (S / G2M
## merged into Cycling), falling back to integrated_snn_res.1_anno for the other cells;
## cell.type.detailed is the same without merging Cycling; patient_id is defined for all cells.
SCORE_OUT_COLS <- c(score_cols, "State2025", "state.score")
full <- base_meta
full[, SCORE_OUT_COLS] <- NA
full[rownames(md), score_cols] <- md[, score_cols]
full[rownames(md), "state.score"] <- md$state.score
full$State2025 <- NA_character_
full[rownames(md), "State2025"] <- as.character(md$State2025)
full$cell.type <- as.character(full[[P$tumor_anno_col]])
full[rownames(md), "cell.type"] <- full[rownames(md), "State2025"]
full$cell.type.detailed <- full$cell.type
full$cell.type <- ifelse(full$cell.type %in% c("S", "G2M"), "Cycling", full$cell.type)
full$patient_id <- sub("^([^_]+_[^_]+)(?:_.*)?$", "\\1", full$orig.ident)
states_file <- file.path(OUT_DIR, "tumor_states.csv")
write.csv(full, states_file)
say("wrote", states_file)
rm(full, base_meta); invisible(gc())

## ---- 4. PL / TL state proportions --------------------------------------------------------
prop_long <- function(d, subset_name) {
  bind_rows(lapply(c("State2025", "cell.type", "State2025_tri"), function(col) {
    lv <- switch(col, State2025 = P$state_order,
                 cell.type = c("Cycling", "OPC-like", "AC-like", "MES-like", "OC-like"),
                 State2025_tri = c("Cycling", "OC/OPC-like", "AC-like", "MES-like"))
    tab <- table(Group = d$Group, state = factor(as.character(d[[col]]), levels = lv))
    data.frame(subset = subset_name, state_col = col, as.data.frame(tab, stringsAsFactors = FALSE)) %>%
      rename(n_cells = Freq) %>% group_by(Group) %>%
      mutate(group_total = sum(n_cells), pct_within_group = round(100 * n_cells / group_total, 2)) %>% ungroup()
  }))
}
props <- bind_rows(prop_long(md[md$included_in_prop, ], "tl6_subset_included_in_prop"),
                   prop_long(md, "all_tumor_cells"))
prop_file <- file.path(OUT_DIR, sprintf("PL_TL_state_proportions_%s.tsv", DATE_TAG))
write.table(props, prop_file, sep = "\t", quote = FALSE, row.names = FALSE)
say("wrote", prop_file)
for (col in c("State2025", "State2025_tri")) {
  cat("\n== % of tumor cells within Group,", col, ", TL_6 subset (included_in_prop) ==\n")
  w <- props %>% filter(subset == "tl6_subset_included_in_prop", state_col == col)
  print(round(xtabs(pct_within_group ~ Group + state, w), 2))
}

if (!DO_PLOT) { say("PLOT=0, done"); quit(save = "no") }
suppressPackageStartupMessages({ library(ComplexHeatmap); library(circlize); library(RColorBrewer) })
my_palette <- rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7",
                    "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"))

## ---- 5. heatmap.jpeg: ComplexHeatmap in blocks by State, colour scale of scale.data +-2 ----
tumor_10k <- subset(tumor_obj, downsample = P$heatmap_downsample)
scaled <- GetAssayData(tumor_10k, layer = "scale.data")
markers <- intersect(unique(unlist(features_present, use.names = FALSE)), rownames(scaled))
mat <- as.matrix(scaled[markers, ])
cell_states <- factor(tumor_10k$State2025, levels = P$state_order)
set.seed(P$heatmap_shuffle_seed)
cell_order <- unlist(lapply(levels(cell_states), function(st) { idx <- which(cell_states == st); idx[sample.int(length(idx))] }))
mat <- mat[, cell_order]; cell_states <- cell_states[cell_order]
state_pal <- setNames(c("#E31A1C", "#FF7F00", "#FDBF6F", "#B15928", "#6A3D9A", "#33A02C")[seq_along(levels(cell_states))],
                      levels(cell_states))
top_anno <- HeatmapAnnotation(State = cell_states, col = list(State = state_pal),
                              annotation_name_gp = gpar(fontsize = 7), show_legend = c(State = FALSE),
                              simple_anno_size = unit(3, "mm"))
ht <- Heatmap(mat, name = "scaled",
              col = colorRamp2(seq(-2, 2, length.out = length(my_palette)), my_palette),
              cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
              top_annotation = top_anno, column_split = cell_states, column_gap = unit(0.3, "mm"),
              column_title_gp = gpar(fontsize = 7), use_raster = TRUE, raster_quality = 3,
              heatmap_legend_param = list(direction = "horizontal"))
f_hm <- file.path(MP_DIR, "heatmap.jpeg")
jpeg(f_hm, width = 12, height = 8, units = "in", res = 300)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
dev.off()
say("saved", f_hm, "(", ncol(mat), "cells x", nrow(mat), "genes )")

## ---- 6. triangle plot (cells kept by the TL_6 downsampling) ------------------------------
m_tri <- md[md$included_in_prop, , drop = FALSE]
m_tri$State2 <- m_tri$State2025_tri
m_tri$AClikeScore  <- m_tri$AC_like
m_tri$MESlikeScore <- m_tri$MES_like
m_tri$cyclingScore <- pmax(m_tri$G2M, m_tri$S) - pmax(m_tri$MESlikeScore, m_tri$AClikeScore)
plot_tri <- function(d) ggplot(d, aes(x = MESlikeScore - AClikeScore, y = cyclingScore)) +
  geom_point(aes(fill = State2, color = State2), size = 0.1, pch = 21, alpha = 0.5) +
  facet_grid(~ Group) +
  scale_fill_manual(values  = c("MES-like" = "#c2a5cf", "AC-like" = "#edc9c566",
                                "Cycling" = "#a94c3c66", "OC/OPC-like" = "#d7927066")) +
  scale_color_manual(values = c("MES-like" = "#c2a5cf", "AC-like" = "#edc9c5",
                                "Cycling" = "#a94c3c", "OC/OPC-like" = "#d79270")) +
  coord_cartesian(xlim = c(-2, 1.5)) +
  theme_classic() + theme(strip.background = element_blank()) +
  labs(x = "AC-like  <-->  MES-like", y = "Cycling score")
## Companion table: cells per Group in the plot and the percentage of each State2 within the Group (wide)
TRI_STATES <- c("Cycling", "OC/OPC-like", "AC-like", "MES-like")
write_tri_prop <- function(d, stem) {
  stopifnot(all(d$State2 %in% TRI_STATES))
  tab <- table(d$Group, factor(d$State2, levels = TRI_STATES))
  pct <- round(100 * unclass(prop.table(tab, 1)), 2)
  out <- data.frame(Group = rownames(tab), n_cells = as.integer(rowSums(tab)),
                    matrix(pct, nrow = nrow(tab), dimnames = list(NULL, paste0(TRI_STATES, "_pct"))),
                    check.names = FALSE)
  f <- file.path(MP_DIR, paste0(stem, "_prop.csv"))
  write.csv(out, f, row.names = FALSE)
  say("wrote", f); print(out, row.names = FALSE)
}
p_tri <- plot_tri(m_tri)
write_tri_prop(m_tri, "fig2g_state_triangle2")
ggsave(file.path(MP_DIR, "fig2g_state_triangle2.jpg"), p_tri, width = 5, height = 2.3)
ggsave(file.path(MP_DIR, "fig2g_state_triangle2.pdf"), p_tri, width = 5, height = 2.3)
write.csv(m_tri, file.path(MP_DIR, "triangle_plot_data.csv"), row.names = TRUE)
say("saved fig2g_state_triangle2.{jpg,pdf} + triangle_plot_data.csv (", nrow(m_tri), "cells )")

## Equal-size version: of the same cells, the larger group is drawn down at random (simple random
## sample within the group, not by sample) to the size of the smaller one, so that the two facets
## have the same number of points and their densities can be compared directly. With the current
## data PL has 54,764 and TL 19,326 cells, so 19,326 PL cells are drawn.
n_min <- min(table(m_tri$Group))
set.seed(P$tri_downsample_seed)
keep_ds <- unlist(lapply(split(rownames(m_tri), m_tri$Group),
                         function(cells) if (length(cells) > n_min) sample(cells, n_min) else cells), use.names = FALSE)
m_tri_ds <- m_tri[rownames(m_tri) %in% keep_ds, , drop = FALSE]      # keeps the row order, so the stacking order of the points is the same as in the full version
stopifnot(all(table(m_tri_ds$Group) == n_min))
p_tri_ds <- plot_tri(m_tri_ds)
write_tri_prop(m_tri_ds, "fig2g_state_triangle2_downsample")
ggsave(file.path(MP_DIR, "fig2g_state_triangle2_downsample.jpg"), p_tri_ds, width = 5, height = 2.3)
ggsave(file.path(MP_DIR, "fig2g_state_triangle2_downsample.pdf"), p_tri_ds, width = 5, height = 2.3)
write.csv(m_tri_ds, file.path(MP_DIR, "triangle_plot_data_downsample.csv"), row.names = TRUE)
say("saved fig2g_state_triangle2_downsample.{jpg,pdf} + triangle_plot_data_downsample.csv (",
    paste(sprintf("%s=%d", names(table(m_tri_ds$Group)), table(m_tri_ds$Group)), collapse = " "), ")")

## ---- 7. heatmap_6state_doheatmap.jpeg ----------------------------------------------------
set.seed(P$doheatmap_seed)
heat_cells <- sample(colnames(tumor_obj), size = P$doheatmap_n, replace = FALSE)
gns <- unique(unlist(features_present))
state_colors <- setNames(c("#d53e4f", "#d53e4f", "#f46d43", "#fdae61", "#fcc5c0", "#fee08b"), P$state_order)
dhp <- DoHeatmap(tumor_obj, features = gns, group.by = "State2025", cells = heat_cells,
                 angle = 30, size = 3, group.colors = state_colors) +
  scale_fill_gradientn(colours = my_palette) +
  theme(axis.text.y = element_blank(), legend.position = "bottom") +
  guides(color = "none")
f_dhp <- file.path(MP_DIR, "heatmap_6state_doheatmap.jpeg")
ggsave(f_dhp, dhp, width = 12, height = 8)
say("saved", f_dhp)
say("done")
