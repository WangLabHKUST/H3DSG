# =============================================================================
# Pseudotime trajectory across the five sampling sites of patient TL_6
#
# TL_6 was sampled at five sites along the cord (C1-C2 dorsal / ventral,
# C2-C3 dorsal / ventral and the C4 edge). This script fits a trajectory on the
# malignant cells of that patient and asks how far along it each site sits.
#
# Input : ../../objects/released/scRNA_final.rds
# Output, in ../../results/trajectory_TL6/:
#         umap_subset_seurat_clusters.jpg   malignant cells, re-embedded and clustered
#         knn10_prop_umap_subset.jpg        local site composition of each cell
#         meanCI_by_loc_pt_mean.jpg         mean pseudotime per site, 95% CI
#         tl6_malignant_cells.csv           per-cell table
#         pt_mean_by_location.csv, lineages.txt
#
# The trajectory is fitted in the CCA corrected PC space (30 PCs) with the res = 0.1
# clusters as slingshot cluster labels. The root cluster is cluster 1, the cluster
# with the lowest pseudotime and the highest share of OC-like cells. Set
# START_CLUS to fit from another root; output files then carry a _root<N> suffix.
#
# Runtime about 6 min. Seurat 5.2.1, slingshot, R 4.5.1.
#
# Paths are relative to this directory: run the script from here.
# =============================================================================

library(Seurat)
library(slingshot)
library(ggplot2)
library(patchwork)

OBJ      <- '../../objects/released/scRNA_final.rds'
OUT_DIR  <- '../../results/trajectory_TL6'
SAMPLES  <- read.delim('../../sample_metadata.tsv', colClasses = 'character')
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SIDS_TL6   <- paste0('TL_6_', 1:5)     # the five sites of patient SP73
N_PCS      <- 30                       # PCs used for the trajectory
UMAP_DIMS  <- 1:20; UMAP_NN <- 100; UMAP_MIN_DIST <- 0.5
CLUSTER_RES <- 0.1
KNN_K      <- 10
START_CLUS <- Sys.getenv('START_CLUS', unset = '1')
SUFFIX     <- if (identical(START_CLUS, '1')) '' else paste0('_root', START_CLUS)
sfx <- function(name) sub('(\\.[^.]+)$', paste0(SUFFIX, '\\1'), name)

# 1. malignant cells of TL_6 ---------------------------------------------------
obj <- readRDS(OBJ)
keep <- obj$in_cleaned_0806 & obj$orig.ident %in% SIDS_TL6 &
  obj$integrated_snn_res.1_anno == 'tumor'
obj <- subset(obj, cells = colnames(obj)[keep])
message('TL_6 malignant cells: ', ncol(obj))

# site labels: "C1-C2-Dorsal" -> section "C1-C2", side "Dorsal", short name "C1-C2-Dor"
sm <- SAMPLES[SAMPLES$sid %in% SIDS_TL6, ]
sm <- sm[order(as.integer(sm$location_sort)), ]          # cranial to caudal
tl <- unname(setNames(sm$tumor_location, sm$sid)[as.character(obj$orig.ident)])
parts <- strsplit(tl, '-', fixed = TRUE)
obj$sec_name  <- vapply(parts, function(p) paste(p[-length(p)], collapse = '-'), character(1))
obj$side_name <- vapply(parts, function(p) p[length(p)], character(1))
abbrev <- function(x) {
  x <- sub('-Dorsal$', '-Dor', x); x <- sub('-Ventral$', '-Ven', x)
  sub('^(C[0-9]+)-\\1-edge$', '\\1-edge', x, perl = TRUE)
}
obj$loc_name <- factor(abbrev(tl), levels = abbrev(sm$tumor_location))
print(table(obj$loc_name))

# 2. re-embed and cluster the malignant cells ----------------------------------
DefaultAssay(obj) <- 'integrated'
obj <- FindNeighbors(obj, reduction = 'pca', dims = 1:N_PCS, verbose = FALSE)
obj <- RunUMAP(obj, reduction = 'pca', dims = UMAP_DIMS, n.neighbors = UMAP_NN,
               min.dist = UMAP_MIN_DIST, reduction.name = 'umap_subset',
               reduction.key = 'UMAPsubset_', verbose = FALSE)
obj <- FindClusters(obj, resolution = CLUSTER_RES, verbose = FALSE)
print(table(Idents(obj)))

if (!nzchar(SUFFIX)) {
  p <- DimPlot(obj, reduction = 'umap_subset', label = TRUE, repel = TRUE)
  ggsave(file.path(OUT_DIR, 'umap_subset_seurat_clusters.jpg'), p, width = 6, height = 5)
}

# 3. site composition among the 10 nearest neighbours of each cell -------------
coords <- Embeddings(obj, 'umap_subset')
nn_idx <- RANN::nn2(coords, k = KNN_K + 1L)$nn.idx[, -1, drop = FALSE]   # column 1 is the cell itself
knn_cols <- character(0)
for (cn in c('side_name', 'sec_name')) {
  lab <- as.character(obj@meta.data[[cn]])
  for (v in sort(unique(lab))) {
    out <- sprintf('nn%d_prop_%s', KNN_K, gsub('[^A-Za-z0-9]+', '_', v))
    obj@meta.data[[out]] <- apply(nn_idx, 1, function(i) mean(lab[i] == v))
    knn_cols <- c(knn_cols, out)
  }
}
if (!nzchar(SUFFIX)) {
  coolwarm <- c('#3B4CC0', '#6788EE', '#9ABBFF', '#C9D7F0', '#EDD1C2', '#F7A889', '#E36745', '#B40426')
  p_knn <- lapply(knn_cols, function(cn)
    FeaturePlot(obj, reduction = 'umap_subset', features = cn, pt.size = 0.2) +
      scale_color_gradientn(colors = coolwarm, limits = c(0, 1), name = 'prop') + ggtitle(cn))
  ggsave(file.path(OUT_DIR, sprintf('knn%d_prop_umap_subset.jpg', KNN_K)), wrap_plots(p_knn, ncol = 3),
         width = 6 * min(3, length(p_knn)), height = 5 * ceiling(length(p_knn) / 3))
}

# 4. trajectory ----------------------------------------------------------------
sds <- slingshot(Embeddings(obj, 'pca')[, 1:N_PCS], clusterLabels = obj$seurat_clusters,
                 start.clus = START_CLUS)
lin <- slingLineages(sds)
print(lin)
writeLines(sprintf('%s: %s', names(lin), vapply(lin, paste, '', collapse = ' -> ')),
           file.path(OUT_DIR, sfx('lineages.txt')))
pt <- slingPseudotime(sds)
obj$pt_mean <- rowMeans(pt, na.rm = TRUE)   # cells belong to several lineages; average them

# 5. mean pseudotime per site --------------------------------------------------
d <- data.frame(loc_name = obj$loc_name, pt_mean = obj$pt_mean)
d <- d[!is.na(d$pt_mean), ]
summ <- do.call(rbind, lapply(split(d$pt_mean, droplevels(d$loc_name)), function(x)
  data.frame(n = length(x), mean = mean(x), se = sd(x) / sqrt(length(x)))))
summ$loc_name <- rownames(summ)
summ$lo <- summ$mean - 1.96 * summ$se
summ$hi <- summ$mean + 1.96 * summ$se
summ <- summ[order(summ$mean), c('loc_name', 'n', 'mean', 'se', 'lo', 'hi')]
write.csv(summ, file.path(OUT_DIR, sfx('pt_mean_by_location.csv')), row.names = FALSE)
print(summ, row.names = FALSE)

summ$loc_name <- factor(summ$loc_name, levels = summ$loc_name)
p_ci <- ggplot(summ, aes(x = mean, y = loc_name)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = 'y') +
  geom_point(aes(color = mean), size = 4) +
  geom_text(aes(label = sprintf('%.2f', mean)), hjust = 0.5, vjust = -0.9, size = 4) +
  scale_color_viridis_c() +
  labs(x = 'mean pseudotime (95% CI)', y = NULL, title = 'PT mean',
       subtitle = if (nzchar(SUFFIX)) paste0('start.clus = ', START_CLUS) else NULL) +
  theme_classic(base_size = 13) + theme(legend.position = 'none')
ggsave(file.path(OUT_DIR, sfx('meanCI_by_loc_pt_mean.jpg')), p_ci, width = 6, height = 4)

# 6. per-cell table ------------------------------------------------------------
cell_tab <- data.frame(cell_id = colnames(obj), orig.ident = obj$orig.ident,
                       loc_name = obj$loc_name, sec_name = obj$sec_name, side_name = obj$side_name,
                       cluster_res0.1 = obj$seurat_clusters, coords,
                       obj@meta.data[, knn_cols], pt, pt_mean = obj$pt_mean, check.names = FALSE)
write.csv(cell_tab, file.path(OUT_DIR, sfx('tl6_malignant_cells.csv')), row.names = FALSE)
message('done -> ', OUT_DIR)
