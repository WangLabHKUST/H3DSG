# =============================================================================
# T/NK cell states in H3K27M diffuse spinal cord glioma (H3-DSG)
#
# Input : scRNA_final.rds, the released integrated object (132,343 cells, 19 samples).
#         T/NK cells are taken from the major cell type annotation column
#         integrated_snn_res.1_anno (4,454 cells).
# Output: scRNA_lymphocyte_final.rds  annotated T/NK object (4,201 cells, 11 states
#                                     in tnk_state). The other 253 T/NK cells are
#                                     dropped here; scRNA_final.rds keeps the state of
#                                     every cell in tnk_state and the reason a cell has
#                                     none in substate_excluded_reason.
#
# Outline
#   1. subset T/NK cells and drop samples with < 30 cells
#   2. per-sample CCA integration, 4000 variable features
#   3. res = 0.3 clustering, drop two clusters that mix glia / neuroblast with T/NK
#   4. res = 0.25 clustering, assign cell types from markers
#   5. split the interferon cluster (res = 0.05) and the large CD8 cluster (res = 0.15)
#   6. final manual states
#
# Cluster ids are hard coded in the maps below. They come from the run described in
# the paper; a re-run that produces a different number of clusters will stop at the
# stopifnot() guards rather than silently relabel cells.
#
# Runtime about 1 h on 40 cores. Seurat 5.2.1, R 4.5.1.
#
# Paths are relative to this directory: run the script from here. The input is the
# deposited object in ../../objects/released/; the output goes to ../../results/, so a
# re-run never overwrites the deposited scRNA_lymphocyte_final.rds it can be compared with.
# =============================================================================

library(Seurat)
library(ggplot2)

OBJ      <- '../../objects/released/scRNA_final.rds'
OUT_DIR  <- '../../results/lymphocyte_states'
OUT_RDS  <- file.path(OUT_DIR, 'scRNA_lymphocyte_final.rds')

ANNO_COL  <- 'integrated_snn_res.1_anno'   # major cell type annotation
MIN_CELLS <- 30                            # samples with fewer T/NK cells are dropped
NFEATURES <- 4000

# clusters at res = 0.3 that mix glia / neuroblast (GFAP, S100B) with T/NK markers
EXCLUDE_RES03 <- c('6', '10')

# res = 0.25 cluster -> cell type
CELLTYPE_MAP <- c('0' = 'CD8 T', '1' = 'memory T', '2' = 'CD56-CD16+ NK', '3' = 'Naive T',
                  '4' = 'IFN', '5' = 'CD56+CD16- NK', '6' = 'Treg', '7' = 'Prolif',
                  '8' = 'Stress T', '9' = 'CD8 T')

# sub-clusters of the interferon cluster -> cell type
IFN_MAP <- c('IFN_0' = 'CD8 T', 'IFN_1' = 'IFN', 'IFN_2' = 'IFN')

# final manual states; keys are res = 0.25 cluster ids, with sub-cluster suffixes
# added for the interferon cluster (4_x) and the large CD8 cluster (0_x)
MANUAL_MAP <- c('9' = 'CD8 Tex', '0_0' = 'CD8 T_m', '0_2' = 'CD8 T_m',
                '1' = 'CD4 T_m', '4_0' = 'CD4 T_m', '0_1' = 'GZMK+ CD8 T')
STATE_LEVELS <- c('GZMK+ CD8 T', 'CD8 Tex', 'CD8 T_m', 'CD4 T_m', 'CD56-CD16+ NK', 'Stress T',
                  'CD56+CD16- NK', 'Prolif', 'Treg', 'Naive T', 'IFN')

# 1. subset T/NK ---------------------------------------------------------------
obj <- readRDS(OBJ)
tnk_cells <- colnames(obj)[!is.na(obj@meta.data[[ANNO_COL]]) & obj@meta.data[[ANNO_COL]] == 'T/NK']
message('T/NK cells: ', length(tnk_cells))
stc <- subset(obj, cells = tnk_cells)
rm(obj); gc()

all_cells  <- colnames(stc)
orig_ident <- as.character(stc$orig.ident)

# only the raw counts are carried over; everything else is recomputed on the subset
DefaultAssay(stc) <- 'RNA'
for (a in setdiff(Assays(stc), 'RNA')) stc[[a]] <- NULL
stc@reductions <- list(); stc@graphs <- list(); stc@neighbors <- list()

n_per_sample <- table(orig_ident)
keep_samples <- names(n_per_sample)[n_per_sample >= MIN_CELLS]
dropped <- setdiff(names(n_per_sample), keep_samples)
if (length(dropped)) message('dropping samples with < ', MIN_CELLS, ' cells: ', paste(dropped, collapse = ', '))
stc <- subset(stc, orig.ident %in% keep_samples)
message('cells after sample filter: ', ncol(stc))

# 2. per-sample CCA integration ------------------------------------------------
stc.list <- SplitObject(stc, split.by = 'orig.ident')
rm(stc); gc()
stc.list <- lapply(stc.list, function(x) NormalizeData(x, verbose = FALSE))
stc.list <- lapply(stc.list, function(x)
  FindVariableFeatures(x, selection.method = 'vst', nfeatures = NFEATURES, verbose = FALSE))

features <- SelectIntegrationFeatures(object.list = stc.list, nfeatures = NFEATURES)
anchors  <- FindIntegrationAnchors(object.list = stc.list, anchor.features = features,
                                   dims = 1:30, k.filter = 10)
rm(stc.list); gc()
stc.combined <- IntegrateData(anchorset = anchors, dims = 1:30, k.weight = 10)
rm(anchors); gc()

DefaultAssay(stc.combined) <- 'integrated'
stc.combined <- ScaleData(stc.combined, verbose = FALSE)
stc.combined <- RunPCA(stc.combined, npcs = 30, seed.use = 42, verbose = FALSE)
stc.combined <- RunUMAP(stc.combined, reduction = 'pca', dims = 1:20, seed.use = 42, verbose = FALSE)
stc.combined <- FindNeighbors(stc.combined, reduction = 'pca', dims = 1:10, k.param = 20, verbose = FALSE)
stc.combined <- FindClusters(stc.combined, resolution = 0.4, random.seed = 0,
                             graph.name = 'integrated_snn', verbose = FALSE)

# 3. drop glia / neuroblast doublet clusters -----------------------------------
# reuses the SNN graph built above, no second FindNeighbors
stc.combined <- FindClusters(stc.combined, resolution = 0.3, random.seed = 0,
                             graph.name = 'integrated_snn', verbose = FALSE)
print(table(stc.combined$integrated_snn_res.0.3))
VlnPlot(stc.combined, features = c('CD3D', 'NKG7', 'GFAP', 'S100B'),
        group.by = 'integrated_snn_res.0.3', stack = TRUE, flip = TRUE)

stopifnot(all(EXCLUDE_RES03 %in% levels(stc.combined$integrated_snn_res.0.3)))
keep <- colnames(stc.combined)[!as.character(stc.combined$integrated_snn_res.0.3) %in% EXCLUDE_RES03]
stc.clean <- subset(stc.combined, cells = keep)
rm(stc.combined); gc()
message('cells after removing glia / neuroblast clusters: ', ncol(stc.clean))

# 4. cluster and annotate ------------------------------------------------------
stc.clean <- FindNeighbors(stc.clean, reduction = 'pca', dims = 1:20, k.param = 20, verbose = FALSE)
stc.clean <- RunUMAP(stc.clean, reduction = 'pca', dims = 1:10, seed.use = 42, verbose = FALSE)
stc.clean <- FindClusters(stc.clean, resolution = 0.25, random.seed = 0,
                          graph.name = 'integrated_snn', verbose = FALSE)
res_col <- 'integrated_snn_res.0.25'
print(table(stc.clean@meta.data[[res_col]]))

stopifnot(setequal(levels(stc.clean@meta.data[[res_col]]), names(CELLTYPE_MAP)))
stc.clean$celltype <- unname(CELLTYPE_MAP[as.character(stc.clean@meta.data[[res_col]])])

# 5. split the interferon cluster, then the large CD8 cluster ------------------
Idents(stc.clean) <- 'celltype'
stc.clean <- FindSubCluster(stc.clean, cluster = 'IFN', graph.name = 'integrated_snn',
                            subcluster.name = 'celltype_IFN_sub', resolution = 0.05)
stopifnot(all(grep('^IFN_', stc.clean$celltype_IFN_sub, value = TRUE) %in% names(IFN_MAP)))
ifn_sub <- as.character(stc.clean$celltype_IFN_sub)
stc.clean$celltype_IFN_sub_anno <- ifelse(ifn_sub %in% names(IFN_MAP), unname(IFN_MAP[ifn_sub]), ifn_sub)

# cluster ids carrying the interferon split, e.g. cluster 4 becomes 4_0 / 4_1 / 4_2
res_id  <- as.character(stc.clean@meta.data[[res_col]])
split_i <- startsWith(ifn_sub, 'IFN_')
res_sub <- res_id
res_sub[split_i] <- paste0(res_id[split_i], '_', sub('^.*_', '', ifn_sub[split_i]))
Idents(stc.clean) <- res_sub
stc.clean <- FindSubCluster(stc.clean, cluster = '0', graph.name = 'integrated_snn',
                            subcluster.name = 'cluster_final', resolution = 0.15, algorithm = 1)
grp <- as.character(stc.clean$cluster_final)
print(table(grp))

# 6. final states --------------------------------------------------------------
# clusters not listed in MANUAL_MAP keep the cell type of step 4/5
lab <- tapply(as.character(stc.clean$celltype_IFN_sub_anno), grp, function(x) names(which.max(table(x))))
stopifnot(all(names(MANUAL_MAP) %in% names(lab)))
final <- lab
final[names(MANUAL_MAP)] <- MANUAL_MAP
stopifnot(setequal(STATE_LEVELS, unique(final)))
stc.clean$tnk_state <- factor(unname(final[grp]), levels = STATE_LEVELS)
print(table(stc.clean$tnk_state))

DefaultAssay(stc.clean) <- 'RNA'
Idents(stc.clean) <- 'tnk_state'
saveRDS(stc.clean, OUT_RDS)

# why the other T/NK cells carry no state; the all-cell object records this in
# substate_excluded_reason
reason <- setNames(rep(NA_character_, length(all_cells)), all_cells)
reason[!(orig_ident %in% keep_samples)] <- 'sample_lt30cells'
reason[setdiff(all_cells[orig_ident %in% keep_samples], colnames(stc.clean))] <- 'res0.3_glia_neuroblast_mix'
print(table(reason, useNA = 'ifany'))
message('wrote ', OUT_RDS)
