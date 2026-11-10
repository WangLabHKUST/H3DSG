# =============================================================================
# Myeloid cell states in H3K27M diffuse spinal cord glioma (H3-DSG)
#
# Input : scRNA_final.rds, the released integrated object (132,343 cells, 19 samples).
#         Myeloid cells are the cells labelled myeloid or neutrophil in
#         integrated_snn_res.1_anno that also passed the global doublet / low quality
#         filter (in_cleaned_0806), 42,232 cells.
# Output: scRNA_myeloid_final.rds  annotated myeloid object (39,573 cells;
#                                  myeloid_state_coarse, 8 states, and
#                                  myeloid_state_fine, 12 states). The other 2,659
#                                  myeloid cells are dropped here; scRNA_final.rds keeps
#                                  the state of every cell in the same two columns and
#                                  the reason a cell has none in substate_excluded_reason.
#
# Outline
#   1. subset myeloid / neutrophil cells and drop samples with < 30 cells
#   2. per-sample CCA integration, 2000 variable features
#   3. res = 0.3 clustering, drop one cluster that mixes tumour cells with myeloid
#   4. re-cluster at res = 0.3, split the monocyte / dendritic cluster (res = 0.1)
#   5. fine and coarse state labels; two fine clusters are technical artefacts and
#      are removed
#
# Cluster ids are hard coded in the maps below. They come from the run described in
# the paper; a re-run that produces a different number of clusters will stop at the
# stopifnot() guards rather than silently relabel cells.
#
# Runtime about 6 h on 40 cores, dominated by the integration step.
# Seurat 5.2.1, R 4.5.1.
#
# Paths are relative to this directory: run the script from here. The input is the
# deposited object in ../../objects/released/; the output goes to ../../results/, so a
# re-run never overwrites the deposited scRNA_myeloid_final.rds it can be compared with.
# =============================================================================

library(Seurat)
library(ggplot2)

OBJ      <- '../../objects/released/scRNA_final.rds'
OUT_DIR  <- '../../results/myeloid_states'
OUT_RDS  <- file.path(OUT_DIR, 'scRNA_myeloid_final.rds')

ANNO_COL  <- 'integrated_snn_res.1_anno'         # major cell type annotation
KEEP      <- c('myeloid', 'neutrophil')
MIN_CELLS <- 30                                  # samples with fewer myeloid cells are dropped
NFEATURES <- 2000

EXCLUDE_ORIG <- '6'          # tumour-mixed cluster in the first res = 0.3 clustering
SUBCLUSTER   <- '9'          # monocyte / dendritic cluster, split at res = 0.1

# res = 0.3 cluster (with the sub-clusters of cluster 9) -> fine state
FINE_MAP <- c('0' = 'CD163+ macrophage', '1' = 'CCL2+ inflammatory macrophage',
              '2' = 'CX3CR1+ microglia-like', '3' = 'C1QC+ microglia-like',
              '4' = 'GPNMB+ TAM', '5' = 'stress-activated myeloid',
              '6' = 'S100A8/A9+ neutrophil', '7' = 'IFN-activated myeloid',
              '8' = 'artifact',
              '9_0' = 'FCN1+ inflammatory monocyte', '9_1' = 'AREG+ LYZ+ myeloid',
              '9_2' = 'CD1C+ dendritic cell',
              '10' = 'proliferating myeloid', '11' = 'artifact')

# the same clusters -> coarse state used in the figures
COARSE_MAP <- c('0' = 'proinflammatory', '1' = 'proinflammatory', '2' = 'proinflammatory',
                '3' = 'microglia-like', '4' = 'proangiogenic', '5' = 'stress-related',
                '6' = 'neutrophil', '7' = 'interferon-activated', '8' = 'artifact',
                '9_0' = 'monocyte', '9_1' = 'proinflammatory', '9_2' = 'proinflammatory',
                '10' = 'proliferative', '11' = 'artifact')

# display order used in the paper figures
STATE_LEVELS <- c('proangiogenic', 'proinflammatory', 'neutrophil', 'interferon-activated',
                  'microglia-like', 'proliferative', 'monocyte', 'stress-related')
ARTIFACT <- 'artifact'

# 1. subset myeloid cells ------------------------------------------------------
obj <- readRDS(OBJ)
myeloid_cells <- colnames(obj)[!is.na(obj@meta.data[[ANNO_COL]]) &
                                 obj@meta.data[[ANNO_COL]] %in% KEEP &
                                 obj$in_cleaned_0806]
message('myeloid cells: ', length(myeloid_cells))
stc <- subset(obj, cells = myeloid_cells)
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
                                   normalization.method = 'LogNormalize', reduction = 'cca',
                                   k.filter = 10, dims = 1:30)
rm(stc.list); gc()
stc.combined <- IntegrateData(anchorset = anchors, k.weight = 10)
rm(anchors); gc()

DefaultAssay(stc.combined) <- 'integrated'
stc.combined <- ScaleData(stc.combined, verbose = FALSE)
stc.combined <- RunPCA(stc.combined, npcs = 30, seed.use = 42, verbose = FALSE)

# 3. drop the tumour-mixed cluster ---------------------------------------------
stc.combined <- FindNeighbors(stc.combined, reduction = 'pca', dims = 1:30, verbose = FALSE)
stc.combined <- FindClusters(stc.combined, resolution = 0.3, random.seed = 0,
                             cluster.name = 'integrated_snn_res.0.3_orig', verbose = FALSE)
print(table(stc.combined$integrated_snn_res.0.3_orig))
VlnPlot(stc.combined, features = c('PTPRC', 'AIF1', 'SOX2', 'PTPRZ1'),
        group.by = 'integrated_snn_res.0.3_orig', stack = TRUE, flip = TRUE)

stopifnot(EXCLUDE_ORIG %in% levels(stc.combined$integrated_snn_res.0.3_orig))
keep <- colnames(stc.combined)[as.character(stc.combined$integrated_snn_res.0.3_orig) != EXCLUDE_ORIG]
stc.kept <- subset(stc.combined, cells = keep)
rm(stc.combined); gc()
message('cells after removing the tumour-mixed cluster: ', ncol(stc.kept))

# 4. re-cluster and split the monocyte / dendritic cluster ---------------------
stc.kept <- FindNeighbors(stc.kept, reduction = 'pca', dims = 1:30, verbose = FALSE)
stc.kept <- FindClusters(stc.kept, resolution = 0.3, random.seed = 0,
                         cluster.name = 'integrated_snn_res.0.3_cleaned', verbose = FALSE)
stc.kept <- RunUMAP(stc.kept, reduction = 'pca', dims = 1:20, seed.use = 42, verbose = FALSE)
print(table(stc.kept$integrated_snn_res.0.3_cleaned))

Idents(stc.kept) <- 'integrated_snn_res.0.3_cleaned'
stc.kept <- FindSubCluster(stc.kept, SUBCLUSTER, graph.name = 'integrated_snn',
                           subcluster.name = 'sub.cluster', resolution = 0.1)
print(table(stc.kept$sub.cluster))

# 5. fine and coarse states, drop the artefact clusters ------------------------
clusters <- as.character(stc.kept$sub.cluster)
stopifnot(setequal(unique(clusters), names(FINE_MAP)))
stc.kept$myeloid_state_fine   <- unname(FINE_MAP[clusters])
stc.kept$myeloid_state_coarse <- unname(COARSE_MAP[clusters])

keep_final <- stc.kept$myeloid_state_coarse != ARTIFACT
stc.clean <- subset(stc.kept, cells = colnames(stc.kept)[keep_final])
rm(stc.kept); gc()
message('cells after removing artefact clusters: ', ncol(stc.clean))

stc.clean$myeloid_state_coarse <- factor(stc.clean$myeloid_state_coarse, levels = STATE_LEVELS)
stopifnot(!anyNA(stc.clean$myeloid_state_coarse))
print(table(stc.clean$myeloid_state_coarse))

DefaultAssay(stc.clean) <- 'RNA'
Idents(stc.clean) <- 'myeloid_state_coarse'
saveRDS(stc.clean, OUT_RDS)

# why the other myeloid cells carry no state; the all-cell object records this in
# substate_excluded_reason
reason <- setNames(rep(NA_character_, length(all_cells)), all_cells)
reason[!(orig_ident %in% keep_samples)] <- 'sample_lt30cells'
reason[setdiff(all_cells[orig_ident %in% keep_samples], keep)] <- 'tumor_mixed_res0.3_orig'
reason[setdiff(keep, colnames(stc.clean))] <- ARTIFACT
print(table(reason, useNA = 'ifany'))
message('wrote ', OUT_RDS)
