# =============================================================================
# CCA integration of the 19 H3-DSG samples and major cell type annotation
#
# Input : ../../objects/per_sample/<sid>_minf400.rds, the deposited per-sample objects
#         (01_preprocessing/01_preprocess_per_sample.R writes the same files to
#         ../../results/per_sample/; point PREP_DIR there to chain your own run)
# Output: ../../results/integration/scRNA_integrated.rds, the integrated object with
#           - assays RNA and integrated
#           - reductions pca, umap (integrated), umap_uncorrected, umap_pc30
#           - integrated_snn_res.1 and the major cell type column
#             integrated_snn_res.1_anno (tumor, myeloid, neutrophil, T/NK, B,
#             oligo, vascular)
#
# The released object scRNA_final.rds is this object with the tumour, T/NK and
# myeloid state columns of the downstream analyses merged in.
#
# Runtime about 14 h on 40 cores, dominated by FindIntegrationAnchors.
# Seurat 5.2.1, R 4.5.1.
#
# Paths are relative to this directory: run the script from here.
# =============================================================================

library(Seurat)
library(ggplot2)
library(future)
options(future.globals.maxSize = 16 * 1024^3)

PREP_DIR <- '../../objects/per_sample'
OUT_DIR  <- '../../results/integration'
OUT_RDS  <- file.path(OUT_DIR, 'scRNA_integrated.rds')
SAMPLES  <- read.delim('../../sample_metadata.tsv', colClasses = 'character')
# Merge in sorted sid order, the order of the released object. Sample order affects the integration, and
# major_cluster_annotation.csv labels clusters by number, so another order could mislabel cell types.
SAMPLES  <- SAMPLES[order(SAMPLES$sid), ]
ANNO_MAP <- read.csv('major_cluster_annotation.csv', colClasses = 'character')

N_PC      <- 30
NFEATURES <- 2000

# 1. merge the per-sample objects ----------------------------------------------
# cell names get the sample id as prefix, e.g. PL_1_AAACCCAAGAAGTATC-1
obj_list <- lapply(SAMPLES$sid, function(sid) readRDS(file.path(PREP_DIR, paste0(sid, '_minf400.rds'))))
names(obj_list) <- SAMPLES$sid
obj <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids = SAMPLES$sid)
DefaultAssay(obj) <- 'RNA'
for (a in setdiff(Assays(obj), 'RNA')) obj[[a]] <- NULL
obj@reductions <- list()
message('merged: ', ncol(obj), ' cells from ', length(obj_list), ' samples')

obj$Group          <- SAMPLES$group[match(as.character(obj$orig.ident), SAMPLES$sid)]
obj$Tumor.Location <- SAMPLES$tumor_location[match(as.character(obj$orig.ident), SAMPLES$sid)]

# 2. uncorrected embedding, kept for comparison --------------------------------
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = 'vst', nfeatures = NFEATURES, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = N_PC, seed.use = 42, verbose = FALSE)
obj <- RunUMAP(obj, reduction = 'pca', dims = 1:N_PC, seed.use = 42,
               reduction.name = 'umap_uncorrected', verbose = FALSE)
DimPlot(obj, reduction = 'umap_uncorrected', group.by = 'orig.ident')

# 3. per-sample CCA integration ------------------------------------------------
obj_list <- SplitObject(obj, split.by = 'orig.ident')
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, selection.method = 'vst', nfeatures = NFEATURES, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = NFEATURES)

anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features,
                                  reduction = 'cca', dims = 1:N_PC, k.filter = 200)
rm(obj_list); gc()
obj_cca <- IntegrateData(anchorset = anchors, dims = 1:N_PC, k.weight = 100)
rm(anchors); gc()

DefaultAssay(obj_cca) <- 'integrated'
obj_cca <- ScaleData(obj_cca, verbose = FALSE)
obj_cca <- RunPCA(obj_cca, npcs = N_PC, seed.use = 42, verbose = FALSE)
obj_cca <- RunUMAP(obj_cca, reduction = 'pca', dims = 1:N_PC, seed.use = 42, verbose = FALSE)

# 4. major cell type annotation ------------------------------------------------
# a second UMAP with a large neighbourhood is used for the figures
obj_cca <- FindNeighbors(obj_cca, reduction = 'pca', dims = 1:N_PC, verbose = FALSE)
obj_cca <- RunUMAP(obj_cca, reduction = 'pca', dims = 1:N_PC, n.neighbors = 150, min.dist = 0.1,
                   seed.use = 42, reduction.name = 'umap_pc30', verbose = FALSE)
obj_cca <- FindClusters(obj_cca, resolution = 1, random.seed = 0, verbose = FALSE)

# cluster -> major cell type, assigned from canonical markers (see markers below)
res_col <- 'integrated_snn_res.1'
stopifnot(setequal(levels(obj_cca@meta.data[[res_col]]), ANNO_MAP$cluster))
obj_cca$integrated_snn_res.1_anno <- ANNO_MAP$anno[match(as.character(obj_cca@meta.data[[res_col]]),
                                                         ANNO_MAP$cluster)]
print(table(obj_cca$integrated_snn_res.1_anno))

markers <- c('SOX2', 'PTPRZ1', 'OLIG2',            # malignant
             'PTPRC',                              # immune
             'AIF1', 'C1QC', 'P2RY12',             # myeloid
             'CSF3R', 'FCGR3B', 'S100A9',          # neutrophil
             'CD3D', 'IL32', 'NKG7',               # T / NK
             'CD79A', 'MS4A1',                     # B
             'PLP1', 'MAG',                        # oligodendrocyte
             'VWF', 'CLDN5')                       # vascular
DotPlot(obj_cca, features = markers, assay = 'RNA', group.by = 'integrated_snn_res.1_anno') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic'))
DimPlot(obj_cca, reduction = 'umap_pc30', group.by = 'integrated_snn_res.1_anno', label = TRUE)

DefaultAssay(obj_cca) <- 'RNA'
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(obj_cca, OUT_RDS)
message('wrote ', OUT_RDS)
