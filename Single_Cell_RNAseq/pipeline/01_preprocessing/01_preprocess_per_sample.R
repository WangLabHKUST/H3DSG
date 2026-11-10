# =============================================================================
# Per-sample preprocessing of the H3-DSG single cell RNA-seq cohort
#
# Input : CellRanger filtered matrices in ../../inputs/raw_counts/<file_id>/, one
#         directory per sample (see ../../sample_metadata.tsv, column file_id)
# Output: ../../results/per_sample/<sid>_minf400.rds, one processed Seurat object per
#         sample. These per-sample objects are the input of
#         02_integration/01_cca_integration.R and of the inferCNV analysis; the
#         deposited copies are in ../../objects/per_sample/.
#
# Paths are relative to this directory: run the script from here.
#
# Per sample: read counts -> drop genes in < 3 cells and cells with < 200 genes ->
# mitochondrial / ribosomal content -> drop unnamed ENSG features -> per-sample
# adaptive QC filter -> normalise, cluster, SCTransform.
#
# The QC filter keeps cells within mean +/- 2 SD of UMI count and gene count, with a
# floor of 400 genes, and below an adaptive mitochondrial percentage capped to 5-12%.
#
# Runtime about 10 min per sample. Seurat 5.2.1, R 4.5.1.
# =============================================================================

library(Seurat)
library(Matrix)

COUNTS_DIR <- '../../inputs/raw_counts'
OUT_DIR    <- '../../results/per_sample'
SAMPLES    <- read.delim('../../sample_metadata.tsv', colClasses = 'character')

MIN_CELLS    <- 3           # a gene must be detected in at least this many cells
MIN_FEATURES <- 200         # a cell must have at least this many genes
MIN_FEAT_FLOOR <- 400       # floor for the adaptive gene count cutoff
MITO_PATTERN <- '^MT-'
RIBO_PATTERN <- '^RPS|^RPL|^MRPS|^MRPL'
ENSG_PATTERN <- '^ENSG[0-9]{11}'   # features without a gene symbol

# adaptive per-sample cutoffs, mean +/- 2 SD
compute_thresholds <- function(meta, min_feats = MIN_FEAT_FLOOR) {
  list(min_umi      = max(0, round(mean(meta$nCount_RNA) - 2 * sd(meta$nCount_RNA))),
       max_umi      =        round(mean(meta$nCount_RNA) + 2 * sd(meta$nCount_RNA)),
       min_features = max(min_feats, round(mean(meta$nFeature_RNA) - 2 * sd(meta$nFeature_RNA))),
       max_features =        round(mean(meta$nFeature_RNA) + 2 * sd(meta$nFeature_RNA)),
       max_mito     = min(12, max(5, round(mean(meta$percent.mito) + 2 * sd(meta$percent.mito)))))
}

for (i in seq_len(nrow(SAMPLES))) {
  sid     <- SAMPLES$sid[i]
  file_id <- SAMPLES$file_id[i]
  message('========== ', sid, ' (', file_id, ') ==========')

  # 1. read counts
  mat <- Read10X(file.path(COUNTS_DIR, file_id))
  seurat <- CreateSeuratObject(counts = mat, min.cells = MIN_CELLS,
                               min.features = MIN_FEATURES, project = sid)
  seurat$file_id <- file_id

  # 2. mitochondrial and ribosomal content
  cm <- GetAssayData(seurat, layer = 'counts')
  cs <- Matrix::colSums(cm)
  seurat$percent.mito <- Matrix::colSums(cm[grep(MITO_PATTERN, rownames(cm)), , drop = FALSE]) / cs * 100
  seurat$percent.ribo <- Matrix::colSums(cm[grepl(RIBO_PATTERN, rownames(cm)), , drop = FALSE]) / cs * 100

  # 3. drop features that have no gene symbol.
  #    nCount_RNA / nFeature_RNA are deliberately NOT recomputed here: the QC cutoffs
  #    below were derived on the pre-ENSG counts, and recomputing them first shifts
  #    the mean / SD enough to move a few cells across the threshold.
  seurat <- subset(seurat, features = grep(ENSG_PATTERN, rownames(seurat), value = TRUE, invert = TRUE))

  # 4. QC filter
  th <- compute_thresholds(seurat@meta.data)
  print(unlist(th))
  keep <- with(seurat@meta.data,
               nFeature_RNA > th$min_features & nFeature_RNA < th$max_features &
               nCount_RNA   > th$min_umi      & nCount_RNA   < th$max_umi &
               percent.mito < th$max_mito)
  message('cells: ', ncol(seurat), ' -> ', sum(keep))
  seurat <- subset(seurat, cells = colnames(seurat)[keep])

  # counts of the retained cells, now over the gene-symbol features only
  cm <- GetAssayData(seurat, layer = 'counts')
  seurat$nCount_RNA   <- Matrix::colSums(cm)
  seurat$nFeature_RNA <- Matrix::colSums(cm > 0)

  # 5. normalise, cluster, SCTransform
  seurat <- NormalizeData(seurat, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
  seurat <- ScaleData(seurat, vars.to.regress = c('nCount_RNA', 'percent.mito'), verbose = FALSE)
  seurat <- RunPCA(seurat, npcs = 100, seed.use = 42, verbose = FALSE)
  seurat <- RunTSNE(seurat, dims = 1:30, seed.use = 42)
  seurat <- RunUMAP(seurat, dims = 1:30, seed.use = 42, verbose = FALSE)
  seurat <- FindNeighbors(seurat, reduction = 'pca', dims = 1:30, nn.eps = 0.5, verbose = FALSE)
  seurat <- FindClusters(seurat, resolution = 0.5, random.seed = 42, n.start = 10, verbose = FALSE)
  seurat <- SCTransform(seurat, vars.to.regress = c('percent.mito', 'percent.ribo'), verbose = FALSE)
  DefaultAssay(seurat) <- 'RNA'

  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(seurat, file.path(OUT_DIR, paste0(sid, '_minf400.rds')))
  rm(seurat, mat, cm); gc()
}
message('done: ', nrow(SAMPLES), ' samples written to ', OUT_DIR)
