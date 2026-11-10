# =============================================================================
# Quality control overview of the integrated cohort
#
# Input : ../../objects/released/scRNA_final.rds
# Output: QC violin plots per sample and the marker panel used to name the major
#         cell types in 02_integration/01_cca_integration.R
#
# Paths are relative to this directory: run the script from here.
# =============================================================================

library(Seurat)
library(ggplot2)

OBJ     <- '../../objects/released/scRNA_final.rds'
SAMPLES <- read.delim('../../sample_metadata.tsv', colClasses = 'character')

so <- readRDS(OBJ)
mt <- so@meta.data
mt$sample <- factor(SAMPLES$renamed_id[match(as.character(mt$orig.ident), SAMPLES$sid)],
                    levels = SAMPLES$renamed_id[order(SAMPLES$group, SAMPLES$sid)])

qc_violin <- function(y, log_scale = TRUE) {
  p <- ggplot(mt, aes(x = sample, y = .data[[y]])) +
    geom_violin(aes(fill = sample), show.legend = FALSE) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (log_scale) p + scale_y_log10() else p
}
qc_violin('nFeature_RNA')
qc_violin('nCount_RNA')
qc_violin('percent.mito', log_scale = FALSE)

# cells per sample after filtering
print(table(mt$sample))

# markers behind the major cell type labels
gns <- c('PTPRZ1', 'SOX2',            # malignant
         'PTPRC',                     # immune
         'AIF1', 'CD14',              # myeloid
         'FCGR3B', 'S100A9', 'CSF3R', # neutrophil
         'CD3D', 'IL32', 'NKG7',      # T / NK
         'CD79A', 'MS4A1',            # B
         'VWF',                       # vascular
         'PLP1', 'MAG')               # oligodendrocyte
DimPlot(so, reduction = 'umap_pc30', group.by = 'integrated_snn_res.1_anno', label = TRUE)
FeaturePlot(so, features = gns, reduction = 'umap_pc30', order = TRUE)
