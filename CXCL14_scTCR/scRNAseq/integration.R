setwd('~/Documents/Projects/H3DSG/Revision/CXCL14_scTCR/reanalyze/scRNAseq/')
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
gcc = c('TCell', 'TCXCL14', 'TNET') 

fls = list.files(path = ".", pattern = 'seurat.Rda',recursive = T, full.names = T)
info_samples = data.frame(Path = fls, Sample = sapply(fls, function(x) strsplit(x, "\\/")[[1]][2]))

info_samples = info_samples[info_samples$Sample %in% gcc,]
seurat_indiv <- list()
for (row in 1:nrow(info_samples)) { 
  seurat <- get(load(info_samples[row, ]$Path))
  
  seurat@project.name <- info_samples[row, ]$Sample
  seurat$Sample <- seurat@project.name
  print(seurat@project.name)
  
  seurat_indiv[[row]] <- seurat
  rm(seurat)
}

# Preprocess data ------------------------------------------------------------
message("@ preprocessing data...")

# merge into a single Seurat object, normalize, scale, and run PCA
seurat_joint <- merge(x = seurat_indiv[[1]],
                      y = seurat_indiv[2:length(seurat_indiv)],
                      merge.data = FALSE)

# clean up
rm(seurat_indiv)

seurat_joint <- seurat_joint %>% 
  # all samples should have been normalized the same way, but just in case, 
  # re-run it here
  Seurat::NormalizeData() %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = .@var.genes, npcs = 30) %>% 
  #RunTSNE(dims = 1:30,  seed.use = 100, check_duplicates = FALSE) %>%
  RunUMAP(dims = 1:30,  seed.use = 100)

# perform clustering
seurat_joint <- seurat_joint %>% 
  FindNeighbors(seurat,
                reduction = "pca",
                dims      = 1:30,
                verbose   = TRUE,
                nn.eps    = 0.5) %>% 
  FindClusters(n.start = 10,
               random.seed = 100, 
               resolution = 0.5)

# Save -----------------------------------------------------------------------

message("@ saving joined data...")

save(seurat_joint, file = "seurat_joint.Rda")
ggsave(plot = VlnPlot(seurat_joint,features = c('nFeature_RNA','nCount_RNA','percent.mito','percent.ribo'),group.by = 'orig.ident',stack = T,flip = T),
       filename = "seurat_QC_bySample.pdf", width = 7, height = 4)

library(scCustomize)
Stacked_VlnPlot(seurat_object = seurat_joint_harmony, features = c('nFeature_RNA','nCount_RNA','percent.mito'), 
                group.by = 'orig.ident',x_lab_rotate = F,colors_use = c("dodgerblue", "gold", "orange"))

#Integration by Harmony
library(harmony)
png(filename = "harmony.convergence.png", width = 500, height = 400)
library(Seurat)
seurat_joint_harmony <- seurat_joint %>% 
  RunHarmony(group.by.vars    = 'orig.ident',
             assay.use = 'SCT',
             #reduction        = "pca",
             #dims.use         = 1:30,
             plot_convergence = TRUE)
dev.off()

# clean up
rm(seurat_joint)

message("@ performing downstream analysis...")

# run tSNE, UMAP, clustering
seurat_joint_harmony <- seurat_joint_harmony %>% 
  RunTSNE(dims = 1:30, verbose = T, seed.use = 100, reduction = "harmony") %>%
  RunUMAP(dims = 1:30, verbose = T, seed.use = 100, reduction = "harmony")

seurat_joint_harmony <- seurat_joint_harmony %>% 
  FindNeighbors(seurat,
                reduction = "harmony",
                dims      = 1:30,
                verbose   = TRUE,
                nn.eps    = 0.5) %>% 
  FindClusters(verbose = T,
               n.start = 10,
               random.seed = 100, 
               resolution = 0.3)

# Save ----

message("@ saving integrated data...")

seurat_joint_harmony_dr <- list(
  "pca"  = seurat_joint_harmony@reductions$pca@cell.embeddings[, c(1, 2)],
  "tsne" = seurat_joint_harmony@reductions$tsne@cell.embeddings[, c(1, 2)],
  "umap" = seurat_joint_harmony@reductions$umap@cell.embeddings[, c(1, 2)])
saveRDS(seurat_joint_harmony_dr, file = "dimred.harmony.Rds")

save(seurat_joint_harmony, file = "seurat_joint.harmony.Rda")
ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "umap",split.by = 'orig.ident'),
       filename = "UMAP_clusters.harmony.pdf", width = 6.5, height = 2.5)

#examine some genes
mks_all = seurat_joint_harmony %>% JoinLayers() %>% FindAllMarkers()
mks_all %>% group_by(cluster) %>%dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%ungroup() -> mks_top10
DotPlot(seurat_joint_harmony, mks_top10$gene) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x = '', y = '', fill = 'avg.exp', size = 'pct.exp')

gns = c('CD3D','CD8A','PDCD1','CTLA4','FOXP3','GNLY','SELL','CXCL10','HIST1H3C','MAF','KLRB1',
        'MCM2','TOP2A','UBE2C','GZMB','GZMA','TNFAIP3','CCR7','IL32')
DotPlot(seurat_joint_harmony, gns) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(x = '', y = '', fill = 'avg.exp', size = 'pct.exp')



