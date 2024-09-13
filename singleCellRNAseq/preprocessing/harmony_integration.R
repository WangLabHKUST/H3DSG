setwd('~/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/')
library(Seurat)
library(ggplot2)
gcc = c( "SP04" ,"SP15" , "SP41" ,"SP204" ,"SP202","SP24")
gnf1 = c("SP05","SP21","SP205","SP11","SP203","SP13","SP201","SP22","SP20")
#X172017 = SP202; X177793 = SP203, X747347 = SP205
fls = list.files(path = ".", pattern = 'seurat.Rda',recursive = T, full.names = T)
info_samples = data.frame(Path = fls, Sample = sapply(fls, function(x) strsplit(x, "\\/")[[1]][2]))
info_samples$Sample[info_samples$Sample=='X172017'] = 'SP202'
info_samples$Sample[info_samples$Sample=='X177793'] = 'SP203'
info_samples$Sample[info_samples$Sample=='X747347'] = 'SP205'
info_samples$Sample[info_samples$Sample=='SP20T5'] = 'SP20';
info_samples$Sample[info_samples$Sample=='SP21T5'] = 'SP21'


seurat_indiv <- list()
for (row in 1:12) { #ignore the three snRNAseq
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
dir.create('integration')
seurat_joint_dr <- list(
  "pca"  = seurat_joint@reductions$pca@cell.embeddings[, c(1, 2)],
  #"tsne" = seurat_joint@reductions$tsne@cell.embeddings[, c(1, 2)],
  "umap" = seurat_joint@reductions$umap@cell.embeddings[, c(1, 2)])
saveRDS(seurat_joint_dr, file = "integration/dimred.Rds")

seurat_joint_meta <- seurat_joint@meta.data
saveRDS(seurat_joint_meta, file = "integration/metadata.Rds")

save(seurat_joint, file = "integration/seurat_joint.Rda")
ggsave(plot = DimPlot(seurat_joint, reduction = "umap", label =T), filename = "integration/joint_UMAP_clusters.pdf", width = 10, height = 8)


#Integration by Harmony
library(harmony)
png(filename = "integration/harmony.convergence.png", width = 500, height = 400)
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
  RunTSNE(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony") %>%
  RunUMAP(dims = 1:info_experiment$n_pcs, verbose = info_experiment$verbose, seed.use = 100, reduction = "harmony")

seurat_joint_harmony <- seurat_joint_harmony %>% 
  FindNeighbors(seurat,
                reduction = "harmony",
                dims      = 1:info_experiment$n_pcs,
                verbose   = TRUE,
                nn.eps    = 0.5) %>% 
  FindClusters(verbose = info_experiment$verbose,
               n.start = 10,
               random.seed = 100, 
               resolution = 0.5)

# Save ----

message("@ saving integrated data...")

seurat_joint_harmony_dr <- list(
  "pca"  = seurat_joint_harmony@reductions$pca@cell.embeddings[, c(1, 2)],
  "tsne" = seurat_joint_harmony@reductions$tsne@cell.embeddings[, c(1, 2)],
  "umap" = seurat_joint_harmony@reductions$umap@cell.embeddings[, c(1, 2)])
saveRDS(seurat_joint_harmony_dr, file = "output/dimred.harmony.Rds")

save(seurat_joint_harmony, file = "output/seurat_joint.harmony.Rda")
ggsave(plot = DimPlot(seurat_joint_harmony, reduction = "umap"), filename = "integration/UMAP_clusters.harmony.pdf", width = 10, height = 8)



