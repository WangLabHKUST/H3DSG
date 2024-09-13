library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(glue)
library(DT)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(Seurat)


config = list(species = 'h_sapiens',min_cells = 3, min_features = 200,max_mito = NA,
              normalization = 'LogNormalize',var_regress = 'nCount_RNA,percent.mito',
              pcs_compute = 100, pcs_keep = 30, clustering_resolution = 1,seed = 42)

#0. specify the sample_id,  input and output folders
# 3' scRNAseq: SP04,SP05,SP11,SP13,SP15,SP22,SP24,SP41; n=8
# 5' scRNAseq: SP20, SP21; n=2
# snRNAseq: X172017, X177793, X747347,__; n=4 (another 1 only has spatial seq)
smplist = c('X172017', 'X177793', 'X747347','SP204') #'SP04','SP05','SP11','SP13','SP15','SP22','SP24','SP41',
for (smp in smplist){
  print(paste('==========',smp,'=========='))
  sample_id = smp
  cellranger_trio_dir = paste0('~/Dropbox/SpinalCordGlioma/scRNAseq/',smp,'.filtered_feature_bc_matrix/')
  
  output_dir = paste0('~/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/',smp)
  dir.create(output_dir,recursive = T)
  #1. read in data
  
  cellranger_matrix <- Seurat::Read10X(cellranger_trio_dir)
  
  #2. pre-filtering: cutoff on min cells and min features to reduce data volume
  seurat_prefilt <- CreateSeuratObject(
    counts       = cellranger_matrix,
    min.cells    = config$min_cells,    # specifies the min number of cells in which each gene must be detected
    min.features = config$min_features, # specifies the min number of genes which must be detected in each cell
    project      = sample_id            # populates the project.name slot in the Seurat object
  )
  
  seurat_prefilt@meta.data$sample.id    <- sample_id
  seurat_prefilt@meta.data$cell.barcode <- paste0(sample_id, "_", rownames(seurat_prefilt@meta.data))
  
  #3. compute mitochondria and ribosome RNA
  mito_genes <- grep("^MT-", rownames(GetAssayData(object = seurat_prefilt)), value = TRUE)
  ribo_genes <- grepl("^RPS|^RPL|^MRPS|^MRPL", rownames(GetAssayData(object = seurat_prefilt)))
  
  percent_mito <- Matrix::colSums(GetAssayData(object = seurat_prefilt)[mito_genes, ]) /
    Matrix::colSums(GetAssayData(object = seurat_prefilt)) * 100
  seurat_prefilt <- AddMetaData(seurat_prefilt, percent_mito, "percent.mito")
  percent_ribo <- Matrix::colSums(GetAssayData(object = seurat_prefilt)[ribo_genes, ]) /
    Matrix::colSums(GetAssayData(object = seurat_prefilt)) * 100
  seurat_prefilt <- AddMetaData(seurat_prefilt, percent_ribo, "percent.ribo")
  
  #4. QC and further filtering
  (thresholds <- data.frame(
    # te minimum number of features will be the greater of:
    # 400, or 2 standard deviations below the mean
    min_features = max(400, round(mean(seurat_prefilt@meta.data$nFeature_RNA) -
                                    2*sd(seurat_prefilt@meta.data$nFeature_RNA))),
    max_features = round(mean(seurat_prefilt@meta.data$nFeature_RNA) +
                           2*sd(seurat_prefilt@meta.data$nFeature_RNA)),
    min_mito     = 0,
    # by default,
    # the max mitochondrial content will be the maximum of:
    # 5%, or 2 standard deviations above the mean
    # the parameter config$max_mito allows to set a hard upper threshold,
    # which takes precedence
    max_mito     = ifelse(!is.na(config$max_mito),
                          config$max_mito,
                          max(5, round(mean(seurat_prefilt@meta.data$percent.mito) +
                                         2*sd(seurat_prefilt@meta.data$percent.mito))
                          )
    ),
    # set a max of 0 in case the value 2 standard deviations below the mean
    # is negative
    min_umi      = max(0, round(mean(seurat_prefilt@meta.data$nCount_RNA) -
                                  2*sd(seurat_prefilt@meta.data$nCount_RNA))),
    max_umi      = round(mean(seurat_prefilt@meta.data$nCount_RNA) +
                           2*sd(seurat_prefilt@meta.data$nCount_RNA))
  ))
  
  keep_cells <- WhichCells(
    object = seurat_prefilt,
    expression = nFeature_RNA>thresholds$min_features & 
      nFeature_RNA<thresholds$max_features & 
      nCount_RNA>thresholds$min_umi &
      nCount_RNA<thresholds$max_umi &
      percent.mito>thresholds$min_mito &
      percent.mito<thresholds$max_mito
  )
  
  
  plot_grid(VlnPlot(seurat_prefilt, c("nFeature_RNA"), pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat_prefilt, c("nCount_RNA"), pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat_prefilt, c("percent.mito"), pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat_prefilt, c("percent.ribo"), pt.size = -1) +
              theme(legend.position = "none"),
            ncol = 4)
  ggsave(filename =file.path(output_dir,"prefiltering.nfeaturencountpmitopribo.pdf"),width = 8,height = 4.5 )
  
  
  seurat <- subset(seurat_prefilt, cells = keep_cells)
  print('Before filtering:');seurat_prefilt
  print('After filtering:');seurat
  
  
  
  #5. generate a summary of filtering metrics
  filtering_criteria <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo")
  
  # compute summary stats for each metric
  plot_grid(VlnPlot(seurat, c("nFeature_RNA"), pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat, c("nCount_RNA"),   pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat, c("percent.mito"), pt.size = -1) +
              theme(legend.position = "none"),
            VlnPlot(seurat, c("percent.ribo"), pt.size = -1) +
              theme(legend.position = "none"),
            ncol = 4)
  ggsave(filename =file.path(output_dir,"postfiltering.nfeaturencountpmitopribo.pdf"),width = 8,height = 4.5 )
  
  filtering_metrics <- sapply(filtering_criteria, function(criterion) {
    
    min_pre   <- round(min(seurat_prefilt@meta.data  %>% pull(criterion)), 2)
    mean_pre  <- mean(seurat_prefilt@meta.data %>% pull(criterion))
    max_pre   <- max(seurat_prefilt@meta.data  %>% pull(criterion))
    sd_pre    <- sd(seurat_prefilt@meta.data   %>% pull(criterion))
    
    min_post  <- min(seurat@meta.data  %>% pull(criterion))
    mean_post <- mean(seurat@meta.data %>% pull(criterion))
    max_post  <- max(seurat@meta.data  %>% pull(criterion))
    sd_post   <- sd(seurat@meta.data   %>% pull(criterion))
    
    return(c("min.preQC"   = min_pre,
             "mean.preQC"  = mean_pre,
             "max.preQC"   = max_pre,
             "sd.preQC"    = sd_pre,
             "min.postQC"  = min_post,
             "mean.postQC" = mean_post,
             "max.postQC"  = max_post,
             "sd.postQC"   = sd_post))
    
  })
  
  # round to 2 decimal places
  filtering_metrics <- apply(filtering_metrics, 2, round, 2)
  
  # transform into a dataframe
  filtering_metrics <- filtering_metrics %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "criterion") %>%
    dplyr::select(criterion, min.preQC, min.postQC, max.preQC, max.postQC, mean.preQC, mean.postQC, sd.preQC, sd.postQC)
  
  # add thresholds
  filtering_metrics$min.threshold <- c(thresholds$min_features,
                                       thresholds$min_umi,
                                       thresholds$min_mito,
                                       # No min threshold used for % ribo
                                       NA)
  
  filtering_metrics$max.threshold <- c(thresholds$max_features,
                                       thresholds$max_umi,
                                       thresholds$max_mito,
                                       # No max threshold used for % ribo
                                       NA)
  filtering_metrics
  
  # compute number of cells before and after filtering
  N_cells_metrics <- data.frame(
    "N_cells_before" = dim(seurat_prefilt@meta.data)[1],
    "N_cells_after"  = dim(seurat@meta.data)[1]) %>%
    mutate(Prop_kept = round(N_cells_after / N_cells_before, 2))
  
  N_cells_metrics
  
  
  #6. Normalization and scaling
  # normalization 1: scale counts to 10000 UMIs per cell, and log2-transform the counts
  # set the default assay to RNA to run the log-normalization & scaling
  seurat <- seurat %>% 
    NormalizeData(normalization.method = "LogNormalize",
                  scale.factor         = 10000) %>% 
    # identify variable genes
    FindVariableFeatures(mean.function       = ExpMean,
                         dispersion.function = LogVMR) %>%
    # regress out variables which are sources of unwanted variation, and z-score data
    ScaleData(vars.to.regress = unlist(str_split(config$var_regress, ",")))
  
  # normalization 2: using SCTransform
  # this command also identifies variable features and produces scaled values
  seurat <- seurat %>% SCTransform(vars.to.regress = unlist(str_split(config$var_regress, ",")), verbose = FALSE)
  
  # choose the normalization method to use for downstream analysis
  DefaultAssay(object = seurat) <- switch(config$normalization,"LogNormalize" = "RNA","SCTransform"  = "SCT")
  
  # confirm the default assay set
  DefaultAssay(seurat)
  
  #7. Dimensionality reduction and clustering
  seurat <- seurat %>%
    # compute PCA, based on scaled data
    RunPCA(pc.genes        = VariableFeatures(.),
           npcs            = config$pcs_compute,
           ndims.print     = 1:5,
           nfeatures.print = 5) %>%
    # compute tSNE embedding, based on retained PCs
    RunTSNE(dims = 1:config$pcs_keep, verbose = FALSE, seed.use = config$seed) %>%
    # compute UMAP embedding, based on retained PCs
    RunUMAP(dims = 1:config$pcs_keep, verbose = FALSE, seed.use = config$seed)
  
  ElbowPlot(seurat, ndims = config$pcs_compute)
  ggsave(filename =file.path(output_dir,"PCA.elbowplot.pdf"),width = 4.5,height = 3.4 )
  
  seurat <- FindNeighbors(seurat,reduction = "pca",dims = 1:config$pcs_keep,verbose   = TRUE,nn.eps    = 0.5)
  
  clustering_fun <- purrr::partial(FindClusters,
                                   seurat,
                                   verbose = FALSE,
                                   n.start = 10,
                                   random.seed = config$seed)
  
  seurat <- clustering_fun(resolution = 0.6)
  
  p1 <- DimPlot(object = seurat, reduction = "tsne") +
    ggtitle("res 0.6")
  
  seurat <- clustering_fun(resolution = 0.8)
  
  p2 <- DimPlot(object = seurat, reduction = "tsne") +
    ggtitle("res 0.8")
  
  seurat <- clustering_fun(resolution = 1)
  
  p3 <- DimPlot(object = seurat, reduction = "tsne") +
    ggtitle("res 1")
  
  seurat <- clustering_fun(resolution = 2)
  
  p4 <- DimPlot(object = seurat, reduction = "tsne") +
    ggtitle("res 2")
  
  seurat <- clustering_fun(resolution = config$clustering_resolution)
  
  p5 <- DimPlot(object = seurat, reduction = "tsne") +
    ggtitle(paste0("chosen resolution: ", config$clustering_resolution))
  
  plot_grid(p1, p2, p3, p4, p5, ncol = 3)
  ggsave(filename =file.path(output_dir,"DimensionReduction.tSNE.pdf"),width =17,height = 8.5 )
  Idents(object = seurat) <- paste0(switch(config$normalization,
                                           "LogNormalize" = "RNA",
                                           "SCTransform"  = "SCT"),
                                    "_snn_res.",
                                    config$clustering_resolution)
  #8. Cell cycle scoring
  s_genes   <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  
  seurat <- CellCycleScoring(seurat,s.features   = s_genes,g2m.features = g2m_genes)
  
  cc.markers <- switch(config$species,
                       "h_sapiens" = c("PCNA", "TOP2A", "MCM6", "MKI67"),
                       "m_musculus" = c("Pcna", "Top2a", "Mcm6", "Mki67"))
  
  RidgePlot(seurat, group.by = "Phase", features = cc.markers, ncol = 2)
  ggsave(filename =file.path(output_dir,"CellCyle.phaseMarkers.pdf"),width =6.8,height =5.8 )
  #9. Cluster-level QC
  table(Idents(object = seurat))
  
  p0 <- DimPlot(object = seurat, reduction = "pca", cols = seurat@misc$colours)
  p1 <- DimPlot(object = seurat, reduction = "tsne", cols = seurat@misc$colours,label = T)
  p2 <- FeaturePlot(object = seurat, reduction = "tsne", features = "nFeature_RNA") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$nFeature_RNA))
  p3 <- FeaturePlot(object = seurat, reduction = "tsne", features = "nCount_RNA") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$nCount_RNA))
  p4 <- FeaturePlot(object = seurat, reduction = "tsne", features = "percent.mito") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$percent.mito))
  p5 <- FeaturePlot(object = seurat, reduction = "tsne", features = "percent.ribo") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$percent.ribo))
  
  plot_grid(p0, p1, p2, p3, p4, p5,ncol = 3)
  ggsave(filename =file.path(output_dir,"QCfeature.tSNE.pdf"),width =17,height = 8.5 )
  
  p0 <- DimPlot(object = seurat, reduction = "pca", cols = seurat@misc$colours)
  p1 <- DimPlot(object = seurat, reduction = "umap", cols = seurat@misc$colours,label = T)
  p2 <- FeaturePlot(object = seurat, reduction = "umap", features = "nFeature_RNA") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$nFeature_RNA))
  p3 <- FeaturePlot(object = seurat, reduction = "umap", features = "nCount_RNA") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$nCount_RNA))
  p4 <- FeaturePlot(object = seurat, reduction = "umap", features = "percent.mito") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$percent.mito))
  p5 <- FeaturePlot(object = seurat, reduction = "umap", features = "percent.ribo") +
    scale_color_gradient2(low = "darkblue", mid = "lightgrey", high = "red",
                          midpoint = mean(seurat@meta.data$percent.ribo))
  
  plot_grid(p0, p1, p2, p3, p4, p5,ncol = 3)
  ggsave(filename =file.path(output_dir,"QCfeature.UMAP.pdf"),width =17,height = 8.5 )
  
  clust_df <- data.frame(table(Idents(object = seurat)))
  clust_df[1, 2] <- paste0("N=", clust_df[1, 2])
  colnames(clust_df) <- c("ident", "N")
  
  vln_fun <- function(criterion) {
    
    # plot the labels at a value slightly below the max, to ensure they're shown
    # within plot limits
    clust_df$y <- max(seurat[[criterion]]) * 0.95
    Seurat::VlnPlot(seurat, criterion, cols = seurat@misc$colours, pt.size = 0.5) +
      theme(legend.position = "none") +
      geom_text(data = clust_df, aes(label = N, x = ident, y = y))
    
  }
  
  vln_fun("nFeature_RNA")
  ggsave(filename =file.path(output_dir,"QCfeature.nFeature_RNA.violin.pdf"),width =17,height = 4.5 )
  vln_fun("nCount_RNA")
  ggsave(filename =file.path(output_dir,"QCfeature.nCount_RNA.violin.pdf"),width =17,height = 4.5 )
  vln_fun("percent.mito")
  ggsave(filename =file.path(output_dir,"QCfeature.percent.mito.violin.pdf"),width =17,height = 4.5 )
  vln_fun("percent.ribo")
  ggsave(filename =file.path(output_dir,"QCfeature.percent.ribo.violin.pdf"),width =17,height = 4.5 )
  
  #10. Identify cluster markers
  cluster_markers <- FindAllMarkers(object = seurat, verbose = FALSE)
  
  # display the top 10 per cluster
  
  top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  DoHeatmap(seurat, features = top10$gene, group.colors = seurat@misc$colours) +
    NoLegend() +
    scale_fill_gradientn(colors = c("#2166AC", "#E5E0DC", "#B2182B"))
  ggsave(filename =file.path(output_dir,"Cluster.top10marker.heatmap.pdf"),width =17,height = 17 )
  #11. save to output folders
  write_tsv(cluster_markers, file.path(output_dir, "cluster_markers.tsv"))
  seurat@misc$params$config     <- config
  seurat@misc$filtering_metrics <- filtering_metrics
  seurat@misc$n_cells           <- N_cells_metrics
  
  # write metrics/thresholds to file
  filtering_metrics_out <- filtering_metrics %>%
    gather("metrics_name", "value", -criterion) %>%
    unite("metrics", criterion:metrics_name, sep = "_") %>%
    spread("metrics", "value") %>% 
    bind_cols(N_cells_metrics)
  
  # include the short version of the most recent git repository SHA
  write.table(as.data.frame(t(filtering_metrics_out)), file = paste0(output_dir, "/seurat_metrics.tsv"),col.names = F, quote = F, sep = "\t")
  save(seurat, file = file.path(output_dir, "seurat.Rda"))
}

