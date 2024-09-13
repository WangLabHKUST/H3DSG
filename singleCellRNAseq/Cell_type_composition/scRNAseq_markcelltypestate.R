setwd('~/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/')
sid = 'SP41'
load(paste0(sid,'/seurat.Rda'))
DimPlot(seurat, label = T)
gns = c('SOX2','PTPRZ1',  #tumor (stem) cells
        'PTPRC', #immune cells: microglia, macrophage, T, B, monocyte, NK cells
        'VWF','CLDN5', #endothellial
        'RGS5','BGN','TAGLN', #mural cells
        'CSF2RA','AIF1','C1QC','P2RY12','CCL3', #myeloid cells (microglia/macrophage/monocyte)
        'CD3D','CD3E','CD2', 'CD247',  #T/NK
        'BANK1','MS4A1', 'CD79A',#B
        'PLP1', 'MOG', #Oligodendrocytes
        'MKI67', #proliferating cells (tumor, normal)
        'SNAP25','MYT1L', #neurons,
        'CD33','FCGR3B' #nrutrophils
)#

VlnPlot(seurat,features = gns, stack = T, flip = T)

gns = read.delim('../rerun_cNMF/filteredmalignant/K27scDMG.cellular.states.txt')
seurat <- AddModuleScore(seurat,features = gns, name = 'P')

df= seurat@meta.data
df$cell.type = ifelse(df$seurat_clusters %in% c(8,17,18,20),'myeloid',
                      ifelse(df$seurat_clusters==14,'neutrophil',
                             ifelse(df$seurat_clusters %in% c(22),'B',
                             #ifelse(df$seurat_clusters==20,'oligodendrocyte',
                                    #ifelse(df$seurat_clusters==17,'unknown',
                                           ifelse(df$seurat_clusters %in% c(16),'T','malignant'))))#))
sts =c('G2M','OC-like','S','AC-like','MES-like','OPC-like')
df$cell.state = sts[apply(df[,c('P1','P2','P3','P4','P5','P6')],1,which.max)]
df$cell.score = apply(df[,c('P1','P2','P3','P4','P5','P6')],1,max)
df$cell.state=ifelse(df$cell.type=='malignant',df$cell.state, df$cell.type)
#df$cell.state[rownames(df) %in% rownames(seurat@reductions$umap@cell.embeddings)[seurat@reductions$umap@cell.embeddings[,2]<0 & seurat$seurat_clusters==16]]='OC-like'

seurat <-AddMetaData(seurat, metadata = df[,c('cell.state','cell.score')])
DimPlot(seurat, group.by  = c('cell.state'), label = T)
ggsave(filename =paste0(sid,'/',sid,".celltypeandstate.pdf"),width =10,height =7 )
FeaturePlot(seurat, features = c('cell.score'),max.cutoff = 1)


#
write.table(df[,c('cell.barcode','cell.type','cell.state','cell.score')],
            file = paste0(sid,'/',sid,".sc.celltype.txt"), row.names = F,quote = F, sep = "\t")
m = as.matrix(seurat@assays$RNA@counts)
write.table(m, file = paste0(sid,'/',sid,".sc.rawcount.txt"), quote = F, sep = "\t")
save(seurat, file = file.path(sid, "seurat.Rda"))

graphics.off(); rm(list=ls())
