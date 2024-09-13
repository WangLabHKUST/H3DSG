

library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scg_myeloid//')
load('~/Dropbox/communter/seurat_joint.Rda')
stc = subset(seurat_joint, subset = cell.state %in% c('myeloid'))
rm(seurat_joint); gc()

table(stc$orig.ident) #note most samples have a pretty large number of cells
#start integration
stc.list <- SplitObject(stc, split.by = "orig.ident") #split by source
stc.list <- lapply(X = stc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = stc.list, nfeatures = 3000)
stc.list <- PrepSCTIntegration(object.list = stc.list, anchor.features = features)

#find integration anchors and perfrom integration
stc.anchors <- FindIntegrationAnchors(object.list = stc.list, anchor.features = features, normalization.method = 'SCT')
stc.combined <- IntegrateData(anchorset = stc.anchors,k.weight = 10,normalization.method = "SCT")

stc.combined <- ScaleData(stc.combined, verbose = FALSE) %>%
  RunPCA( npcs = 30, verbose = FALSE) %>%
  RunUMAP( reduction = "pca", dims = 1:30) %>%
  RunTSNE( reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

#DefaultAssay(stc.combined) <- "integrated"
p1 <- DimPlot(stc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(stc.combined, reduction = "umap",  label = TRUE,repel = TRUE)
p1 + p2
#ggsave(file = 'Integrated.UMAP.pdf', width = 8, height = 3.5)
p1 <- DimPlot(stc.combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(stc.combined, reduction = "tsne",  label = TRUE,repel = TRUE)
p1 + p2

# identify markers for each subgroup

DimPlot(stc.combined, label = T,reduction = 'tsne')
FeaturePlot(stc.combined, features = c('SOX2','PTPRZ1'),reduction = 'tsne')
stc.combined = subset(stc.combined, subset = seurat_clusters!=7) #7 may be malignant cells

DefaultAssay(stc.combined) <- "integrated"
stc.combined <- ScaleData(stc.combined, verbose = FALSE) %>%
  RunPCA( npcs = 30, verbose = FALSE) %>%
  RunUMAP( reduction = "pca", dims = 1:30) %>%
  RunTSNE( reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = 0.2)

stc.combined <- stc.combined %>% FindClusters(resolution = 0.3)

#DimPlot(stc.combined, label = T,reduction = 'tsne')
DimPlot(stc.combined, label = T,reduction = 'umap')
DefaultAssay(stc.combined) <- "SCT"
DotPlot(stc.combined,features = c('TOP2A','FCN1','IFIT1','CD163','P2RY12','CX3CR1','CCL4'))+coord_flip()
FeaturePlot(stc.combined,features = c('TOP2A','FCN1','IFIT1','CD163','P2RY12','CX3CR1','CCL4'))

mks <-stc.combined %>% PrepSCTFindMarkers() %>% FindAllMarkers()
mks$pctdiff = mks$pct.1-mks$pct.2
mks$logP = sign(mks$avg_log2FC)* -log10(mks$p_val)
mks %>%group_by(cluster) %>% slice_max(n = 3, order_by = pctdiff+logP) ->mks.top2
print(mks.top2, n = nrow(mks.top2))

VlnPlot(stc.combined, features = unique(mks.top2$gene)[1:5],stack = T, flip = T)
save(stc.combined, file = 'scg.scRNAseq.myeloid.Rda')
write.table(mks.top2, file = 'scg.scRNAseq.myeloid.top2markers.txt',row.names = F, quote = F, sep = "\t")
write.table(mks, file = 'scg.scRNAseq.myeloid.allmarkers.txt',row.names = F, quote = F, sep = "\t")


DefaultAssay(stc.combined) <-'SCT'
DimPlot(stc.combined, label = T,reduction = 'tsne')
DimPlot(stc.combined, label = T,reduction = 'umap')
FeaturePlot(stc.combined, features = c('CD3D','SOX2','PTPRZ1'),reduction = 'tsne')
gns = c('SPP1',#0, SPP1,
        'CD163',#1, CD163
        'PTCH2',#2, PTCH2
        'JUN',#3, AP-1 TF, mesenchymal
        'IL1A','CCL4',#4, pro-inflammatory
        'AIF1',#5, seems not really exist
        'HSPA1B',#6, hot-shock proteins, stress
        'LGALS1', #7, LGALS1
        'IFITM1', #8, IFITM1
        'MKI67',#9, MKI67, proliferative
        'EREG','FCER1A','FCN1', #10,TAMo: EREG,S100A12
        'LGALS3',#LGALS1
        'SERPINE1', #12
        'SOX2','PTPRZ1'
)
DotPlot(stc.combined, features = gns)+coord_flip()
fts_zemin = c('FCN1','S100A8','S100A9',#mono_CD14
              'FCGR3A','LST1','LILRB2',#mono_CD16
              'INHBA','CCL4','IL1RN',#macro_INHBA
              'NLRP3','EREG','IL1B', #macro_NLRP3
              'LYVE1','PLTP','SPP1', #macro_LYVE1
              'C1QC','C1QA','APOE' #macro_C1QC
)
DotPlot(stc.combined, features = fts_zemin)+coord_flip()

fts_xiaoqun = c('FCN1','S100A8','S100A9','CD52',# monocyte
              'C1QC','CX3CR1', #microglia
              'IFIT3','IFIT1','ISG15', #interferon activated
              'CCL2','CCL3','CCL4','NEDD9','IL1B', #pro-inflammatory
              'CD163','CD44','VEGFA',#mesenchymal,pro-angio
              'MKI67','TOP2A','TYMS' #proliferative
)
DotPlot(stc.combined, features = fts_xiaoqun)+coord_flip()

fts_diaz = read.delim('myeloid.66lineagemarkers.txt')
fts_diaz = c(fts_diaz$MG_Markers, fts_diaz$Mac_Markers)
fts_diaz = fts_diaz[fts_diaz!='']
DotPlot(stc.combined, features = fts_diaz)+coord_flip()

fts_yychen =c('LYZ','FCN1','VCAN',#blood-monocyte
              'EREG','AREG','CXCL2','MARCO','LRG1','VEGFA', #TAMO
              'BHLHE41','P2RY12','BIN1','NAV3','TAL1','SALL1','SIGLEC8','ADRB2',#microglia
              'HLA-DQA1','DSE','METRNL','FPR3','ARG2','IL10','CREM' #macrophage
              )

DotPlot(stc.combined, features = fts_yychen)+coord_flip()

DefaultAssay(stc.combined) <- "integrated"
stc.combined<- BuildClusterTree(stc.combined)
PlotClusterTree(stc.combined)

##compare by tumor group
propdf = as.data.frame(prop.table(table(stc.combined$orig.ident[stc.combined$seurat_clusters!=5],
                                        stc.combined$seurat_clusters[stc.combined$seurat_clusters!=5] ), margin = 1))
propdf$group = ifelse(propdf$Var1 %in% c('SP04','SP15','SP41','SP24'),'A','B')
library(ggpubr)
ggplot(propdf, aes(x = Var2, y = 100*Freq, fill = group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_dodge(width = 0.75))+
  scale_fill_manual(values = c('#C2A5F9','#F2A460'))+
  stat_compare_means(label = 'p.format')+
  theme_classic()+theme(axis.text.x = element_text(angle=30,hjust = 1))+
  labs(x = '', y = 'Proportion of myeloid cells (%)')

stc.combined$myloid.cell.state = stc.combined$cell.state
stc.combined$myloid.cell.state[stc.combined$seurat_clusters %in%c(0,1)]='myeloid_proinflammatory'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters %in%c(3,6)]='myeloid_stress'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==2]='myeloid_microglia-like'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==4]='myeloid_proangiogenic'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==5]='myeloid_neutrophil'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==7]='myeloid_monocyte'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==8]='myeloid_interferon-activated'
stc.combined$myloid.cell.state[stc.combined$seurat_clusters==9]='myeloid_proliferative'
DimPlot(stc.combined, group.by = 'myloid.cell.state', reduction = 'umap')
save(stc.combined, file = 'scg.scRNAseq.myeloid.withState.Rda')

load('scg.scRNAseq.myeloid.withState.Rda')
propdf = as.data.frame(prop.table(table(stc.combined$orig.ident,
                                        stc.combined$myloid.cell.state ), margin = 1))
propdf$group = ifelse(propdf$Var1 %in% c('SP04','SP15','SP41','SP24'),'A','B')
library(ggpubr)
propdf$Var2 = gsub('myeloid_','',propdf$Var2)
ggplot(propdf, aes(x = Var2, y = 100*Freq, fill = group))+
  geom_boxplot(outlier.shape = NA,show.legend = F)+
  geom_point(position = position_dodge(width = 0.75),cex=0.5,show.legend = T)+
  scale_fill_manual(values = c('#8da0cb','#fc8d62'))+
  stat_compare_means(label = 'p.signif',hide.ns = T,label.y.npc = 0.9)+
  scale_y_continuous(limits = c(0,60))+
  theme_classic()+theme(axis.text.x = element_text(angle=50,hjust = 1))+
  labs(x = '', y = 'Percentage (%)')

fts_myeloid_state = c('HSPA1A','HSPH1','HSPA6', #stress
                      'MKI67','TOP2A','TYMS', #proliferative,
                'IL1B','NEDD9','CCL3','CCL4','CD68', #pro-inflammatory
                'CD163','VEGFA',#mesenchymal,pro-angio #'CD44',
                'CSF3R','FCGR3B',#neutrophils
                'FCN1','CD52',# monocyte
                'C1QC','CX3CR1','P2RY12', #microglia
                'IFIT3','IFIT1','ISG15'#interferon activated
)
DotPlot(stc.combined, features =fts_myeloid_state,
        group.by = 'myloid.cell.state')+
  #coord_flip()+
  theme(axis.text.x = element_text(angle = 90,hjust=1))+
  labs(x = '',y='')
