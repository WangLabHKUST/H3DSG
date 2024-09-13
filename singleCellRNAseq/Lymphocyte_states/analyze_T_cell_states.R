library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)

setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/Tcells/')
load('~/Dropbox/communter/seurat_joint.Rda')
stc = subset(seurat_joint, subset = cell.state=='T')
rm(seurat_joint); gc()

table(stc$orig.ident) #note most samples have small number of cells!

stc.list <- SplitObject(stc, split.by = "orig.ident") #split by source

# normalize and identify variable features for each dataset independently
stc.list <- lapply(X = stc.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = stc.list)

#find integration anchors and perfrom integration
stc.anchors <- FindIntegrationAnchors(object.list = stc.list, anchor.features = features, k.filter = 10)
stc.combined <- IntegrateData(anchorset = stc.anchors,k.weight = 10)
DefaultAssay(stc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
stc.combined <- ScaleData(stc.combined, verbose = FALSE) %>% 
  RunPCA( npcs = 30, verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:20) %>% 
  #RunTSNE( reduction = "pca", dims = 1:10) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.4)
DefaultAssay(stc.combined) = 'RNA'
p1 <- DimPlot(stc.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(stc.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave(file = 'Integrated.UMAP.pdf', width = 8, height = 3.5)

stc.combined = subset(stc.combined, subset = seurat_clusters!=3)
stc.combined = subset(stc.combined, subset = seurat_clusters!=5)
DefaultAssay(stc.combined) <- "integrated"
stc.combined <- FindVariableFeatures(stc.combined) %>% 
  ScaleData( verbose = FALSE) %>% 
  RunPCA( npcs = 30, verbose = FALSE) %>% 
  RunUMAP( reduction = "pca", dims = 1:20) %>% 
  #RunTSNE( reduction = "pca", dims = 1:10) %>% 
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.4)
DefaultAssay(stc.combined) = 'RNA'
p1 <- DimPlot(stc.combined, reduction = "umap", group.by = "orig.ident", pt.size = 0.1)
p2 <- DimPlot(stc.combined, reduction = "umap", label = TRUE, pt.size = 0.1)+
  scale_color_brewer(type = 'qual', palette = 3)
p1 + p2
ggsave(file = 'Integrated.UMAP.pdf', width = 8, height = 3.5)


# identify markers for each subgroup
mks <-FindAllMarkers(stc.combined)
mks$pctdiff = mks$pct.1-mks$pct.2
mks$logP = sign(mks$avg_log2FC)* -log10(mks$p_val)
mks %>%group_by(cluster) %>% slice_max(n = 3, order_by = pctdiff+logP) ->mks.top2
print(mks.top2, n = nrow(mks.top2))
gns = mks.top2$gene[!duplicated(mks.top2$gene)]
gns = gns[!gns %in% c('TPT1','TRDC','LYPD3')]
VlnPlot(stc.combined, features = gns, stack = T, flip = T,group.by = 'seurat_clusters')+NoLegend()
ggsave(file = 'Violin.top3markers.pdf', width = 3.5, height =5.5)
FeaturePlot(stc.combined, features = c('CD8A',
                                       'CD4',
                                       'CD3D',
                                       'NKG7',
                                       'FOXP3',
                                       'CCR7',
                                       'CD40LG',
                                       'HSPA1B',
                                       'ISG15',
                                       'KLF2',
                                       'C1QA',
                                       'PTN',
                                       'STMN1',
                                       'GNLY',
                                       'KLRC1',
                                       #'KLR1' #CD161
                                       'FGFBP2'
                                       ),
            pt.size = 0.1)
ggsave(file = 'FeaturePlot.markers.pdf', width = 11.5, height = 7.5)

FeaturePlot(stc.combined, features = c('TIGIT',
                                       'PDCD1',
                                       'FAM3C',
                                       'CD80',
                                       'CD86',
                                       'CTLA4',
                                       #'TNFRSF14',
                                       'TNFRSF9', #=CD137
                                       'NCAM1','FCGR3A',
                                       'FGFBP2','CD8A','FOXP3'#,'KLRB1'
),pt.size = 0.1)
ggsave(file = 'FeaturePlot.immunoinhibitor.pdf', width = 11.5, height = 7.5)

stc.combined <- RenameIdents(stc.combined, `0` = "CD8 T", `1` = "CD4 memory T", `2` = "CD56dimCD16+ NK",
                                `3` = "stress T", `4` = "CD56+CD16- NK", `5` = "proliferative T", 
                             `6` = "regulatory T", `7` = "naive T", `8` = "interferon T")
DimPlot(stc.combined, label = T, repel = T,pt.size = 0.1, label.size = 3)+scale_color_brewer(type = 'qual', palette = 3)
ggsave(file = 'Celltype.labeled.pdf', width = 5.3, height = 3.1)

table(stc.combined$orig.ident)
save(stc.combined, file = 'scg.scRNAseq.Tcells.Rda')

load('scg.scRNAseq.Tcells.Rda')
propdf = as.data.frame(prop.table(table(stc.combined$orig.ident,Idents(stc.combined) ), margin = 1))
propdf$group = ifelse(propdf$Var1 %in% c('SP04','SP15','SP41','SP24'),'A','B')
library(ggpubr)
ggplot(propdf, aes(x = Var2, y = 100*Freq, fill = group))+
  geom_boxplot(outlier.size = 0.1)+scale_fill_manual(values = c('#8da0cb','#fc8d62'))+
  stat_compare_means(label = 'p.signif',hide.ns = T,label.y.npc = 0.98,col='red')+
  theme_classic()+theme(axis.text.x = element_text(angle=30,hjust = 1))+
  labs(x = '', y = 'Proportion of T cells (%)')
ggsave( file = 'Tcell.proportion.compareAB.pdf', width = 4.8, height = 3.6)

gns = c("CD8A","CD8B","GZMK","CD40LG","IL7R","FGFBP2","FCGR3A","KLRF1","HSPA1A","HSPA1B","CD160","KRT81","TMIGD2","NCAM1","ZWINT","TYMS","MKI67","FOXP3","IL2RA","TNFRSF4","CCR7","LEF1","IFIT2","IFIT1","ISG15")
VlnPlot(stc.combined, features = gns, stack = T, flip = T)+NoLegend()
DotPlot(stc.combined, features = gns)+theme(axis.text.x = element_text(angle = 60,hjust=1))+labs(x = '',y='')

ggsave(file = 'Violin.celltypesandtheirmarkers.pdf', width = 6.8, height =3.6)

#exhausiton markers
ems = c('PDCD1','TOX','CXCL13','TIGIT','CTLA4','TNFRSF9','HAVCR2','LAG3')
DotPlot(stc.combined, features = ems)+theme(axis.text.x = element_text(angle = 60,hjust=1))+labs(x = '',y='')
hmp <-DoHeatmap(subset(stc.combined, subset = orig.ident %in% c('SP05','SP11','SP13','SP20T5','SP21T5','SP22')),
                 features = ems,slot = 'data')
stc.combined$group = ifelse(stc.combined$sample.id %in% c('SP04','SP15','SP41','SP24'),'Pons-like','Thalamus-like')
FeaturePlot(stc.combined,features = c('CD8A','GZMK'),split.by = 'group', order = T, pt.size = 0.2)
DotPlot(stc.combined,features = c('CD8A','GZMK'),group.by = 'group')
stc.combined$tcellstate = Idents(stc.combined)
gedf = FetchData(stc.combined, vars = c('CD8A','GZMB','group','tcellstate','sample.id'),assay = 'RNA' ) 
gedf$group = factor(gedf$group, levels = c('Thalamus-like','Pons-like'))
ggplot(gedf, aes(x = CD8A, y = GZMB))+
  geom_point(alpha = 0.5,aes(color = group),show.legend = F)+
  scale_color_manual(values = c('#ba461e','#2a65a9'))+
  facet_wrap(~group,ncol=1)+
  theme_classic()
library(dplyr)
tmpdf1 = as.data.frame(prop.table(table( gedf$sample.id,gedf$tcellstate=='CD8 T'),margin=1))
tmpdf1$Group = ifelse(tmpdf1$Var1%in% c('SP04','SP15','SP41','SP24'),'Pons-like','Thalamus-like')
library(ggbeeswarm)
ggplot(tmpdf1[tmpdf1$Var2==TRUE,], aes(x = Group, y = 100*Freq))+
  geom_boxplot(aes(fill = Group),show.legend = F,outlier.shape = NA)+
  geom_quasirandom(show.legend = F)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons-like','Thalamus-like')),label.y = 65)+
  scale_fill_manual(values = c('#ba461e','#2a65a9'))+
  theme_classic()+labs(x = '',y='cytotxic CD8 T cells (%)')+lims(y=c(0,70))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))

tmpdf1 = as.data.frame(prop.table(table( gedf$sample.id,gedf$tcellstate=='CD8 T'&gedf$GZMB>0),margin=1))
tmpdf1$Group = ifelse(tmpdf1$Var1%in% c('SP04','SP15','SP41','SP24'),'Pons-like','Thalamus-like')
library(ggbeeswarm)
ggplot(tmpdf1[tmpdf1$Var2==TRUE,], aes(x = Group, y = 100*Freq))+
  geom_boxplot(aes(fill = Group),show.legend = F,outlier.shape = NA)+
  geom_quasirandom(show.legend = F)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons-like','Thalamus-like')),label.y = 65)+
  scale_fill_manual(values = c('#ba461e','#2a65a9'))+
  theme_classic()+labs(x = '',y='CD8 T cells (%)')+lims(y=c(0,70))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))



ggplot(tmpdf1[tmpdf1$Var2==TRUE,], aes(x = Group, y = 100*Freq))+
  geom_boxplot(aes(fill = Group),show.legend = F,outlier.shape = NA)+
  geom_quasirandom(show.legend = F)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons-like','Thalamus-like')),label.y = 65)+
  scale_fill_manual(values = c('#8da0cb','#fc8d62'))+
  theme_classic()+labs(x = '',y='CD8 T cells (%)')+lims(y=c(0,70))+
  theme(axis.text.x = element_text(angle = 30,hjust=1))



ggplot(hmp$data, aes(x= Cell, y = Feature))+
  geom_tile(aes(fill = Expression))+
  scale_fill_gradient2(high = 'red',low = 'blue',mid = 'blue',midpoint = 0)+
  #scale_fill_viridis_c()+
  theme_classic()+
  theme(axis.text.x = element_blank())
stc.combined <-AddModuleScore(stc.combined,features = list(ems))
FeaturePlot(stc.combined, features = 'Cluster1',order = T)

dt = hmp$data
ddt = reshape2::dcast(hmp$data,formula = Cell ~ Feature, value.var = 'Expression')
ddt$type = dt$Identity[match(ddt$Cell, dt$Cell)]

library(pheatmap)
rownames(ddt) = ddt$Cell
ddt = ddt[order(ddt$type, -ddt$PDCD1),]
pheatmap(t(ddt[,2:9]), annotation_col = subset(ddt, select = 'type'),
         cluster_rows = F, cluster_cols = F, show_colnames = F)

DotPlot(subset(seurat_joint,subset = cell.state !='unknown'), features = ems,group.by = 'cell.state')+theme(axis.text.x = element_text(angle = 60,hjust=1))+labs(x = '',y='')
