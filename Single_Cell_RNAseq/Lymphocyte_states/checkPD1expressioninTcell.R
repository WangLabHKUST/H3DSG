setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/Tcells/')

library(Seurat)
load('scg.scRNAseq.Tcells.Rda')

stc.combined$Group = ifelse(stc.combined$orig.ident %in% c('SP04','SP15','SP24','SP41'),'Pons-like','Thalamus-like')
DimPlot(stc.combined,split.by = 'Group')
FeaturePlot(stc.combined,features = c('CD8A','PDCD1','TGFB1','CTLA4'), order = T,pt.size = 0.2,ncol=4)
FeaturePlot(stc.combined,features = c('CD8A','PDCD1','TGFB1','CTLA4'), order = T,pt.size = 0.2,ncol=4)

FeaturePlot(stc.combined,features = c('CD8A','PDCD1','LAG3'),
            order = T,pt.size = 0.2,ncol=2, split.by = 'Group')

FeaturePlot(stc.combined,features = c('TGFB1','CTLA4'),
            order = T,pt.size = 0.2,ncol=4, split.by = 'Group')


FeaturePlot(stc.combined,features = c('CD8A','TIGIT'),
            order = T,pt.size = 0.2,ncol=4, split.by = 'Group')
# load('~/Dropbox/communter/seurat_integrated.Rda')
# #FeaturePlot(seurat_integrated,features = c('CXCL14','CD274','HAVCR2'), order = T)
# #DotPlot(seurat_integrated,features = c('CXCL14','CD274','HAVCR2','MPO'),group.by = 'cell.state')
# DimPlot(seurat_integrated, group.by = 'cell.state')
