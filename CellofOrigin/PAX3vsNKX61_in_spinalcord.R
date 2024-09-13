setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/')

#directly start from here
library(Seurat)
library(ggplot2)
g1 = c('SP04','SP15','SP24','SP41','X172017','SP204')
g2 = c('SP05','SP11','SP13','SP20T5','SP21T5','SP22','X177793','X747347')
load('~/Dropbox/communter/seurat_joint.Rda')
so = subset(seurat_joint, subset = cell.state %in% c('AC-like','OC-like','MES-like','OPC-like','S','G2M'))
so$cell.state=factor(so$cell.state, levels = c('S','G2M','OPC-like','AC-like','MES-like','OC-like'))
so1 = subset(so, subset = sample.id %in% g1)
so2 = subset(so, subset = sample.id %in% g2)

DotPlot(so1, features = c('NKX6-1','PAX3'), group.by = 'cell.state')+
     coord_flip()+theme(axis.text.x = element_text(angle=45,hjust = 1))

DotPlot(so2, features = c('NKX6-1','PAX3'), group.by = 'cell.state')+
  coord_flip()+theme(axis.text.x = element_text(angle=45,hjust = 1))
