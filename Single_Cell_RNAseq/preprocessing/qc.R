library(Seurat)
setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/CellTypeandStates/')
so = readRDS('Spinalcord.10scRNAseq.integrated.withnewstated.RDS')
#VlnPlot(so,features = 'nCount_RNA', group.by = 'sample.id')+
mt = so@meta.data
library(ggplot2)
ggplot(mt, aes(x = sample.id, y = nFeature_RNA))+
  geom_violin(aes(fill = sample.id), show.legend = F)+scale_y_log10()+theme_classic()
ggplot(mt, aes(x = sample.id, y = nCount_RNA))+
  geom_violin(aes(fill = sample.id), show.legend = F)+scale_y_log10()+theme_classic()
ggplot(mt, aes(x = sample.id, y = percent.mito))+
  geom_violin(aes(fill = sample.id), show.legend = F)+scale_y_log10()+theme_classic()
DimPlot(so)
gns =c('PTPRZ1','SOX2','PTPRC','AIF1','CD14','FCGR3B','S100A9',
       'CSF3R','CD3D','IL32','NKG7','CD79A','MS4A1',
       'VWF','PLP1','MAG')
FeaturePlot(so, features = gns, order = T)
f = read.delim('~/Documents/Projects/H3DSG/mal.prop.txt')
f$mal_prop = f$mal_a/f$mal_n
f$nmal_prop = f$nmal_a/f$nonmal_n
library(reshape2)
library(ggsci)
f2 = melt(f[,c('sample','protocol','mal_prop','nmal_prop')], id.vars = c('sample','protocol'))
ggplot(f2, aes(x = variable, y = 100*value ))+
  geom_boxplot(aes(color = variable),show.legend = F)+
  geom_jitter(aes(shape = protocol),size=2,alpha=0.7)+
  scale_color_aaas()+
  theme_classic()+
  labs(x = '', y ='%cells carrying H3 K27M' )
