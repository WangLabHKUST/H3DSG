setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/')
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)



fls = list.files('scK27Mmutation',pattern = 'K27M.reads.txt')
for (i in 1:length(fls)){
  tmp = read.delim(paste0('scK27Mmutation/',fls[i]), sep = ",")
  tmp$SID = strsplit(fls[i],"\\.")[[1]][1]
  if (i ==1){
    k27 = tmp
  }else{
    k27 = rbind(k27, tmp)
  }
}

tmp = read.delim('scK27Mmutation/SP22.H3F3B.K27I.reads.txt', sep = ",")
tmp$SID = 'SP22'
k27 = rbind(k27, tmp)
k27$af = ifelse(k27$depth_Total>=3,k27$depth_Alt/k27$depth_Total,NA)
k27$cell_ID = paste0(k27$SID,"_",k27$cell_ID, "-1")

#seurat@meta.data$depth_Alt = log2(k27$depth_Alt[match(rownames(seurat@meta.data), k27$cell_ID)]+1)
#seurat@meta.data$freq_Alt = k27$af[match(rownames(seurat@meta.data), k27$cell_ID)]

#seurat1 = subset(seurat,cells = rownames(seurat@meta.data)[!is.na(seurat$freq_Alt)])
#p0 <- DimPlot(seurat1, label = T)
#p1<- FeaturePlot(object = seurat1, order = T,features = c('depth_Alt','freq_Alt','PTPRZ1','MOG','PLP1','MSR1'))+scale_color_gradientn(colours = c('white','blue'), na.value = "#dddddd")
#cowplot::plot_grid(p0,p1)

library(scCustomize)


#
s = readRDS('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/CellTypeandStates/Spinalcord.10scRNAseq.integrated.withnewstated.RDS')
#DimPlot(s, group.by = 'cellType')
s$depth_Tot = k27$depth_Total[match(s$cell.barcode, k27$cell_ID)]
s$depth_Alt = k27$depth_Alt[match(s$cell.barcode, k27$cell_ID)]
s$freq_Alt = k27$af[match(s$cell.barcode, k27$cell_ID)]

s$logDALT = log2(s$depth_Alt + 1)
FeaturePlot_scCustom( s, features = c('freq_Alt'),order = T,pt.size = 0.1,na_cutoff = 0.4)

m = s@meta.data
m$depth_Alt[is.na(m$depth_Alt)]=0
m$freq_Alt[is.na(m$freq_Alt)]=0
table(m$orig.ident[m$celltype=='malignant'], m$depth_Alt[m$celltype=='malignant']>1)

m$detected = ifelse( m$depth_Alt>=1  ,'yes','no')
sm = as.data.frame(prop.table(table(m$orig.ident[m$celltype=='malignant'], m$detected[m$celltype=='malignant']),margin=1))
sm$cellType = 'malignant'
sm2 = as.data.frame(prop.table(table(m$orig.ident[m$celltype!='malignant'], m$detected[m$celltype!='malignant']),margin=1))
sm2$cellType = 'non-malignant'

sm = rbind(sm, sm2)
sm$tech = ifelse(sm$Var1 %in% c('SP20T5','SP21T5'),'T5','T3')
library(ggsci)
library(ggpubr)
ggplot(sm[sm$Var2=='yes',], aes(x = cellType, y = Freq))+
  geom_boxplot(outlier.shape = NA,aes(color =cellType ),show.legend = F, width = 0.5)+
  stat_compare_means(comparisons = list(c('malignant','non-malignant')),label.y = 0.85)+
  geom_jitter(width = 0.15,show.legend =F,alpha = 0.5)+
  scale_color_npg()+
  theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x='',y = 'proportion of cells with\n H3 K27 alteration',shape = 'Protocol')
