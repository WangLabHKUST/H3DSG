library(scRepertoire)
setwd('~/Documents/Projects/H3DSG/Revision/CXCL14_scTCR/reanalyze/scTCRseq//')

#load data
r0 <- read.csv("Tcell.filtered_contig_annotations.csv")
r1 <- read.csv("TCXCL14.filtered_contig_annotations.csv")
r2 <- read.csv("TNET.filtered_contig_annotations.csv")

contig_list <- list(r0,r1,r2)

#combine TCR contigs into clones
combined.TCR <- combineTCR(contig_list, samples = c("Tcell",'TCXCL14','TNET'),
                           removeNA = FALSE, removeMulti = FALSE,  filterMulti = FALSE)
head(combined.TCR[[1]])

exportClones(combined.TCR, write.file = TRUE, file.name = "TCR.clones.csv")
             
#basic clone analysis
clonalQuant(combined.TCR,  cloneCall="aa", chain = "both",  scale = TRUE)
clonalCompare(combined.TCR,  cloneCall = "aa",top.clones = 20)+guides(fill = 'none')

library(ggplot2)
tb = clonalCompare(combined.TCR,cloneCall = 'aa',top.clones = 1000,exportTable = T)#,  top.clones = 10, cloneCall="aa",  graph = "alluvial")
tb = tb[order(tb$Sample, tb$Proportion, decreasing = T),]

tb$idx[tb$Sample=='Tcell'] = 1:sum(tb$Sample=='Tcell')
tb$idx[tb$Sample=='TCXCL14'] = 1:sum(tb$Sample=='TCXCL14')
tb$idx[tb$Sample=='TNET'] = 1:sum(tb$Sample=='TNET')
tb$Sample = factor(tb$Sample, levels = c('Tcell','TNET','TCXCL14'))
tb1 = tb[tb$clones %in% tb$clones[tb$Sample=='TCXCL14' & tb$idx<=20 ],]
ggplot(tb, aes(x = idx, y = Proportion))+
  geom_point(aes(color = Sample))+
  #geom_line(aes(group = Sample))+
  geom_line(aes(group = clones))+
  theme_classic()
library(ggsci)
library(ggbeeswarm)
ggplot(tb1, aes(x = Sample, y = 100*Proportion))+
  geom_quasirandom(aes(color = Sample),width = 0.2, alpha=0.75,show.legend = F,size=2)+
  #geom_line(aes(group = Sample))+
  geom_line(aes(group = clones),color = '#cccccc',linewidth=0.5,alpha = 0.8)+
  scale_color_npg()+
  stat_compare_means(comparisons = list(c('Tcell','TCXCL14'),c('TNET','TCXCL14')))+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'))+
  labs(y='clonal proportion (%)')

library(ggalluvial)
tb1$clones = factor(tb1$clones, levels = tb1$clones[tb1$Sample=='TCXCL14'])
library(viridis)
ggplot(tb1, aes(x = Sample, fill = clones, group = clones, 
                   stratum = clones, alluvium = clones, y = 100*Proportion, 
                   label = clones)) + 
  scale_y_continuous(expand = c(0,0))+theme_classic() + 
  theme(axis.title.x = element_blank(), axis.text = element_text(color = 'black'),
        legend.text = element_text(size = rel(0.5)), legend.key.size = unit(0.5, "line")) +
  geom_stratum(size=0.1,color = 'black') +  geom_flow(stat = "alluvium")+
  scale_fill_viridis_d(option="inferno",direction = -1)+
  #guides(fill = 'none')+
  labs(y='Proportion of CD8+ T cells (%)')

clonalHomeostasis(combined.TCR,  cloneCall = "aa")
clonalProportion(combined.TCR, cloneCall = "aa") 

#sequence analysis 
percentAA(combined.TCR, chain = "TRB", aa.length = 20)
positionalEntropy(combined.TCR,  chain = "TRA", aa.length = 20)

vizGenes(combined.TCR, x.axis = "TRBV",y.axis = NULL, plot = "barplot",  scale = TRUE)

percentGenes(combined.TCR,  chain = "TRA", gene = "Vgene")

percentKmer(combined.TCR, cloneCall = "aa",chain = "TRB", motif.length = 3,  top.motifs = 25)

#diversity analysis
clonalDiversity(combined.TCR, cloneCall = "aa")
clonalRarefaction(combined.TCR, plot.type = 1,hill.numbers = 0, n.boots = 2)
clonalOverlap(combined.TCR,  cloneCall = "aa",  method = "morisita")

#combine with scRNAseq data
library(SingleCellExperiment)
library(Seurat)
load('../scRNAseq/seurat_joint.harmony.Rda')
sce <- Seurat::as.SingleCellExperiment(JoinLayers(seurat_joint_harmony))
sce$cell.barcode = gsub('TCell','Tcell',sce$cell.barcode)
rownames(colData(sce)) <- sce$cell.barcode
sce.combined <- combineExpression(combined.TCR, sce, cloneCall="gene",  proportion = TRUE)
sce.combined$orig.ident = factor(sce.combined$orig.ident, levels = c('TCell','TNET','TCXCL14'))
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
library(scater)
plotUMAP(sce.combined, colour_by = "clonalFrequency",size_by = "clonalFrequency") #+
  #scale_color_manual(values=rev(colorblind_vector[c(4,3,2,2)]))
as.data.frame(colData(sce.combined)) ->cdt


dr = as.data.frame(reducedDim(sce.combined,type = 'UMAP'))
cdt = merge(cdt,dr,by=0)
cdt = cdt[order(cdt$clonalFrequency),]
cdt$Freq2 = ifelse(cdt$orig.ident=='TCell',cdt$clonalFrequency/8849,
                   ifelse(cdt$orig.ident=='TCXCL14',cdt$clonalFrequency/8571,cdt$clonalFrequency/11094))
cdt$orig.ident = factor(cdt$orig.ident, levels = c('TCell','TNET','TCXCL14'))

ggplot()+
  geom_point(cdt[is.na(cdt$Freq2),],mapping = aes(x = umap_1, y = umap_2),alpha = 0.9,size=0.1,color = '#cccccc')+
  geom_point(cdt[!is.na(cdt$Freq2),],mapping = aes(x = umap_1, y = umap_2,color  = 100*Freq2,size= 100*Freq2),alpha = .5)+
  scale_color_gradientn(colours = rev(c("#040404","#3E134F","#851170","#C53270","#F36E35","#F8B83C","#FFFE9E")),na.value='#cccccc')+
  scale_size_continuous(range = c(0.1,2),limits = c(0.025,0.125),breaks = 0.25*c(1,3,5))+
  facet_wrap(~orig.ident)+ 
  theme_classic()+#guides(color = 'none')+
  theme(strip.background = element_blank())+
  labs(color = 'clonal proportion (%)')


#
library(scCustomize)
seurat_joint_harmony$cellState = ifelse(seurat_joint_harmony$seurat_clusters %in% c(2,3,4,7),'proliferative',
                                        ifelse(seurat_joint_harmony$seurat_clusters %in% c(0,6),'memory-like',
                                               ifelse(seurat_joint_harmony$seurat_clusters==5,'effector-like','activated')))
seurat_joint_harmony$orig.ident = factor(seurat_joint_harmony$orig.ident, levels = c('TCell','TNET','TCXCL14'))
seurat_joint_harmony$cellState = factor(seurat_joint_harmony$cellState,levels = c('effector-like','memory-like','proliferative','activated'))
DimPlot_scCustom(seurat_object = seurat_joint_harmony, group.by = 'cellState',
                 colors_use = pal_npg()(4),figure_plot = TRUE,label = F,alpha = 0.5)

DimPlot_scCustom(seurat_object = seurat_joint_harmony, group.by = 'cellState',split.by = 'orig.ident',
                 colors_use = pal_npg()(4),figure_plot = TRUE,label = F,alpha = 0.5)

gns = c('CD3D','CD8A','PDCD1','CTLA4','FOXP3','GNLY','CXCL10','HIST1H3C','MAF','KLRB1',
        'MCM2','TOP2A','UBE2C','GZMB','GZMA','TNFAIP3','CCR7','IL32','TCF1','KLRG1','IL7R','SELL')
gns = c('CD3D','CD8A','KLRB1','EGLN3','SELL','TCF7','CCR7','IL7R','MKI67','TOP2A','IL2RA','CXCR6')

DotPlot_scCustom(seurat_joint_harmony, gns, group.by = 'cellState',scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #coord_flip()+
  labs(x = '', y = '', fill = 'avg.exp', size = 'pct.exp')
gns = c('CD3D','CD8A','KLRB1','EGLN3','CTLA4','PDCD1','SELL','TCF7','CCR7','IL7R','MKI67','TOP2A','IL2RA','CXCR6')

DotPlot_scCustom(seurat_joint_harmony, gns, group.by = 'cellState',scale = T) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  #coord_flip()+
  labs(x = '', y = '', fill = 'avg.exp', size = 'pct.exp')



dtt = as.data.frame(table(seurat_joint_harmony$cellState,seurat_joint_harmony$orig.ident))
library(ggalluvial)
ggplot(dtt, aes(x = Var2, y = Freq, fill = Var1))+
  geom_bar(stat = 'identity',position = 'fill',width = 0.5)+
  scale_fill_brewer(type = 'qual',palette = 6)+ theme_bw()
dtt2 = as.data.frame(prop.table(table(seurat_joint_harmony$cellState,seurat_joint_harmony$orig.ident),margin = 2))
dtt2$Var1 = factor(dtt2$Var1, levels = c('effector-like','memory-like','proliferative','activated'))
ggplot(dtt2,aes(x = Var2, y = Freq, stratum = Var2, alluvium = Var1)) +
  geom_bar(stat = 'identity',position = 'fill',width = 0.35,aes(fill = Var1),color = '#222222',size=0.2)+
  geom_alluvium(aes(fill = Var1), width = 1/3,alpha=0.5) +
  scale_y_continuous(expand = c(0,0),breaks = c(0,1,2,3,4)/4,labels = 25*c(0,1,2,3,4))+
  scale_fill_npg()+theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = '',y='Proportion of CD8 T cells (%)', fill = '')


cdt$cellState = ifelse(cdt$seurat_clusters %in% c(2,3,4,7),'proliferative',
                       ifelse(cdt$seurat_clusters %in% c(0,6),'memory-like',
                              ifelse(cdt$seurat_clusters==5,'effector-like','activated')))
library(ggpubr)
library(ggbeeswarm)
cdt$orig.ident = factor(cdt$orig.ident, levels = c('TCell','TNET','TCXCL14'))
ggplot(cdt[!is.na(cdt$clonalFrequency),], aes(y = orig.ident, x = 100*clonalProportion ))+
  geom_bar(stat = 'identity',position = 'stack', aes(fill = Freq2),show.legend = T,width = 0.5)+
  #geom_boxplot(outlier.shape = NA,aes(color = cellState))+
  #geom_jitter(position = position_dodge(width = 0.75),aes(color = cellState))+
  #geom_quasirandom(aes(color = cellState),width = 0.1,dodge.width = 0.1,alpha = 0.2,orientation = 'x')+
  facet_wrap(cellState~.,scales = 'free_x')+
  scale_x_continuous(expand =c(0,0))+
  scale_fill_viridis_c(option="inferno",direction = -1)+
  #stat_compare_means(comparisons = list(c('TCell','TCXCL14'),c('TNET','TCXCL14')))+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),strip.background = element_blank(),axis.title.y = element_blank())+
  labs(x = 'Proportion of CD8 T cells (%)',fill = 'clone size')


ggplot(cdt[!is.na(cdt$clonalFrequency),], aes(x = orig.ident, y = clonalProportion ))+
  #geom_boxplot(outlier.shape = NA,aes(color = cellState))+
  geom_quasirandom(aes(color = cellState),width = 0.1,dodge.width = 0.1,alpha = 0.2,orientation = 'x')+
  #facet_grid(~cellState)+
  stat_compare_means(comparisons = list(c('TCell','TCXCL14'),c('TNET','TCXCL14')))+
  theme_classic()


