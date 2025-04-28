library(Seurat)
library(ggplot2)
library(dplyr)
#1-load our knowledge
kb = read.delim('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/HOX/HumanSpinalCordCellTypeKnowledge.txt')
kb = kb[kb$Neural_pop%in% c(paste0('dp',1:6),paste0('p',c(0:3,'MN')),paste0('dl',1:6),'V0','V1','V2a','V2b','V3','MN','FP','RP'),]
for (i in 1:nrow(kb)){
  tmp0 = strsplit(kb$Genes_map_step1[i],", ")[[1]]
  if(kb$Genes_map_step2[i]!=''){
    tmp = list(strsplit(kb$Genes_map_step2[i],", ")[[1]])
    names(tmp) = kb$Neural_pop[i]
  }
  if(i==1){
    mklist=tmp
    mklist0 = tmp0
  }else{
    mklist = c(mklist, tmp)
    mklist0 = c(mklist0, tmp0)
  }
}

mklist0 = unique(mklist0); mklist0 = mklist0[!is.na(mklist0)]
unique(unlist(mklist))
gns = unique(c(mklist0,unique(unlist(mklist))))

#2, load cell annotations
mt =read.delim('~/Dropbox/communter/HumanSpinalGW4toGW7_GSE171890.MetaData.csv',sep = ",")
mt$cellid = paste0(mt$orig.ident,'_',substr(mt$X,1,16),"-1")
mt$Type_step2[is.na(mt$Type_step2)]=mt$Type_step1[is.na(mt$Type_step2)]
mt$Type_step2 = gsub(' ','_', mt$Type_step2)

# recluster the cells
so = readRDS('~/Documents/Projects/H3DSG/Revision//HumanSpinalGW4toGW7_GSE171890.integrated.RDS')
so$Type_step1 = mt$Type_step1[match(rownames(so@meta.data),mt$cellid)]
so$Type_step2 = mt$Type_step2[match(rownames(so@meta.data),mt$cellid)]
Idents(so) <-so$Type_step2
#mksall = FindAllMarkers(object = so)

refso1 = subset(so, features = gns)
refso1 = subset(refso1 ,subset = Type_step2 %in% c(paste0('dp',1:6),paste0('dl',1:6),paste0('p',c(0:3,'MN')),'V0','V1','V2a','V2b','V3','MN','FP','RP'))
#DefaultAssay(refso1) = 'integrated'
refso1 <- refso1 %>% 
  NormalizeData(assay = 'RNA',normalization.method = 'CLR') %>% 
  # regress out variables which are sources of unwanted variation, and z-score data
  ScaleData(features = gns,assay = 'RNA')

refso1 <- refso1 %>%
  # compute PCA, based on scaled data
  RunPCA(assay = 'RNA',features = gns,
         npcs            = 50) 
  #RunTSNE(check_duplicates = FALSE,dims = 1:30) %>% 

refso1 <- refso1 %>% RunUMAP(dims = 1:30, verbose = FALSE, seed.use = 42)

refso1$Type_step2 = factor(refso1$Type_step2, levels = rev(c('RP',paste0('dp',1:6),paste0('p',c(0,1,2,'MN',3)),'FP',paste0('dl',1:6),'V0','V1','V2a','V2b','MN','V3')))
DimPlot(refso1, label = T, group.by = 'Type_step2',reduction = 'umap',
        cols=c('#66c2a5',"#edf8e9","#c7e9c0","#a1d99b","#74c476","#31a354","#006d2c",
               "#feedde","#fdbe85","#fd8d3c","#d94701","#de2d26",'#ffd92f',
               "#f1eef6","#d0d1e6","#a6bddb","#74a9cf","#2b8cbe","#045a8d",
               "#f2f0f7","#cbc9e2","#9e9ac8","#756bb1","#54278f",'#df65b0'))
library(ggplot2)
DotPlot(refso1, features = gns,
        cluster.idents = F, group.by = 'Type_step2')+
  theme(axis.text.x = element_text(angle=90,hjust=1))


DotPlot(refso1, features = c('SOX2','PAX3','NKX6-1','ELAVL3'),
        cluster.idents = F, group.by = 'Type_step2')+coord_flip()+
  theme(axis.text.x = element_text(angle=90,hjust=1))+labs(x='',y='')

table(refso1$Type_step2)
refso1$seurat_clusters_bckp = refso1$seurat_clusters
Idents(refso1) <-refso1$Type_step2
#mks = FindAllMarkers(refso1)

DotPlot(refso1, features = unique(mks$gene),
        cluster.idents = F, group.by = 'Type_step2')+
  theme(axis.text.x = element_text(angle=90,hjust=1))


refge1 = AverageExpression(so, group.by = 'Type_step1')
refge2 = AverageExpression(so, group.by = 'Type_step2',assays = 'integrated',slot = 'data')

pax3nkx6.1 = FetchData(so, vars = c('PAX3','NKX6-1','PLP1'),slot = 'data')
pax3nkx6.1$celltype1 = mt$Type_step1[match(rownames(pax3nkx6.1),mt$cellid)]
pax3nkx6.1$celltype2 = mt$Type_step2[match(rownames(pax3nkx6.1),mt$cellid)]
pax3nkx6.1$celltype1[which(pax3nkx6.1$PLP1>1& pax3nkx6.1$celltype1=='Progenitor')]='Oligodendrocyte'

pax3nkx6.1_neuron = pax3nkx6.1[which(pax3nkx6.1$celltype1=='Neuron'),]
pax3nkx6.1_neuron = pax3nkx6.1_neuron[order(pax3nkx6.1_neuron$PAX3, -pax3nkx6.1_neuron$`NKX6-1`,decreasing = T),]
library(pheatmap)
pheatmap(t(pax3nkx6.1_neuron[pax3nkx6.1_neuron$PAX3>0 | pax3nkx6.1_neuron$`NKX6-1`>0,1:2]),
         show_colnames = F, cluster_rows = F, cluster_cols = F, legend = F)


pax3nkx6.1_neuron = pax3nkx6.1[which(pax3nkx6.1$celltype1=='Progenitor'),]
pax3nkx6.1_neuron = pax3nkx6.1_neuron[order(pax3nkx6.1_neuron$PAX3, -pax3nkx6.1_neuron$`NKX6-1`,decreasing = T),]
pheatmap(t(pax3nkx6.1_neuron[pax3nkx6.1_neuron$PAX3>0 | pax3nkx6.1_neuron$`NKX6-1`>0,1:2]),
         show_colnames = F, cluster_rows = F, cluster_cols = F, legend = F)



pax3nkx6.1_neuron = pax3nkx6.1[which(pax3nkx6.1$celltype1=='Peripheral neuron'),]
pax3nkx6.1_neuron = pax3nkx6.1_neuron[order(pax3nkx6.1_neuron$PAX3, -pax3nkx6.1_neuron$`NKX6-1`,decreasing = T),]
pheatmap(t(pax3nkx6.1_neuron[pax3nkx6.1_neuron$PAX3>0 | pax3nkx6.1_neuron$`NKX6-1`>0,1:2]),
         show_colnames = F, cluster_rows = F, cluster_cols = F, legend = F)

pax3nkx6.1_neuron = pax3nkx6.1[which(pax3nkx6.1$celltype1=='Oligodendrocyte'),]
pax3nkx6.1_neuron = pax3nkx6.1_neuron[order(pax3nkx6.1_neuron$PAX3, -pax3nkx6.1_neuron$`NKX6-1`,decreasing = T),]
pheatmap(t(pax3nkx6.1_neuron[pax3nkx6.1_neuron$PAX3>0 | pax3nkx6.1_neuron$`NKX6-1`>0,1:2]),
         show_colnames = F, cluster_rows = F, cluster_cols = F, legend = F)

#
load('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/SP04/seurat.Rda')
sp4 = seurat
sp4ae = AverageExpression(sp4, group.by = 'cell.state',slot = 'data')

load('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/SP05/seurat.Rda')
sp5 = seurat
sp5ae = AverageExpression(sp5, group.by = 'cell.state')


df1 = as.data.frame(refge2$integrated)
df1$gene = rownames(df1)

df2 = as.data.frame(sp4ae$SCT)
df2$gene = rownames(df2)
tmp = merge(df1, df2, by = 'gene')
library(pheatmap)
pheatmap(cor(tmp[,2:52], tmp[,53:ncol(tmp)],method = 'pearson'), display_numbers = T)

df3 = as.data.frame(sp5ae$RNA)
df3$gene = rownames(df3)
tmp2 = merge(df1, df3, by = 'gene')
pheatmap(cor(tmp2[,2:52], tmp2[,53:ncol(tmp2)],method = 'spearman'),display_numbers = T)


#
refso1$Type_step1 = mt$Type_step1[match(rownames(refso1@meta.data),mt$cellid)]
refso1$Type_step2 = mt$Type_step2[match(rownames(refso1@meta.data),mt$cellid)]
DimPlot(refso1, label = T, group.by = 'Type_step1')

DimPlot(so, label = T, group.by = 'Type_step2')


DimPlot(refso1, label = T, group.by = 'Type_step2')

##
so1 = subset(so, subset = Type_step2 %in% c(paste0('dp',1:6),paste0('p',c(0:3,'MN')),'FP','RP'))

fds = list.dirs(path = '~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/')
fds = fds[3:length(fds)]
fds = fds[!grepl('B',fds)]
for (fd in fds){
  print(fd)
  load(paste0(fd,'/seurat.Rda'))
  DefaultAssay(so1) <- "integrated"
  DefaultAssay(seurat) <- "SCT"
  proj.anchors <- FindTransferAnchors(reference = so1, query = seurat,normalization.method = 'SCT',
                                      dims = 1:30, reference.reduction = "pca")
  predictions <- TransferData(anchorset = proj.anchors, refdata = so1$Type_step2,
                              dims = 1:30)
  seurat <- AddMetaData(seurat, metadata = predictions)
  tmp = as.data.frame(table(seurat$predicted.id, seurat$cell.state))
  tmp$SID = basename(fd)
  if(fd==fds[1]){
    resdf = tmp
  }else{
    resdf = rbind(resdf, tmp)
  }
}

g1 = c('SP04','SP15','SP24','SP41')
g2 = c('SP05','SP11','SP13','SP20T5','SP21T5','SP22')

resdf$g2 = ifelse(resdf$Var2 %in%c('AC-like','OPC-like','MES-like','OC-like','S','G2M'),'tumor cells','normal' )
resdf = resdf[resdf$SID %in% c(g1,g2),]
resdf$tumorGroup=ifelse(resdf$SID %in% g1,'GroupA','GroupB')
resdf$Var1 = factor(resdf$Var1, levels = c('RP','dp1','dp2','dp3','dp4','dp5','dp6','p0','p1','p2','pMN','p3','FP'))
resdf[resdf$g2=='tumor cells',]%>% 
  group_by( tumorGroup, Var2) %>% 
  mutate(Perc = Freq/sum(Freq)) %>% 
  group_by( tumorGroup, Var1) %>% 
  summarise(Frequency = sum(Perc)) ->tmp
library(reshape2)
tmp = rbind(tmp, data.frame(tumorGroup = c('GroupA','GroupA','GroupB'),Var1 = c('dp6','dp5','dp5'), Frequency = c(0,0,0)))
tmp$Var1 = factor(tmp$Var1, levels = rev(c('RP','dp1','dp2','dp3','dp4','dp5','dp6','p0','p1','p2','pMN','p3','FP')))
tmp %>% 
  ggplot(aes(y = Var1, fill=Frequency/6, x = tumorGroup ))+
  #geom_bar(position = 'fill', stat='identity',width = 0.75)+
  geom_tile(color = 'white',linewidth = 0.5)+
  scale_fill_viridis_c()+
  #scale_fill_gradientn(colours = rev(c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0")))+
  theme_classic()+theme(axis.line = element_blank(), axis.ticks = element_blank(),axis.text.x = element_text(angle=45,hjust=1))+
  scale_x_discrete(labels =c('Pons-like','Thalamus-like'))+
  #scale_y_continuous(expand = c(0,0))+
  #+
  theme(axis.text.y = element_text(color = rev(c("#8C1A11","#C6402D","#DE6E51","#ED9364","#ED9364","#F3BE8C","#FCF1DC",
                                             "#F3C7C2","#EDA3B5","#E671A0","#CC4494","#9F207B","#6F1473")),face="bold"))+
  labs(x = '', fill = 'Probability', y = '')

# resdf[resdf$g2=='tumor cells'&resdf$tumorGroup=='GroupB',]%>% 
#   group_by( SID, Var2) %>% 
#   mutate(Perc = Freq/sum(Freq)) %>% 
#   ggplot(aes(x = SID, y=Perc, fill = Var1 ))+
#   geom_bar(stat='identity')+facet_grid(~Var2)
# 
# #
# ggplot(resdf, aes(x = SID, y = Var1))+
#   geom_point(aes(size = Freq))+
#   theme_bw()
