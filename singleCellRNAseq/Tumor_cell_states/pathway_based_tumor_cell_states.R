library(yaGST)
library(pheatmap)
library(Seurat)
#
setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/CellTypeandStates///')
#1. extract tumor cell subpopulations and their average expression
    #"SP04","SP05","SP11","SP13","SP15","SP20T5","SP21T5","SP22","SP24","SP41"

load('../scRNA/rerun_preprocessing/SP04/seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat,label =T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('6','13','1','8','10','2','3','0','12','7')]
write.table(gesp04_mal, file = 'SP04.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP04', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP04.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")

load('../scRNA/rerun_preprocessing/SP05/seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('5','0','17','13','11','2','7','3','10','8')]
write.table(gesp04_mal, file = 'SP05.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP05', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP05.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")

load('../scRNA/rerun_preprocessing/SP11/seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('12','8','7','0','6','14','15')]
write.table(gesp04_mal, file = 'SP11.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP11', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP11.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP13/seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('14','13','15','8','11','4','18','3','2','1','6','5','9')]
write.table(gesp04_mal, file = 'SP13.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP13', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP13.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP15/seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('10','5','8','4','1','0','9','13','2','6')]
write.table(gesp04_mal, file = 'SP15.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP15', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP15.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP20T5//seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('10','7','14','3','2','4','13','6')]
write.table(gesp04_mal, file = 'SP20.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP20', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP20.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP21T5//seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('8','12','15','4','14')]
write.table(gesp04_mal, file = 'SP21.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP21', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP21.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP22//seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('7','6','11','18','0','4','8','17','9','12')]
write.table(gesp04_mal, file = 'SP22.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP22', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP22.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP24//seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('9','0','13','4','2','15','14','12','16','6','3','1','7','11','10','5')]
write.table(gesp04_mal, file = 'SP24.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP24', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP24.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


load('../scRNA/rerun_preprocessing/SP41//seurat.Rda')
DotPlot(seurat, features = c('PTPRC','SOX2','PTPRZ1','OLIG2','CD44','MKI67','MAG'),cluster.idents = T)
DimPlot(seurat, label = T)
gesp04 <- AverageExpression(seurat,features = rownames(seurat@assays$RNA),assays = 'RNA')$RNA
gesp04_mal = gesp04[,c('13','12','21','6','2','5','9','0','10','7','1','4','3','11')]
write.table(gesp04_mal, file = 'SP41.tumor.subpopulation.RNA.avgexp.txt',quote= F, sep = "\t")
gesp04_mal_cellcount = data.frame(Sample = 'SP41', Cluster = dimnames(gesp04_mal)[[2]])
gesp04_mal_cellcount$CellCount = table(seurat$seurat_clusters)[gesp04_mal_cellcount$Cluster]
gesp04_mal_cellcount$TotalCount = sum(gesp04_mal_cellcount$CellCount)
write.table(gesp04_mal_cellcount, file = 'SP41.tumor.subpopulation.RNA.cellcount.txt',quote= F,row.names = F, sep = "\t")


rm(seurat,gesp04,gesp04_mal,gesp04_mal_cellcount)

#2. assess pathway activities of each subpopulation
fls = list.files(pattern = 'tumor.subpopulation.RNA.avgexp.txt')
for (i in 1:length(fls)){
  tmp =read.delim(fls[i],check.names = F)
  names(tmp) = paste0(strsplit(fls[i],"\\.")[[1]][1],"_",names(tmp))
  tmp$gene = rownames(tmp)
  if(i==1){
    ge= tmp
  }else{
    ge = merge(ge, tmp,by = 'gene')
  }
}

rownames(ge) = ge$gene
ge = ge[,names(ge)!='gene']
gelog2 = log2(ge+1) 

library(yaGST)
hm = gmt2GO('h.all.v2023.2.Hs.symbols.gmt')
go= gmt2GO('c5.go.v2023.2.Hs.symbols.gmt')
kegg= gmt2GO('c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt')

hmgokegg = c(hm,kegg)

ssMwwGST <- function(geData, geneSet, nCores = 8){
  library(yaGST)
  library(doMC)
  
  means <- rowMeans(geData)
  sds <- apply(geData, 1, sd)
  
  registerDoMC(nCores)
  ans <- foreach(ss = 1:ncol(geData)) %dopar% {
    currentSample <- (geData[, ss] - means)/sds
    rankedList <- sort(currentSample, decreasing = T)
    
    aMwwGST <- lapply(geneSet, function(x) mwwGST(rankedList, geneSet = x, minLenGeneSet = 15, alternative = "two.sided", verbose = F))
    aMwwGST <- aMwwGST[sapply(aMwwGST, length) != 0]
    tmp_NES <- sapply(aMwwGST, function(x) x$log.pu)
    tmp_pValue <- sapply(aMwwGST, function(x) x$p.value)
    
    ans <- list(tmp_NES = tmp_NES, tmp_pValue = tmp_pValue)
    print(ss)
    return(ans)
  }
  NES <- sapply(ans, function(x) x$tmp_NES)
  pValue <- sapply(ans, function(x) x$tmp_pValue)
  colnames(NES) <- colnames(pValue) <- colnames(geData)
  
  FDR <- t(apply(pValue, 1, function(x) p.adjust(x, method = "fdr")))
  res <- list(NES = NES, pValue = pValue, FDR = FDR)
  return(res)
}

ans <- ssMwwGST(geData = gelog2,geneSet = c(hm,kegg)) #go

NES = ans$NES
FDR = ans$FDR

NES = NES[apply(NES,1,max)>0.3 & apply(FDR,1,min)<0.0001,]
#NES[NES>3]=3
library(pheatmap)
ph = pheatmap(NES,scale = 'row',clustering_method = 'ward.D2',
         fontsize_row = 4,fontsize_col = 5,
         col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75))

#now put cells into different states
fls = list.files(pattern = 'tumor.subpopulation.RNA.cellcount.txt')
for (i in 1:length(fls)){
  tmp =read.delim(fls[i],check.names = F)
  if(i==1){
    tcc= tmp
  }else{
    tcc = rbind(tcc, tmp)
  }
}
tcc$id = paste(tcc$Sample, tcc$Cluster,sep = '_')

tcc = merge(tcc,as.data.frame(cutree(ph$tree_col,k = 3)),by.x='id',by.y=0)
names(tcc)[6]='Group'
tcc$StateName = ifelse(tcc$Group==1,'MTC',ifelse(tcc$Group==2,'GPM','PPR'))
rownames(tcc)=tcc$id
pheatmap(NES,scale = 'row',clustering_method = 'ward.D2',annotation_col = subset(tcc,select = 'StateName'),#tcc[,c(2,7)],
              fontsize_row = 4,legend = T,ontsize_col = 5,show_colnames = F,show_rownames = F,border_color = NA,
         treeheight_row = 15,treeheight_col = 15,
              col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75),
         annotation_colors = list(StateName=c('GPM'='#962f23','MTC'='#f56942','PPR'='#f6c6ae')))
 
tcc$Prop = tcc$CellCount/tcc$TotalCount
tcc$SampleGroup = ifelse(tcc$Sample %in% c('SP04','SP15','SP24','SP41'),'Pons-like','Thalamus-like')

library(dplyr)
library(ggplot2)
tcc %>% group_by(SampleGroup, StateName) %>% 
  summarise(S = sum(CellCount)) %>% 
  ggplot(aes(x = SampleGroup, y = S,fill = StateName))+
  geom_bar(stat = 'identity',position = 'fill') +
  theme_classic()+scale_y_continuous(expand = c(0,0))+
  labs(x ='',y='Proportion' ,fill = 'Cell state')



library(ggpubr)
tcc_sample_state = expand.grid(Sample = unique(tcc$Sample),State = unique(tcc$StateName),Prop=0)
for (i in 1:nrow(tcc_sample_state)){
  tcc_sample_state$Prop[i] = sum(tcc$Prop[tcc$Sample==tcc_sample_state$Sample[i] & tcc$StateName==tcc_sample_state$State[i]])
}
tcc_sample_state$SampleGroup = ifelse(tcc_sample_state$Sample %in% c('SP04','SP15','SP24','SP41'),'Pons-like','Thalamus-like')
tcc_sample_state%>% 
  ggplot( aes(x = State, y = 100*Prop, fill = SampleGroup))+
  geom_boxplot()+geom_point(position = position_dodge(width=0.75))+
  stat_compare_means(label = 'p.format')+
  theme_classic()+
  labs(x = 'State', y = 'Proportion (%)',fill='')

ct =read.delim('celltypeandstates_bygroup.txt')
ct$StateName=factor(ct$StateName, levels = c('GPM','MTC','PPR','B','T/NK','myeloid','neutrophil','Oligo','vascular'))
ggplot(ct, aes(x = SampleGroup, y = S, fill =StateName))+
  geom_bar(stat = 'identity',position = 'fill',width=0.75 )+
  scale_fill_manual(values = c('#c53932','#ef8a62','#fdae61','#519e3e','#3c75af','#66c2a5','#a6dba0','#ffffb3','#bebada'))+
  scale_y_continuous(expand = c(0,0),breaks = 0.25*0:4,labels = 25*0:4)+
  labs(x = '',y='Cell proportion (%)',fill='')+theme_classic()+
  theme(axis.text = element_text(color = 'black'), axis.text.x = element_text(angle = 30,hjust=1))

load('../scRNA/rerun_preprocessing/SP41/seurat.Rda')
seurat$clusterid = paste('SP41',seurat$seurat_clusters,sep = '_')
seurat$newcellstate = ifelse(seurat$clusterid %in% tcc$id,tcc$StateName[match(seurat$clusterid ,tcc$id)],seurat$cell.state)
oldsts = c('S','G2M','AC-like','OC-like','MES-like','OPC-like')
seurat = subset(seurat, cells = rownames(seurat@meta.data)[!seurat$newcellstate %in% oldsts])
DimPlot(seurat,group.by = 'newcellstate', label = T)
#DotPlot(subset(seurat,subset = newcellstate %in%c('MTC','GPM','PPR')), features = c('SOX2','OLIG2','CD44','MKI67'),group.by = 'newcellstate')
write.table(seurat@meta.data[,c('cell.barcode','newcellstate')], file = 'SP41.cell.antoniostate.txt',row.names = F, quote = F, sep = "\t")
saveRDS(seurat, file = 'SP41.addnewstatelabel.Rda')
#now calculate spatial transcriptomics
sp2 = read.delim('ST/sp22.exp.txt')
ans_sp2 <- ssMwwGST(geData = sp2,geneSet = c(hm,kegg)) #go
NES = ans_sp2$NES
FDR = ans_sp2$FDR

NES = NES[apply(NES,1,max)>0.3 & apply(FDR,1,min)<0.0001,]
#NES[NES>3]=3
library(pheatmap)
ph = pheatmap(NES,clustering_method = 'ward.D',scale = 'row',
              fontsize_row = 4,fontsize_col = 5,show_colnames = F,
              col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75))

sp24 = read.delim('ST/')
ans_sp2 <- ssMwwGST(geData = sp2,geneSet = c(hm,kegg)) #go
NES = ans_sp2$NES
FDR = ans_sp2$FDR

NES = NES[apply(NES,1,max)>0.3 & apply(FDR,1,min)<0.0001,]
#NES[NES>3]=3
library(pheatmap)
ph = pheatmap(NES,clustering_method = 'ward.D',scale = 'row',
              fontsize_row = 4,fontsize_col = 5,show_colnames = F,
              col = colorRampPalette(c("#053061", "#4393C3", "cornsilk", "#D6604D", "#67001F"))(75))



#try a different method
mks = read.delim('Antonio_4state_markers.txt')
smps = c('SP04','SP05','SP11','SP13','SP15','SP20T5','SP21T5','SP22','SP24','SP41')
for(i in 1:length(smps)){
  load(paste0('../scRNA/rerun_preprocessing/',smps[i],'/seurat.Rda'))
  seurat = AddModuleScore(seurat, features = as.list(mks),name = 'AntonioState',search = T)
  ai_states = names(mks)
  seurat$AntonioState1_label = ai_states[apply(seurat@meta.data[,paste0('AntonioState',1:4)],1,which.max)]
  mals = c('S','G2M','AC-like','OC-like','MES-like','OPC-like')
  if(i==1){
    tb = table(seurat$AntonioState1_label[seurat$cell.state%in%mals])
  } else{
    tb = rbind(tb,table(seurat$AntonioState1_label[seurat$cell.state%in%mals]))
  }
  
}
rownames(tb)=smps
tb = as.data.frame(tb)
tb$Group = ifelse(rownames(tb) %in% c('SP04','SP15','SP24','SP41'),'Pons-like','Thalamus-like')


