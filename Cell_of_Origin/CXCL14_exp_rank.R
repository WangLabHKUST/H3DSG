#

setwd('~/Dropbox/Mac (2)//Documents/Research/SpinalcordGlioma/HOX/')
df0 = read.delim('Raw_SP_AB_RNA_count_data_with_gene_length.csv', sep = ",")
gty = read.delim('~/Dropbox/2020-PanGliomaEvolution/Collaboration_Shengshuo/Homo_sapiens.GRCh37.75.genes.typeNlengthNCoord.txt.gz')
library(edgeR)
gl = df0$Genelength
df = as.data.frame(rpkm(df0[,3:42],gene.length = gl))
rownames(df) = df0$Geneid
#write.table(df,file = 'SPAB_bulkRNA.RPKM.tsv',sep = "\t")

df1 = df[grepl('^HOX', rownames(df)),]
df1 = df1[!grepl('-', rownames(df1)),]


library(pheatmap)
pheatmap(log2(df1+1))

df2 =df[rownames(df) %in%c('NKX6-1','PAX3'),]

df2 = as.data.frame(t(df2))
boxplot(df2$PAX3, df2$`NKX6-1`)

df2$Patient = rownames(df2)
library(ggplot2)
ggplot(df2, aes(x = PAX3, y =`NKX6-1`))+
  geom_point()+lims(x = c(0,60), y =c(0,60))+
  geom_abline(slope = 1, intercept = 0, lty=2)+
  theme_classic()

cl = read.delim('H3K27altered_project_metadata_updated.tsv')
cl$ID = paste0('P', cl$ID2)
df2$Location = cl$Location[match(df2$Patient, cl$ID)]
df2$Cluster = cl$TSNE.cluster[match(df2$Patient, cl$ID)]

ggplot(df2, aes(x = PAX3, y =`NKX6-1`, color = Cluster))+
  geom_point()+lims(x = c(0,60), y =c(0,60))+
  geom_abline(slope = 1, intercept = 0, lty=2)+
  theme_classic()

ann = data.frame(ID = names(df))
ann = merge(ann, cl, by = 'ID')
rownames(ann) = ann$ID
ann$Location = substr(ann$Location,1,1)
pheatmap(log2(df1+1), annotation_col = ann[,3:6],cluster_rows = F,
         fontsize_row = 6,fontsize_col = 6,border_color = 'white')
tmp = data.frame(gene = rownames(df1)); tmp$Family = substr(tmp$gene,1,4);
tmp$int = as.integer(gsub('HOX[ABCD]','',tmp$gene));tmp = tmp[order(tmp$Family, tmp$int),]

pheatmap(log2(df1[tmp$gene,]+1), annotation_col = ann[,3:6],clustering_method = 'complete',
         cluster_rows = F,gaps_row = c(11,21,30),show_colnames = F,
         fontsize_row = 6,fontsize_col = 6,border_color = '#333333')


dft =as.data.frame(t(df))
dft$group = cl$TSNE.cluster[match(rownames(dft), cl$ID)]
deg = data.frame(gene = names(dft)[1:(ncol(dft)-1)], logFC = 0, pval = 1)
for (i in 1:nrow(deg)){
  x = log2(dft[dft$group=='A',i]+1)
  y = log2(dft[dft$group=='B',i]+1)
  
  deg$logFC[i]=median(y)-median(x)
  deg$pval[i]=wilcox.test(x,y)$p.value
}

deg$geneType = gty$geneType[match(deg$gene, gty$geneName)]

degpc = deg[deg$gene %in% gty$geneName[gty$geneType=='protein_coding'],]
degpc$lab = ifelse(degpc$pval<0.01 & abs(degpc$logFC)>2,degpc$gene,NA)
degpc$type = ifelse(degpc$pval<0.05 & degpc$logFC>1,'up',ifelse(degpc$pval<0.05 & degpc$logFC< -1, 'dn','nosig'))
library(ggrepel)
ggplot(degpc, aes(x = logFC, y = -log10(pval), label = lab))+
  geom_point(aes(color = type, size= type), show.legend = F)+
  geom_text_repel(size=2,segment.size=0.25,min.segment.length = 0,max.overlaps = 10)+
  geom_vline(xintercept = 0,lty=2,col='#888888')+
  geom_hline(yintercept = -log10(0.05),lty=2,col='#888888')+
  scale_color_manual(values = c('blue','#cccccc','red'))+
  scale_size_manual(values = c(1,0.5,1))+
  theme_classic()+
  labs(x = 'log2 (Thalamus-like/Pons-like)', y = '-log10(P-value)')

degpc$rnk =rank(degpc$logFC,ties.method = 'random')
ggplot()+
  geom_point(data = degpc[degpc$logFC>0,], mapping = aes(x = rnk, y = logFC, color = type, size= -log10(pval) ))+
  geom_text_repel(data = degpc[degpc$logFC>0,], mapping = aes(x = rnk, y = logFC, label = lab, size= -log10(pval) ),
                  size=2.5,segment.size=0.25,min.segment.length = 0,max.overlaps = 15,xlim = c(22000,NA),direction = 'both',show.legend = F)+
  geom_point(data = degpc[degpc$logFC<0,], mapping = aes(x = rnk, y = logFC, color = type, size= -log10(pval) ))+
  geom_text_repel(data = degpc[degpc$logFC<0,], mapping = aes(x = rnk, y = logFC, label = lab,size= -log10(pval) ),
                  size=2.5,segment.size=0.25,min.segment.length = 0,max.overlaps = 15,xlim = c(1500,NA),direction = 'both',show.legend = F)+
  #geom_vline(xintercept = 0,lty=2,col='#888888')+
  geom_hline(yintercept = 0,lty=2,col='#888888')+
  scale_color_manual(values = c('blue','#cccccc','red'))+
  scale_size_continuous(range = c(0.01,2),guide='none')+
  #scale_x_log10()+
  #scale_color_viridis_c()+
  #scale_size_manual(values = c(1,0.5,1))+
  scale_x_continuous(limits = c(0,28000),breaks = c(0,10000,20000))+
  scale_y_continuous(limits = c(-4,4),breaks = c(-4,-2,2,4))+
  theme_classic()+
  labs(y = 'log2 (Thalamus-like/Pons-like)', x = 'rank by fold change')

library(ggbeeswarm)
library(ggpubr)
dft$group2 = ifelse(dft$group=='A','Pons-like','Thalamus-like')
ggplot(dft, aes(x = group2 , y = log2(CXCL14+1)))+
  geom_boxplot(aes(fill = group2),outlier.shape = NA,show.legend = F )+
  scale_fill_manual(values = c("#75a5c9","#e2906a"))+
  geom_quasirandom(width = 0.15, alpha=0.7)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons-like','Thalamus-like')))+
  theme_classic()+labs(x ='')

ggplot(dft, aes(x = group2 , y = log2(ARSF+1)))+
  geom_boxplot(aes(fill = group2),outlier.shape = NA,show.legend = F )+
  scale_fill_manual(values = c("#75a5c9","#e2906a"))+
  geom_quasirandom(width = 0.15, alpha=0.7)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons-like','Thalamus-like')))+
  theme_classic()+labs(x ='')


chemokine_genes = degpc$gene[grepl("^CCL|^CXCL|^CX3CL|^XCL",degpc$gene)]
chemokine_genes = chemokine_genes[chemokine_genes!='CCL15-CCL14']
chemokine_genes = c(chemokine_genes,c('PF4','PF4V1','PPBP','IL8'))
chemokine_genes = c('CCL1','CCL2','CCL3','CCL3L1',"CCL3L3","CCL4","CCL4L1", "CCL4L2" ,"CCL5",
                    "CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19",
                    "CCL20","CCL21" , "CCL22",  "CCL23" , "CCL24" , "CCL25",  "CCL26" , "CCL27" , "CCL28",
                    "CXCL1","CXCL2" , "CXCL3" ,"PF4","PF4V1", "CXCL5","CXCL6" ,'PPBP' ,"IL8" ,"CXCL9",
                    "CXCL10" ,"CXCL11", "CXCL12", "CXCL13", "CXCL14" ,"CXCL16" ,"CXCL17",
                    "CX3CL1","XCL1","XCL2")
deg_chemokine = degpc[degpc$gene %in% chemokine_genes,]
rownames(deg_chemokine) = deg_chemokine$gene
deg_chemokine = deg_chemokine[rev(chemokine_genes),]
deg_chemokine$padj = p.adjust(deg_chemokine$pval)
deg_chemokine = deg_chemokine[order(deg_chemokine$padj,decreasing = T),]
deg_chemokine$rank = 1:nrow(deg_chemokine)
deg_chemokine$lab = ifelse(deg_chemokine$type!='nosig',deg_chemokine$gene,NA)

ggplot(deg_chemokine, aes(x = 49-rank, y = -log10(padj) , label = lab ))+
  geom_line(color = '#888888')+geom_point(aes(size=logFC, fill = -log10(padj)),pch=21,show.legend = T)+
  #geom_text_repel(size=3,segment.size=0.25,min.segment.length = 0,max.overlaps =15,xlim = c(5,NA),nudge_y = 0.25,direction = 'x')+
  scale_fill_gradient2(low = '#cccccc',mid = '#eeeeee',high = 'red',midpoint = 1)+
  scale_x_continuous(limits = c(0.25,nrow(deg_chemokine)+0.5),breaks = 1:nrow(deg_chemokine), labels = rev(deg_chemokine$gene), expand = c(0,0))+
  labs(x = 'Chemokines', y = '-log10 (P-adj)')+lims(y=c(-1,4))+
  scale_size_continuous(range = c(0.1,3.5))+
  geom_hline(yintercept = -log10(0.05),lty=2,col='#dddddd')+
  #coord_flip()+
  theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position="right",
                         axis.text.x.bottom = element_blank())+
  theme(legend.position = "right") +
  guides(
    fill = 'none',
    size = guide_legend(ncol = 1)    # 将 size 图例显示为一列
  )

dft_chemokine = log2(dft[,chemokine_genes]+1)

dft_chemokine = dft_chemokine[,deg_chemokine$gene]
dft_chemokine = dft_chemokine[rownames(dft)[order(dft$group,dft$CXCL14,decreasing = T)],]
pheatmap(dft_chemokine[,chemokine_genes],cluster_cols = F, cluster_rows = F,show_rownames   = F,
         color = colorRampPalette(rev(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6")))(50),
         annotation_row  = subset(dft, select = 'group'), border_color = 'white',gaps_row = c(24),
         annotation_colors = list(group=c(A='#75a5c9',B='#e2906a')),fontsize_col = 7.5)


#focus on HOX
library(reshape2)
df1$gene = rownames(df1)
#df1$gene = factor(df1$gene, levels = c(paste0('HOXA',1:13),paste0('HOXB',1:13),paste0('HOXC',1:13),paste0('HOXD',1:13)))
df1 = df1[order(df1$gene),]
df1mt = melt(df1)
df1mt$family = substr(df1mt$gene,1,4)
df1mt$gene = factor(df1mt$gene, levels = c(paste0('HOXA',13:1),paste0('HOXB',13:1),paste0('HOXC',13:1),paste0('HOXD',13:1)))
df1mt$Location = cl$Location[match(df1mt$variable, cl$ID)]
df1mt$Cluster = cl$TSNE.cluster[match(df1mt$variable, cl$ID)]
df1mt$family2 = as.integer(gsub('HOX[ABCD]','',df1mt$gene))
write.table(df1mt, file = 'spinalDMG.HOX.expression.txt',row.names = F, quote = F, sep = "\t")
ggplot(df1mt, aes(x = gene, y = log2(value+1)))+
  geom_boxplot(aes(color = family2, fill = family2), alpha=0.5,
               show.legend = F, width=0.5,outlier.shape = NA,outlier.size = 0.5)+
  geom_jitter(aes(color = family2),pch=20,cex=0.1,width=0.1, show.legend = F)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lwd=0.25,lty=2)+
  facet_grid(Cluster~.)+
  #scale_x_discrete(position = "top") +
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

ggplot(df1mt, aes(x = gene, y = log2(value+1)))+
  geom_boxplot(aes(color = family))+
  facet_grid(Cluster~.)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),panel.grid = element_blank())+
  labs(x = '', y = 'expression (log2 RPKM)')
##
spcmks = c("PAX6","NKX2-2","NKX2-8","FOXA2","LMX1A","IRX3","OLIG2","NEUROG3","PTF1A","DBX2","NKX6-1","GDF7","NEUROG2","DBX1","NKX6-2","PAX7","LHX3","PAX3","MSX1","NEUROG1","OLIG3","GSX1","NEUROG2","ATOH1","NEUROG1","ASCL1","GATA2","BMPR1A","BMPR1B","GSX2","FOXN4","POU4F1","LBX1","EVX1","EVX2","EN1","DLL4","NOTCH1","SOX1","MNX1","SIM1","BARHL1","FOXD3","TLX3","PAX2","TLX1","VSX2","TAL1","ISL1","ISL2","LHX9","LHX1","LHX5","LMX1B","SOX14","GATA2","LHX2","DRGX","LHX5","POU4F1","SOX21","GATA3","PTF1A","GSX1","GSX2","BHLHE22","PITX2","DLL1","PLCG1","WT1","FOXP2")
spcmks = unique(spcmks)
resdeg = data.frame(gene = rownames(df),grpA = 0, grpB = 0, pval = 1, log2FC = 0)
for (i in 1:nrow(df)){
  x = as.numeric(df[i,names(df) %in% cl$ID[cl$TSNE.cluster=='A']])
  y = as.numeric(df[i,names(df) %in% cl$ID[cl$TSNE.cluster=='B']])
  resdeg$grpA[i] = median(x)
  resdeg$grpB[i] = median(y)
  resdeg$pval[i] = wilcox.test(x,y)$p.value
}
resdeg$log2FC = log2(resdeg$grpB+1) - log2(resdeg$grpA+1)
resdeg = resdeg[resdeg$gene %in% gty$geneName[gty$geneType=='protein_coding'],]
humtfs = read.delim('http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt', header = F)

resdeg$isTF = ifelse(resdeg$gene %in% humtfs$V1, 'TF','non-TF')
resdeg$isSPMK = ifelse(resdeg$gene %in% spcmks, 'TF','non-TF')
resdeg$lab = ifelse(resdeg$isTF=='TF' &resdeg$pval<0.05 & abs(resdeg$log2FC)>1, resdeg$gene,NA)
resdeg$lab2 = ifelse(resdeg$isSPMK=='TF' &resdeg$pval<0.05 , resdeg$gene,NA)
resdeg$lab3 = ifelse(resdeg$pval<0.05 & abs(resdeg$log2FC)>2, resdeg$gene,NA)


library(ggrepel)
ggplot(resdeg, aes(x = log2FC ,y=-log10(pval), color = isTF, label = lab))+
  geom_point(cex = 0.1)+geom_text_repel(size=3,color = 'black', segment.size = 0.25)+
  theme_classic() + scale_color_manual(values = c('#dddddd','red'))

ggplot()+
  geom_point(resdeg[resdeg$isSPMK=='non-TF',], mapping=aes(x = log2FC ,y=-log10(pval)), color = '#dddddd')+
  geom_point(resdeg[resdeg$isSPMK=='TF',], mapping=aes(x = log2FC ,y=-log10(pval), label = lab2), color = 'red')+
  geom_text_repel(resdeg[resdeg$isSPMK=='TF',], mapping=aes(x = log2FC ,y=-log10(pval), label = lab2),size=3,color = 'black', segment.size = 0.25)+
  theme_classic() + scale_color_manual(values = c('#dddddd','red'))

ggplot(resdeg, aes(x = log2FC ,y=-log10(pval), label = lab3))+
  geom_point(cex = 0.1)+geom_text_repel(size=3,color = 'black', segment.size = 0.25)+
  theme_classic() + scale_color_manual(values = c('#dddddd','red'))


tmp = as.data.frame(t(df[c('CXCL14','PRAME','DBX2'),]))
tmp$group = cl$TSNE.cluster[match(rownames(tmp), cl$ID)]
library(ggbeeswarm)
ggplot(tmp, aes(x = group, y = log2(CXCL14)))+
  geom_boxplot(aes(color = group))+
  geom_quasirandom()+
  scale_color_brewer(palette=6,type = 'qual')+
  theme_classic()

ggplot(tmp, aes(x = group, y = log2(DBX2)))+
  geom_boxplot(aes(color = group))+
  geom_quasirandom()+
  scale_color_brewer(palette=6,type = 'qual')+
  theme_classic()
###
dfspmk = df[rownames(df) %in% spcmks,]
pheatmap(log2(dfspmk[rownames(dfspmk)!='SOX14',]+1), annotation_col = ann[,3:7],border_color = 'white',
         fontsize_row = 4,fontsize_col = 4,clustering_method = 'ward.D',
         cluster_rows = T,cluster_cols = T, scale = 'row')

df3 =df[rownames(df) %in%c('NKX6-1','PAX3',
                           'DBX2',#'SOX4',
                           'IRX1','CBX2',
                           'HES5','SOX11'),]
pheatmap(log2(df3+1), annotation_col = ann[,3:7])
gns = c('LMX1A','MSX1','MSX2','PAX3','WNT1',#RP
        'OLIG3','IRX3','IRX5','PAX6', #dp1-2
        'PAX7','GSX2','GBX2','GSX1',#dp3-4 'ASCL1',
        'DBX2','DBX1','SP8', #dp5-6
        'NKX6-2','PRDM12',#p0-p1
        'NKX6-1','FOXN4', #p2
        'NKX2-2',#'NKX2-9', #p3, pMN 'OLIG2',
        'FOXA2','FERD3L','ARX','SHH','LMX1B' #FP
)

df4 =df[gns,]
pheatmap(log2(df4+1), annotation_col = ann[,3:7], cluster_rows = F)

df5 = df[c(gns, rownames(df1)),]
pheatmap(log2(df5+1), annotation_col = ann[,3:7], cluster_rows = F, clustering_method = 'ward.D')


#immune genes
immunegns = c("CTLA4","PDCD1","IDO1","TGFB1","TGFBR1","BTLA","CD160","CD274","LAG3","PDCD1LG2",
              "TNFRSF14","IL10","IL10RB","CD80","ICOSLG","ICOS","TNFSF9","TNFRSF9","TNFRSF4",
              "TNFSF4","CD70","TNFSF18","IL6","TMEM173","CD27","TNFRSF13B","TNFRSF17","IL6R",
              "CD86","CD28","TNFRSF18","TNFSF13B","TNFRSF13C","TNFSF13","CD40LG","CD40",
              "HLA-A","TAP1","TAP2","HLA-B","HLA-C","B2M","HLA-DRB6","HLA-DQB2","HLA-DPB2",
              "HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DPA1","HLA-DPB1",
              "HLA-G","HLA-E","HLA-F",'CD47','SIRPA')
df6 = df[immunegns,]
df6 = df6[!is.na(df6[,1]),]
ann.row = data.frame(Type = c(rep('Immunoinhibitor',13),
                              rep('Immunosimulator',23),
                              rep('MHC class I',6), 
                              rep('MHC class II',10),
                              rep('MHC non-class',3), 
                              rep('Dont eat me!',2)), 
                     row.names = immunegns)
library(pheatmap)
pheatmap(log2(df6+1), annotation_col = ann[,3:7], annotation_row = ann.row,cluster_rows = F,
         scale = 'row', clustering_method = 'ward.D', show_colnames = F)

immunegns_compare = data.frame(gene = rownames(df6))
for (ix in 1:nrow(df6)){
  tmp = as.data.frame(t(df6[ix,]))
  names(tmp) = 'genei'
  tmp$ID = rownames(tmp)
  tmp = merge(tmp, ann)
  #p<-ggplot(tmp, aes(x = TSNE.cluster, y = genei))+geom_boxplot()+geom_quasirandom()+stat_compare_means()+theme_classic()
  #print(p+labs(y=rownames(df6)[ix]))
  immunegns_compare$groupA[ix] = median(tmp$genei[tmp$TSNE.cluster=='A'])
  immunegns_compare$groupB[ix] = median(tmp$genei[tmp$TSNE.cluster=='B'])
  immunegns_compare$pval[ix] = wilcox.test(genei~TSNE.cluster, data = tmp)$p.value
  
}
immunegns_compare$Type = ann.row$Type[match(immunegns_compare$gene, rownames(ann.row))]
cig = read.delim('cancer_immune_genes.txt')
immunegns_compare$alias = cig$alias[match(immunegns_compare$gene, cig$gene)]
immunegns_compare$alias[is.na(immunegns_compare$alias)]=immunegns_compare$gene[is.na(immunegns_compare$alias)]
immunegns_compare$lab = ifelse(immunegns_compare$pval<0.01, immunegns_compare$alias,NA)
ggplot(immunegns_compare, aes(x = log2(groupB/groupA), y = -log10(pval), label = lab))+
  geom_point(aes(color = Type),size=2)+
  geom_vline(xintercept = 0,lty=2)+
  geom_text_repel(size=2, segment.size=0.25,min.segment.length = 0)+
  lims(x = c(-3.8,3.8))+
  theme_classic()+labs(color = '')
ggsave(file = 'immune_cancer_genes.comparison.pdf', width = 4.8, height = 3.2)
write.table(immunegns_compare, file = 'immune_cancer_genes.comparison.txt', row.names = F, quote = F, sep = "\t")
