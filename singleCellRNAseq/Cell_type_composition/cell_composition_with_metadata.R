setwd('~/Documents/Research/SpinalcordGlioma/scRNA/rerun_preprocessing/')
fls = list.files(pattern = 'celltype.txt', path = '.',full.names = T,recursive = T)
for (i in 1:length(fls)){
  tmp = read.delim(fls[i])
  tb1 = as.data.frame(table(tmp$cell.type))
  tb2 = as.data.frame(table(tmp$cell.state[tmp$cell.type=='malignant']))
  
  tb1$Sample = tb2$Sample = strsplit(basename(fls[i]),"\\.")[[1]][1]
  
  if (i ==1){
    res1 = tb1; res2 = tb2
  }else{
    res1 = rbind(res1,tb1); res2 = rbind(res2,tb2)
  }
}
res1$Var1[res1$Var1=='plasma'] ='B'
res1 = res1[!res1$Var1 %in% c('mural','unknown','neuron','B'),]
res1$Prop = 0
for (ele in unique(res1$Sample)){
  sele = sum(res1$Freq[res1$Sample==ele])
  res1$Prop[res1$Sample==ele] = res1$Freq[res1$Sample==ele]/sele
}
library(ggplot2)
res1$Sample = gsub('T5','',res1$Sample)
res1$Sample[res1$Sample=='X172017'] = 'SP202'
res1$Sample[res1$Sample=='X177793'] = 'SP203'
res1$Sample[res1$Sample=='X747347'] = 'SP205'

res1$Var1 = factor(res1$Var1, levels = rev(c('malignant','myeloid','T','endothellial','neutrophil','oligodendrocyte')))
res1$Sample = factor(res1$Sample, levels = c("SP04" ,"SP15"  ,'SP24', "SP41" , "SP202","SP204",
                                             "SP05"  , "SP11" ,"SP13" ,  "SP20" ,"SP21" ,"SP22" , "SP201", "SP203" ,"SP205" ))
res1$smgroup = ifelse(res1$Sample %in% c("SP04" ,"SP15"  ,'SP24', "SP41" , "SP202","SP204"),'A','B')
p1<-ggplot(res1[res1$smgroup=='A',], aes(x = Sample, y = 100*Prop, fill = Var1))+
  geom_bar(position = "stack", stat = 'identity', width = 0.75, show.legend = F)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(type = 'qual', palette = 3)+
  labs(x = '', fill = '', y = 'Proportion (%)')
p2<- ggplot(res1[res1$smgroup=='B',], aes(x = Sample, y = 100*Prop, fill = Var1))+
  geom_bar(position = "stack", stat = 'identity', width = 0.75)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y = element_blank(), 
                        axis.line.y = element_blank(),
                        axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank())+
  scale_fill_brewer(type = 'qual', palette = 3)+
  labs(x = '', fill = '', y = 'Proportion (%)')
cowplot::plot_grid(p1,p2, rel_widths = c(0.35,0.65))
ggsave(file = 'celltype.composition.pdf', width = 6.36, height = 2.38)


res2$Prop = 0
for (ele in unique(res2$Sample)){
  sele = sum(res2$Freq[res2$Sample==ele])
  res2$Prop[res2$Sample==ele] = res2$Freq[res2$Sample==ele]/sele
}
res2$Sample = gsub('T5','',res2$Sample)
res2 = res2[!startsWith(res2$Sample,'X'),]
res2$Sample = factor(res2$Sample, levels = c("SP04" ,"SP15"  ,'SP24', "SP41" ,# "SP202","SP204",
                                             "SP05"  , "SP11" ,"SP13" ,  "SP20" ,"SP21" ,"SP22" #, "SP201", "SP203" ,"SP205"
                                             ))
res2$smgroup = ifelse(res2$Sample %in% c("SP04" ,"SP15"  ,'SP24', "SP41" , "SP202","SP204"),'A','B')
res2$Var1 = factor(res2$Var1, levels = c('AC-like','S','G2M','MES-like','OC-like','OPC-like'))
p1<-ggplot(res2[res2$smgroup=='A',], aes(x = Sample, y = 100*Prop, fill = Var1))+
  geom_bar(position = "stack", stat = 'identity', width = 0.75, show.legend = F)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c( '#377eb8','#e41a1c','#e41a1c','#984ea3','#4daf4a','#ff7f00'))+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x = '', fill = '', y = 'Proportion (%)')
p2<-ggplot(res2[res2$smgroup=='B',], aes(x = Sample, y = 100*Prop, fill = Var1))+
  geom_bar(position = "stack", stat = 'identity', width = 0.75)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c( '#377eb8','#fb8072','#e41a1c','#984ea3','#4daf4a','#ff7f00'))+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        axis.text.y = element_blank(), 
                        axis.line.y = element_blank(),
                        axis.ticks.y = element_blank(), 
                        axis.title.y = element_blank())+
  labs(x = '', fill = '', y = 'Proportion (%)')
cowplot::plot_grid(p1,p2, rel_widths = c(0.36,0.62))
ggsave(file = 'cellstate.composition.pdf', width =7.2, height = 2.38)

#
mtd = read.delim('../../Metadata/spinalcordglioma_singlecell_metadata.txt', row.names = 1)
p1<-ggplot(mtd[!mtd$orig.ident %in%c('SP20B5','SP21B5'),], aes(x = orig.ident, y = log10(nCount_RNA)))+geom_violin(fill = '#fed9a6')+theme_classic()+lims(y=c(0,5))+labs(y='nCount\n(log10)')
p2<-ggplot(mtd[!mtd$orig.ident %in%c('SP20B5','SP21B5'),], aes(x = orig.ident, y = log10(nFeature_RNA)))+geom_violin(fill = '#ccebc5')+theme_classic()+lims(y=c(0,5))+labs(y='nFeature\n(log10)')
p3<-ggplot(mtd[!mtd$orig.ident %in%c('SP20B5','SP21B5'),], aes(x = orig.ident, y = percent.mito))+geom_violin(fill = '#fbb4ae')+theme_classic()+lims(y=c(0,15))+labs(x = '')
cowplot::plot_grid(p1+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank()),
                   p2+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.title.x = element_blank()),
                   p3,ncol=1, rel_heights = c(0.3,0.3,0.4))
ggsave(file = 'DataQuality.perSample.pdf', width = 5.4, height = 3)
mtd = mtd[,!names(mtd) %in% c('Censorship','OS')]
mtd = mtd[order(rownames(mtd)=='SP20', mtd$MethylationGroup,mtd$H3.status=='H3.3 K27I',mtd$TP53, mtd$NF1, mtd$CDK4amp, mtd$CDK6amp, mtd$CDKN2Adel ),]

m1 = mtd[,1:20]
m1[is.na(m1)]=0
m1[m1== -2] = -4
m1[m1==-999]=NA

dan = mtd[,31:21]
dan[dan==1]='available'
dan[dan=='']=NA

library(pheatmap)
#
ph <- pheatmap(t(m1), cluster_rows = F, cluster_cols = F,legend = F,
                   annotation_col = dan, breaks = c(-4,-1,0,1,4),
                   gaps_col = c(6),fontsize_row = 9,
                   color = c('blue','#dddddd','darkgreen','red'),
                   border_color = 'black', na_col = '#888888', lwd = 0.5)



##
library(ggpubr)
res1$Subgroup = dan$MethylationGroup[match(res1$Sample, rownames(dan))]
res1$Subgroup[res1$Sample=='SP24'] = 'A'
res1$Var1 = factor(res1$Var1, levels = c('malignant','myeloid','endothellial','T','neutrophil','oligodendrocyte'))
ggplot(res1, aes(x = Subgroup, y = 100*Prop, fill = Subgroup))+
  facet_wrap(facets = ~Var1,scales = 'free_y', nrow = 1)+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75))+
  scale_fill_manual(values = c('#C2A5F9','#F2A460'))+
  stat_compare_means(label = 'p.format', comparisons = list(c('A','B')))+
  theme_classic()+theme( strip.background = element_blank())+
  labs(x = 'Subgroup', y = 'Proportion')
ggsave(file = 'Celltype.proportion.pdf', width = 9, height = 3.6)

res2$Subgroup = dan$MethylationGroup[match(res2$Sample, rownames(dan))]
res2$Subgroup[res2$Sample=='SP24'] = 'A'
res2$Var1 = factor(res2$Var1, levels = c('S','G2M','OPC-like','MES-like','AC-like','OC-like'))
ggplot(res2[!is.na(res2$Subgroup),], aes(x = Subgroup, y = 100*Prop, fill = Subgroup))+
  facet_wrap(facets = ~Var1,scales = 'free_y', nrow = 1)+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75), size =2)+
  scale_fill_manual(values = c('#C2A5F9','#F2A460'))+
  stat_compare_means(label = 'p.format', comparisons = list(c('A','B')))+
  theme_classic()+theme( strip.background = element_blank())+
  labs(x = 'Subgroup', y = 'Proportion')
ggsave(file = 'Cellstate.proportion.pdf', width = 9, height = 3.6)

##finally save the heatmap
pdf(file = 'spinalcordK27.singlecell.metadata2.pdf', width = 7,height = 6.55)
grid::grid.newpage()
grid::grid.draw(ph$gtable)
dev.off()
