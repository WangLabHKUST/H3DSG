setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/HOX/')
library(ggplot2)
hox = read.delim('GrCh3775_HOX.genes.txt', header = F)
flis = list.files('bulkATAC',pattern = 'bedGraph')
fljs = list.files('bulkH3K27Me3/',pattern = 'bedGraph')
for (i in 1:8){
  fli = flis[i]
  fac = read.delim(paste0('bulkATAC/',fli), header = F)
  #limit to the 500bp upstream of the promoter
  # for (gni in unique(fac$V5)){
  #   tssi = ifelse(hox$V5[hox$V7==gni]=='+',hox$V3[hox$V7==gni],hox$V4[hox$V7==gni] )
  #   tssj = ifelse(hox$V5[hox$V7==gni]=='+',tssi-500,tssi+500 )
  #   faci = fac[fac$V5==gni & (fac$V2-tssi)*(fac$V3-tssj)<0,]
  #   if (gni==unique(fac$V5)[1]){
  #     facnew = faci
  #   }else{
  #     facnew = rbind(facnew,faci)
  #   }
  # }
  # fac = facnew
  #</>limit to 500 bp
  dffac = aggregate(V4 ~ V5, FUN = mean, data = fac)
  dffac$HOX = substr(dffac$V5,1,4)
  dffac$int = as.integer(substr(dffac$V5,5,7))
  dffac = dffac[order(dffac$HOX, dffac$int),]
  dffac$V5 = factor(dffac$V5, levels = dffac$V5)
  dffac$type = 'ATAC'
  dffac$ID = strsplit(fli, "_")[[1]][1]
  
  flj = fljs[i]
  fme = read.delim(paste0('bulkH3K27Me3/',flj), header = F)
  # #limit to 500bp upstream of tss
  # for (gni in unique(fme$V5)){
  #   tssi = ifelse(hox$V5[hox$V7==gni]=='+',hox$V3[hox$V7==gni],hox$V4[hox$V7==gni] )
  #   tssj = ifelse(hox$V5[hox$V7==gni]=='+',tssi-500,tssi+500 )
  #   fmei = fme[fme$V5==gni & (fme$V2-tssi)*(fme$V3-tssj)<0,]
  #   if (gni==unique(fme$V5)[1]){
  #     fmenew = fmei
  #   }else{
  #     fmenew = rbind(fmenew,fmei)
  #   }
  # }
  # fme = fmenew
  # #</>limit to 500 bp
  dffme = aggregate(V4 ~ V5, FUN = mean, data = fme)
  dffme$HOX = substr(dffme$V5,1,4)
  dffme$int = as.integer(substr(dffme$V5,5,7))
  dffme = dffme[order(dffme$HOX, dffme$int),]
  dffme$V5 = factor(dffme$V5, levels = dffme$V5)
  dffme$type = 'H3K27Me3'
  dffme$ID = strsplit(flj, "_")[[1]][1]
  
  dfi = rbind(dffac, dffme)
  if (i ==1){
    df = dfi
  }else{
    df = rbind(df, dfi)
  }
}


names(df)[1:2]=c('gene','value')
#df$type = factor(df$type, levels = c('H3K27Me3','ATAC'))
df$group = ifelse(df$ID =='H721088','H3wt',
                  ifelse(df$ID %in% c('H164921','H177793'),'H3_thala','H3_pon'))
df$ID = factor(df$ID, levels = c("H164921","H177793","H160971","H170008","H172017","H179186","H182545","H721088"))
library(ggplot2)
ggplot(df, aes(x = gene, y = value))+
  geom_boxplot(aes(fill = HOX),outlier.shape = NA)+
  geom_jitter(width = 0.1,cex=0.1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+
  facet_grid(type~., scales = 'free_y')+
  #scale_fill_gradient2(low = 'blue', high = 'red',midpoint = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(df[df$type=='ATAC',], aes(x = gene, y = ID, fill = value))+
  geom_tile()+theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = .1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+geom_hline(yintercept = c(2.5,7.5),lty=3)+
  labs(title = 'ATACseq',y='',x='')
ggplot(df[df$type=='H3K27Me3',], aes(x = gene, y = ID, fill = log2(value+1)))+
  geom_tile()+theme(axis.text.x = element_text(angle = 90,hjust = 1)) +
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+geom_hline(yintercept = c(2.5,7.5),lty=3)+
  labs(title = 'H3K27Me3',y='',x='')

#write.table(df, file= 'spinalDMG.HOX.ATACH3K27me3.txt', row.names = F, quote = F, sep = "\t")
df = read.delim('spinalDMG.HOX.ATACH3K27me3.txt')
df$gene = factor(df$gene, levels = df$gene[df$type=='H3K27Me3'&df$ID=='H160971'])
df$tumorgroup = ifelse(df$ID =='H721088','H3wt',
                  ifelse(df$ID %in% c('H164921','H177793'),'H3_thala','H3_pon'))

ge  = read.delim('spinalcordDMG_scRNA.HOX.txt') # 
ge$gene = factor(ge$gene, levels = levels(df$gene))
names(ge)[6]='ID'
ge$type = 'scRNAseq'; ge$value = log2(ge$CPM+1)
ge$tumorgroup = ifelse(ge$Group=='groupA','H3_pon', ifelse(ge$Group=='groupB','H3_thala','others'))

ggplot(ge, aes(x = gene, y = value))+
  geom_boxplot(aes(fill = HOX),outlier.shape = NA)+
  geom_jitter(width = 0.1,cex=0.1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+
  #facet_grid(type~., scales = 'free_y')+
  #scale_fill_gradient2(low = 'blue', high = 'red',midpoint = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

df = rbind(df[,c('gene','value','HOX','ID','type','tumorgroup')], ge[,c('gene','value','HOX','ID','type','tumorgroup')])

df$type = factor(df$type, levels = c('scRNAseq','ATAC','H3K27Me3'))
ggplot(df[df$tumorgroup%in% c('H3_pon','H3_thala'),], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = HOX),outlier.shape = NA)+
  geom_jitter(width = 0.1,cex=0.1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+
  facet_grid(type~tumorgroup, scales = 'free_y')+
  #scale_fill_gradient2(low = 'blue', high = 'red',midpoint = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid = element_blank())

dfacmean = aggregate(value ~ gene, FUN = mean, data = df[df$type=='ATAC',])
dfk27me3mean = aggregate(value ~ gene, FUN = mean, data = df[df$type=='H3K27Me3',])


df$int = as.integer(substr(df$gene,5,7))
df = df[df$ID!='H721088',]
library(cowplot)
p10 = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXA',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))
  
p11 = ggplot(dfacmean[grepl('HOXA',dfacmean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p12=ggplot(dfk27me3mean[grepl('HOXA',dfk27me3mean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())

p20 = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXB',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p21 = ggplot(dfacmean[grepl('HOXB',dfacmean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p22=ggplot(dfk27me3mean[grepl('HOXB',dfk27me3mean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


p30 = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXC',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p31 = ggplot(dfacmean[grepl('HOXC',dfacmean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.2)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p32=ggplot(dfk27me3mean[grepl('HOXC',dfk27me3mean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())



p40 = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXD',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p41 = ggplot(dfacmean[grepl('HOXD',dfacmean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.2)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p42=ggplot(dfk27me3mean[grepl('HOXD',dfk27me3mean$gene),], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 2)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


plot_grid(p10,p20,p30,p40,
          p12,p22,p32,p42,
          NULL,NULL,NULL,NULL,
          p11,p21,p31,p41,
          ncol = 4,
          rel_heights = c(0.35,0.1,-0.05,0.1),
          rel_widths = c(1/4,1/4,1/4,1/4))

#separate group A and group B
dfacmean = aggregate(value ~ gene+tumorgroup, FUN = mean, data = df[df$type=='ATAC',])
dfk27me3mean = aggregate(value ~ gene+tumorgroup, FUN = mean, data = df[df$type=='H3K27Me3',])

p10a = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXA'&df$tumorgroup=='H3_pon',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p11a = ggplot(dfacmean[grepl('HOXA',dfacmean$gene)&dfacmean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p12a=ggplot(dfk27me3mean[grepl('HOXA',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())

p20a = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXB'&df$tumorgroup=='H3_pon',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p21a = ggplot(dfacmean[grepl('HOXB',dfacmean$gene)&dfacmean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p22a=ggplot(dfk27me3mean[grepl('HOXB',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = .3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


p30a = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXC'&df$tumorgroup=='H3_pon',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p31a = ggplot(dfacmean[grepl('HOXC',dfacmean$gene)&dfacmean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p32a=ggplot(dfk27me3mean[grepl('HOXC',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())



p40a = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXD'&df$tumorgroup=='H3_pon',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p41a = ggplot(dfacmean[grepl('HOXD',dfacmean$gene)&dfacmean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p42a=ggplot(dfk27me3mean[grepl('HOXD',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_pon',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = .3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


plot_grid(p10a,p20a,p30a,p40a,
          p12a,p22a,p32a,p42a,
          NULL,NULL,NULL,NULL,
          p11a,p21a,p31a,p41a,
          ncol = 4,
          rel_heights = c(0.35,0.1,-0.05,0.1),
          rel_widths = c(1/4,1/4,1/4,1/4))


p10b = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXA'&df$tumorgroup=='H3_thala',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p11b = ggplot(dfacmean[grepl('HOXA',dfacmean$gene)&dfacmean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p12b=ggplot(dfk27me3mean[grepl('HOXA',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())

p20b = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXB'&df$tumorgroup=='H3_thala',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p21b = ggplot(dfacmean[grepl('HOXB',dfacmean$gene)&dfacmean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p22b=ggplot(dfk27me3mean[grepl('HOXB',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = .3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


p30b = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXC'&df$tumorgroup=='H3_thala',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p31b = ggplot(dfacmean[grepl('HOXC',dfacmean$gene)&dfacmean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p32b=ggplot(dfk27me3mean[grepl('HOXC',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())



p40b = ggplot(df[df$type=='scRNAseq'&df$HOX=='HOXD'&df$tumorgroup=='H3_thala',], aes(x = gene, y = value))+
  geom_boxplot(aes(fill = int),outlier.shape = NA, show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  geom_jitter(width = 0.1,cex=1,aes(color=int),alpha=0.75)+scale_x_discrete(position = "top") +
  scale_color_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b',midpoint = 7)+
  theme_bw()+theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank(),
                   axis.text.x = element_text(angle=90,hjust=1))

p41b = ggplot(dfacmean[grepl('HOXD',dfacmean$gene)&dfacmean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = 0.1)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())
p42b=ggplot(dfk27me3mean[grepl('HOXD',dfk27me3mean$gene)&dfk27me3mean$tumorgroup=='H3_thala',], aes(x = gene, y = 1))+
  geom_tile(aes(fill = value),color = 'black',show.legend = F)+
  scale_fill_gradient2(low = '#f7fcfd',mid = '#8c96c6',high = '#4d004b', midpoint = .3)+
  theme_bw()+theme(axis.ticks = element_blank(),axis.text = element_blank(),rect = element_blank(),
                   panel.grid = element_blank(), axis.title = element_blank())


plot_grid(p10b,p20b,p30b,p40b,
          p12b,p22b,p32b,p42b,
          NULL,NULL,NULL,NULL,
          p11b,p21b,p31b,p41b,
          ncol = 4,
          rel_heights = c(0.35,0.1,-0.05,0.1),
          rel_widths = c(1/4,1/4,1/4,1/4))


#
gebulk = read.delim('spinalDMG.HOX.expression.txt')
gebulk = gebulk[gebulk$Group!='others',]
gebulk$gene = factor(gebulk$gene,levels = c(paste0('HOXA', 1:13),paste0('HOXB', 1:13),paste0('HOXC', 1:13),paste0('HOXD', 1:13)))
ggplot(gebulk, aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(aes(fill = HOX),outlier.shape = NA)+
  geom_jitter(width = 0.1,cex=0.1)+
  geom_vline(xintercept = c(11.5,21.5,30.5),lty=3)+
  #facet_grid(type~., scales = 'free_y')+
  #scale_fill_gradient2(low = 'blue', high = 'red',midpoint = 3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

