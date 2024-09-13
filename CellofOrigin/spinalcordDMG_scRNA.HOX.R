setwd('~/Dropbox/SpinalCordGlioma/scRNAseq/')
fds = list.files(pattern = 'bc_matrix')
for (fd in fds){
  print(fd)
  cellranger_matrix <- Seurat::Read10X(fd)
  
  seurat_prefilt <- CreateSeuratObject(
    counts       = cellranger_matrix,
    min.cells    =0,    # specifies the min number of cells in which each gene must be detected
    min.features = 200, # specifies the min number of genes which must be detected in each cell
    project      = 'test'            # populates the project.name slot in the Seurat object
  )
  gns = rownames(seurat_prefilt@assays$RNA@counts)
  gns = gns[grepl('^HOX',gns)]
  gns = gns[!grepl('-',gns)]
  
  as.matrix(seurat_prefilt@assays$RNA@data[gns,])->dt
  df = as.data.frame(apply(dt, 1, sum))
  df$HOX = substr(rownames(df),1,4)
  df$int = as.integer(substr(rownames(df),5,7))
  names(df)[1]='count'
  df = df[order(df$HOX, df$int),]
  df$gene = rownames(df)
  df$gene = factor(df$gene, levels = df$gene)
  df$CPM = 10^6*df$count/sum(seurat_prefilt@assays$RNA@counts)
  df$SID = fd
  if(fd==fds[1]){
    res = df
  }else{
    res = rbind(res, df)
  }
}

res$SID= gsub("\\.filtered_feature_bc_matrix",'',res$SID)
g1 = c('SP04','SP15','SP24','SP41','X172017','SP204')
g2 = c('SP05','SP11','SP13','SP20T5','SP21T5','SP22','X177793','X747347')
res$Group = ifelse(res$SID %in% g1, 'groupA',
                   ifelse(res$SID %in% g2, 'groupB','others'))
write.table(res, file = 'spinalcordDMG_scRNA.HOX.txt', row.names = F, quote =F, sep = "\t")


#read in and plot
res = read.delim('spinalcordDMG_scRNA.HOX.txt')
res = res[order(res$HOX, res$int),]
res$gene = factor(res$gene, levels = res$gene[!duplicated(res$gene)])



library(ggplot2)
library(cowplot)
p1 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXA',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  #facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p2 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXB',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  #facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p3 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXC',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  #facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p4 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXD',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  #facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))

plot_grid(p1,p2,p3,p4,ncol=4)

p1 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXA',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p2 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXB',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p3 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXC',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
p4 = ggplot(res[res$Group%in% c('groupA','groupB')&res$HOX=='HOXD',], aes(x = gene, y = log2(CPM+1)))+
  geom_boxplot(width = 0.5, aes(color = int,fill=int),outlier.shape = NA, show.legend = F,alpha=0.5)+
  geom_jitter(aes(color = int),pch=20,cex=0.1,width=0.1, show.legend = F)+
  facet_wrap(~Group,ncol = 1)+
  scale_color_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  scale_fill_gradient2(high = '#f7fcfd',mid = '#8c96c6',low = '#4d004b',midpoint = 7)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text  = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))

plot_grid(p1,p2,p3,p4,ncol=4)


