setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/HOX/')

ge = read.delim('CGGA.H3K27M.hg38.expression.RPKM.txt',row.names = 1)
df = read.delim('Spinalcord.H3K27M.expression.RPKM.txt')
df.cl = read.delim('H3K27altered_project_metadata_updated.tsv')

df.pons = df[,names(df) %in% paste0('P',df.cl$ID2[df.cl$TSNE.cluster=='A'])]
df.thalamus = df[,names(df) %in% paste0('P',df.cl$ID2[df.cl$TSNE.cluster=='B'])]

ge.thalamus = ge[,names(ge) %in% c("CGGA_791","CGGA_P266","CGGA_1185","CGGA_1320","CGGA_1735")]
shrdGns = unique(intersect(rownames(df.pons), rownames(ge.thalamus)))
deg.thalamus= data.frame(gene = shrdGns, sp = 0, cg = 0 , pval=1)
for (i in 1:nrow(deg.thalamus)){
  x = as.numeric(df.thalamus[deg.thalamus$gene[i],]) #spine thalamus-like
  y = as.numeric(ge.thalamus[deg.thalamus$gene[i],]) #brain thalamus
  deg.thalamus$sp[i] = median(x)
  deg.thalamus$cg[i] = median(y)
  deg.thalamus$pval[i] = wilcox.test(x,y)$p.value
}
deg.thalamus$lfc = log2(deg.thalamus$sp+1) - log2(deg.thalamus$cg+1)
gns = read.delim('~/Dropbox/AYAglioma/AYA_H3F3A/data_summary/Homo_sapiens.GRCh37.75.chr.genes.txt',header = F)
deg.thalamus = deg.thalamus[deg.thalamus$gene %in% gns$V7[gns$V2=='protein_coding'],]
deg.thalamus$lab = ifelse(deg.thalamus$pval<0.01 & deg.thalamus$lfc>1&startsWith(deg.thalamus$gene,'HOX'),deg.thalamus$gene,NA)
deg.thalamus$type = ifelse(deg.thalamus$pval<0.01 & deg.thalamus$lfc>1, 'up',
                  ifelse(deg.thalamus$pval<0.01 & deg.thalamus$lfc< -1, 'dn','nosig'))
library(ggplot2)
library(ggrepel)
ggplot(deg.thalamus, aes(x = lfc, y = -log10(pval), label = lab))+
  geom_point(aes(color = type,size= type))+
  scale_color_manual(values = c('blue','#cccccc','red'))+
  scale_size_manual(values = c(1,0.2,1))+
  geom_text_repel(segment.size=0.25,size=3)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()

ge.pons = read.delim('StJude.DIPG.RPKM.txt')
ge.pons = ge.pons[!is.na(ge.pons$gene),]
shrdGns = unique(intersect(rownames(df.pons), ge.pons$gene))
deg.pons= data.frame(gene = shrdGns, sp = 0, cg = 0 , pval=1)

for (i in 1:nrow(deg.pons)){
  gni = deg.pons$gene[i]
  x = as.numeric(df.pons[gni,]) #spine thalamus-like
  tmpy = ge.pons[ge.pons$gene==gni,]
  tmpy = tmpy[1,]
  y = as.numeric(tmpy[,1:5]) #brain pons
  deg.pons$sp[i] = median(x)
  deg.pons$cg[i] = median(y)
  deg.pons$pval[i] = wilcox.test(x,y)$p.value
}
deg.pons$lfc = log2(deg.pons$sp+1) - log2(deg.pons$cg+1)
deg.pons = deg.pons[deg.pons$gene %in% gns$V7[gns$V2=='protein_coding'],]

#
mg = data.frame(gene = intersect(deg.pons$gene, deg.thalamus$gene))
mg$logFC_pons = deg.pons$lfc[match(mg$gene, deg.pons$gene)]
mg$logFC_thalamus = deg.thalamus$lfc[match(mg$gene, deg.thalamus$gene)]
mg$lab = ifelse(startsWith(mg$gene,'HOX'), mg$gene,NA)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
mg$density = get_density(mg$logFC_pons, mg$logFC_thalamus, n = 100)

ggplot()+
  geom_point(data = mg[!(mg$logFC_thalamus==0 & mg$logFC_pons< -5),],show.legend = F,
             mapping=aes(color = density,x = logFC_pons, y = logFC_thalamus),alpha=0.5,size=0.1)+
  geom_point(mg[!is.na(mg$lab),],mapping=aes(x = logFC_pons, y = logFC_thalamus),size=0.5,col='red')+
  geom_vline(xintercept = 0,lty=2,col='#dddddd')+
  geom_hline(yintercept = 0,lty=2,col='#dddddd')+
  geom_abline(slope = 1,intercept  = 0,lty=2,col='#dddddd')+
  #lims(x = c(-4,4), y =c(-4,4))+
  scale_x_continuous(limits = c(-4,4),expand = c(0,0))+
  scale_y_continuous(limits = c(-4,4),expand = c(0,0))+
  scale_color_viridis_c()+
  geom_text_repel(data = mg[(mg$logFC_pons>1&mg$logFC_pons<3.5)| (mg$logFC_thalamus>1&mg$logFC_thalamus<3.5),], mapping = aes(x = logFC_pons, y = logFC_thalamus,label = lab),
                  bg.color = "white", bg.r = 0.25,
                  col='red',min.segment.length = 0,segment.size=0.25,size=1.5,max.overlaps = 20)+
  theme_classic()+theme(axis.text = element_text(color = 'black'))+
  labs(x = ' ', y = ' ')

#plot HOXA2 and A9 as an example to show the values in each sample
gni = 'HOXA2'
df_hoxd8 = data.frame(Gene = c(as.numeric(df.thalamus[gni,]),
                                as.numeric(df.pons[gni,]),
                                as.numeric(ge.thalamus[gni,c('CGGA_1320', 'CGGA_1735', 'CGGA_791')]),
                                as.numeric(ge.pons[ge.pons$gene==gni,1:5])),
                      Group = c(rep('Thalamus-like spinal DMG',ncol(df.thalamus)),
                                rep('Pons-like spinal DMG',ncol(df.pons)),
                                rep('Thalamus DMG',3),
                                rep('Pons DMG',5))
                                )
df_hoxd8$Group=factor(df_hoxd8$Group, levels = c('Thalamus DMG','Pons DMG','Thalamus-like spinal DMG','Pons-like spinal DMG'))
library(ggpubr)
library(ggbeeswarm)
ggplot(df_hoxd8, aes(x = Group, y = Gene))+
  geom_boxplot(aes(fill = Group),outlier.shape = NA,show.legend = T )+
  scale_fill_manual(values = c('#d1e5f0',"#ffd92f",'#377eb8','#ff7f00'))+
  geom_quasirandom(width = 0.15, alpha=0.7)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons DMG','Thalamus DMG'),
                                                            c('Pons-like spinal DMG','Thalamus-like spinal DMG'),
                                                            c('Pons-like spinal DMG','Pons DMG'),
                                                            c('Thalamus DMG','Thalamus-like spinal DMG') ))+
  theme_classic()+labs(x ='', y = paste(gni,' (RPKM)'))+
  theme(axis.text.x = element_blank())

ggplot(df_hoxd8, aes(x = Group, y = log2(Gene+1)))+
  geom_boxplot(aes(fill = Group),outlier.shape = NA,show.legend = T )+
  scale_fill_manual(values = c('#d1e5f0',"#ffd92f",'#377eb8','#ff7f00'))+
  geom_quasirandom(width = 0.15, alpha=0.7)+
  stat_compare_means(label = 'p.format', comparisons = list(c('Pons DMG','Thalamus DMG'),
                                                            c('Pons-like spinal DMG','Thalamus-like spinal DMG'),
                                                            c('Pons-like spinal DMG','Pons DMG'),
                                                            c('Thalamus DMG','Thalamus-like spinal DMG') ))+
  theme_classic()+labs(x ='', y = paste(gni,' (log2 RPKM)'))+
  theme(axis.text.x = element_blank())


gns = ge.pons$gene[grepl('HOX[ABCD]',ge.pons$gene)]
df_hox_all = data.frame(Gene = c(as.numeric(df.thalamus[gni,]),
                               as.numeric(df.pons[gni,]),
                               as.numeric(ge.thalamus[gni,c('CGGA_1320', 'CGGA_1735', 'CGGA_791')]),
                               as.numeric(ge.pons[ge.pons$gene==gni,1:5])),
                      Group = c(rep('Thalamus-like spinal DMG',ncol(df.thalamus)),
                                rep('Pons-like spinal DMG',ncol(df.pons)),
                                rep('Thalamus DMG',3),
                                rep('Pons DMG',5))
)
hox.thalamus_like = df.thalamus[gns,]
hox.pons_like = df.pons[gns,]
hox.thalamus_brain = ge.thalamus[gns,c('CGGA_1320', 'CGGA_1735', 'CGGA_791')]
hox.pons_brain = ge.pons[ge.pons$gene %in% gns,1:6]
rownames(hox.pons_brain) = hox.pons_brain$gene
hoxgns = c(paste0('HOXA',c(1:7,9:11,13)),paste0('HOXB',c(1:9,13)),
           paste0('HOXC',c(4:6,8:13)),paste0('HOXD',c(1,3,4,8:13)))
hox.thalamus_like = hox.thalamus_like[hoxgns,]
hox.pons_like = hox.pons_like[hoxgns,]
hox.thalamus_brain = hox.thalamus_brain[hoxgns,]
hox.pons_brain = hox.pons_brain[hoxgns,c(4,5)]
library(pheatmap)
m = cbind(hox.thalamus_brain,hox.thalamus_like,hox.pons_brain, hox.pons_like )

m_ann = data.frame(row.names = names(m), type = c(rep('Brain thalamus',3),
                                                  rep('Thalamus-like',24),
                                                  rep('Brain pons',2),
                                                  rep('Pons-like',16)))
pheatmap(log2(m+1),color = colorRampPalette(colors = c("cornflowerblue","#ffff90",  "yellow", "red"))(50),
         cluster_rows = F, gaps_row = c(11,21,30),border_color = '#dddddd',
         annotation_col = m_ann, clustering_method = 'single',show_colnames = T,
         fontsize_row = 8, treeheight_col = 10,#gaps_col = c(2,5),
         annotation_colors = list(type = c('Brain pons'='#f1a23a','Brain thalamus'='#95add2','Pons-like'='#86180f','Thalamus-like'='#3477b0')))

tumorlocation = read.delim('SP45_ID_location_detail.txt')
tumorlocation$id = paste0('P',tumorlocation$ID2)
m_ann$location = tumorlocation$Location[match(rownames(m_ann), tumorlocation$id)]
m_ann$location2 = tumorlocation$Location_2[match(rownames(m_ann), tumorlocation$id)]
m_ann$location[m_ann$type=='Brain pons']='Brain_pons'
m_ann$location[m_ann$type=='Brain thalamus']='Brain_thalamus'
m_ann$location2[m_ann$type=='Brain pons']='Brain_pons'
m_ann$location2[m_ann$type=='Brain thalamus']='Brain_thalamus'
m_ann$location2 = factor(m_ann$location2,levels = c('Brain_thalamus','Brain_pons',
                                                    'M-C','C','C-T','M-T1','T','T-L'))
m_ann = m_ann[order(m_ann$location2,m_ann$location,row.names(m_ann)),]
m_ann$idx = 1:nrow(m_ann)
m_ann = m_ann[c(1:14,17,18,15:16,19:45),]
names(m_ann)[2:3]=c('location2','location')
pheatmap(log2(m[,rownames(m_ann)]+1),color = colorRampPalette(colors = c("darkblue",'#ffff90',"#ffff50",  "yellow",'orange', "red"))(50),
         cluster_rows = F, cluster_cols = F,show_colnames = T,
         gaps_row = c(11,21,30),border_color = '#dddddd',legend = T,
         annotation_col = subset(m_ann,select = c('type','location')),# clustering_method = 'single',show_colnames = F,
         fontsize_row = 6, treeheight_col = 10,gaps_col = c(3,5),
         annotation_colors = list(type = c('Brain pons'='#f1a23a','Brain thalamus'='#95add2','Pons-like'='#86180f','Thalamus-like'='#3477b0')))
