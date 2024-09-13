library(imcRtools)
library(cytomapper)

library(pheatmap)
library(viridis)
library(lisaClust)


library(ggplot2)
library(viridis)
library(plyr)
library( gtools) 
library(dplyr)
library(rstatix)
library(circlize)
out_path = r"(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\0523)"
meta = read.csv(r'(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\objects\meta.csv)')
rownames(meta) = meta$ID
### make dir
obj_dir = 'objects_filtered_final'
heatmap_dir = 'comb_filtered_new_type_final'
dir.create(file.path(out_path), showWarnings = FALSE)
dir.create(file.path(out_path,obj_dir), showWarnings = FALSE)
dir.create(file.path(out_path,'vis'), showWarnings = FALSE)
dir.create(file.path(out_path,'community'), showWarnings = FALSE)
setwd(out_path)
name_lst = c("X158665","X158785","X161865",'X164398',"X164921","X165114",'X179186','X179980')#
group_lst = c('')

x_name = 'Centroid.X.µm'
y_name = 'Centroid.Y.µm'

spe_all = list()
for (name in name_lst) {
  spe_all[[name]] = readRDS (file.path(out_path,obj_dir,paste0(name,'_clustered_new_celltype.rds'))) ##_clustered_new_celltype.rds
}
for(name in name_lst){
  print(name)
  print(dim(spe_all[[name]]))
}
all_df = data.frame()
for (name in name_lst) {
  spe = spe_all[[name]]
  prev_len = dim(all_df)[1]
  this_df = spe$aggregatedNeighbors
  this_df$X = spe$X
  this_df$name = name
  this_df$group = meta[name,'Group']
  this_df$x = spe$coord_x
  this_df$y = spe$coord_y

  for (n_cl in c(3,4,5,6,7,8,10,12,14,16,18)) {
    this_df[[paste0('cluster_',n_cl)]] = spe[[paste0('cluster_',n_cl)]]
  }
  all_df = smartbind(all_df, this_df)
  later_len = dim(all_df)[1]
}
head(all_df)
tail(all_df)
all_df[is.na(all_df)] = 0
lst = sapply(rownames(all_df), function(x) {
  split_element <- strsplit(x, ":")
  split_element[[1]][2]
})

rownames(all_df) = paste0(all_df$name,"_",lst)
rownames(all_df)[1]
write.csv(all_df, file.path(obj_dir,'all_df.csv'))
all_df_clustering = subset(all_df, select = -c(name,X))
colnames(all_df)
all_df$Group = meta[all_df$name,'Group']
df_used = all_df
out_folder_name = 'final'
head(df_used)
for (n_cl in c(3,4,5,6,7,8,10,12,14,16,18)) {
  for_plot_bak <- table(df_used[[paste0('cluster_',n_cl)]], df_used$X)
  for_plot = for_plot_bak
  
  p = pheatmap(for_plot,
               color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
               show_rownames = T, cluster_cols = FALSE, cluster_rows = F,
               scale = "column", number_format = "%.0f")

  ggsave(file.path(heatmap_dir,out_folder_name,paste0(n_cl,"_kmeans.jpg")),p, width = 10, height = n_cl*0.5)
  ggsave(file.path(heatmap_dir,out_folder_name,paste0('cluster_',n_cl),paste0(n_cl,'_kmeans.jpg')), p, width = 10, height = n_cl*0.5)
  ##### plot group comp #####
  # for_plot_Group <- data.frame(prop.table(table(df_used[[paste0('cluster_',n_cl)]], df_used$name)))
  for_plot_Group <- table(df_used[[paste0('cluster_',n_cl)]], df_used$name)
  # for_plot_Group = for_plot_Group / colsum(for_plot_Group)
  for_plot_Group <- t(t(for_plot_Group) / colSums(for_plot_Group))
  for_plot_Group = as.data.frame(for_plot_Group)
  colnames(for_plot_Group) = c('Niche','Patient','Proportion')
  for_plot_Group$Group = meta[for_plot_Group$Patient,'Group']
  for_plot_Group$Niche = as.numeric(for_plot_Group$Niche)
  for_plot_Group$Group = factor(for_plot_Group$Group)
  for_plot_Group$Group_num = as.numeric(for_plot_Group$Group)
  for_plot_Group$Proportion = as.numeric(for_plot_Group$Proportion)
  ggplot(for_plot_Group, aes(x=Niche, y=Proportion, fill=Group)) +
    geom_boxplot() +
    facet_wrap(~Niche, scale="free") +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black",  position = position_dodge(width = 0.75)) +

    ggtitle(paste("Num cluster:",n_cl)) + theme_classic()+
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background  = element_blank()) + scale_fill_manual(values=c('#cc4c02',"#0570b0"))

  ggsave(file.path(out_path,heatmap_dir,out_folder_name,paste0('comp_',n_cl,'.pdf')), width = 8, height = 6, dpi = 400)
  ####### plot by enrichment #####
  ## get the median prop of each Niche for each group
  median_df = data.frame()
  median_ratio_df = data.frame()
  for_plot_Group = for_plot_Group[order(for_plot_Group$Niche),]
  for (n in unique(for_plot_Group$Niche)) {
    for (g in unique(for_plot_Group$Group)) {
      median_prop = median(for_plot_Group[for_plot_Group$Niche == n & for_plot_Group$Group == g,]$Proportion)
      median_df = smartbind(median_df, data.frame(Niche = n, Group = g, median_prop))
    }
    median_ratio = median_df[median_df$Niche == n & median_df$Group=='Thalamus-like',]$median_prop / median_df[median_df$Niche == n & median_df$Group=='Pons-like',]$median_prop
    median_ratio_df = smartbind(median_ratio_df, data.frame(Niche = n, median_ratio))
  }

  rownames(median_ratio_df) = median_ratio_df$Niche
  rownames(for_plot) = as.numeric(rownames(for_plot))


  row_anno = as.data.frame(median_ratio_df$median_ratio)

  row_anno$Niches = rownames(row_anno)
  colnames(row_anno) = 'Enrichment'
  row_anno[row_anno$Enrichment > 2, 'Enrichment'] = 2

  index_order = order(row_anno$Enrichment)
  row_anno$Niche = rownames(row_anno)
  celltype_order = c('PPR','MTC','GPM','Myeloid','Neutrophil','T','NK','B','Endothelial','OPC','Neuron','Astrocyte')
  for_plot = for_plot[,celltype_order]
  for_plot = for_plot[sort(index_order),]
  for_plot = for_plot[index_order,]


  col_fun = colorRamp2(c(min(row_anno$Enrichment), 1, max(row_anno$Enrichment)), c("#FF3300", "white", "#2b8cbe"), space = "LAB" )
  row_anno_color = col_fun(row_anno[index_order,'Enrichment'])
  p = pheatmap(for_plot,
               color = colorRampPalette(c("dark blue", "white", "dark red"))(100),
               show_rownames = T, cluster_cols = FALSE, cluster_rows = F,
               scale = "column", number_format = "%.0f",
               gaps_col = c(3,8,9),
               annotation_names_col=T,
               annotation_row = data.frame(Enrichment = row_anno$Enrichment),
               annotation_colors = list(Enrichment = row_anno_color)
               )
  ggsave(file.path(out_path,heatmap_dir,out_folder_name,paste0('comp_',n_cl,'_enrichment.pdf')), p,
         width = 5, height = 3, dpi = 300)
  ggsave(file.path(out_path,heatmap_dir,out_folder_name,paste0('comp_',n_cl,'_enrichment.png')), p,
         width = 5, height = 3, dpi = 300)
}

##### test significance #####
df_used = all_df_filtered

for (n_cl in c(3,4,5,6,7,8,10,12,14,16,18)) {
  n_cl = 10
  # for_plot <- table(all_df_filtered[[paste0('cluster_',n_cl)]], all_df_filtered$X)
  for_plot <- table(df_used[[paste0('cluster_',n_cl)]], df_used$X)
  for_plot_Group <- data.frame(prop.table(table(df_used[[paste0('cluster_',n_cl)]], df_used$name)))
  colnames(for_plot_Group) = c('Niche','Patient','Proportion')
  for_plot_Group$Group = meta[for_plot_Group$Patient,'Group']
  for_plot_Group$Niche = as.numeric(for_plot_Group$Niche)
  for_plot_Group$Group = factor(for_plot_Group$Group)
  for_plot_Group$Group_num = as.numeric(for_plot_Group$Group)
  for_plot_Group$Proportion = as.numeric(for_plot_Group$Proportion)
  #### test for every Niche whether there is a significant difference between groups ####
  # group by Niche, use wilcox test to test the difference between groups
  wilc_test_res = data.frame()
  for (n in unique(for_plot_Group$Niche)) {
    test_res = wilcox.test(Proportion ~ Group, data = for_plot_Group[for_plot_Group$Niche == n,], exact=F)
    test_res_df = data.frame(Niche = n, test_res$statistic, test_res$p.value)
    wilc_test_res = smartbind(wilc_test_res,test_res_df)
  }
  wilc_test_res
  write.csv(wilc_test_res, file.path(out_path,heatmap_dir,out_folder_name,paste0('comp_',n_cl,'_wilcox_test.csv')))
}
