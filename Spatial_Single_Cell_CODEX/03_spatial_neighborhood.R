library(imcRtools)
library(cytomapper)

library(pheatmap)
library(viridis)
library(lisaClust)


library(ggplot2)
library(viridis)
library(plyr)
library( gtools) 
out_path = r"(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\0623)"


setwd(out_path)
name_lst = c("X158665","X158785","X161865",'X164398',"X164921","X165114",'X179186','X179980')#
group_lst = c('')


x_name = 'Centroid_X_µm'
y_name = 'Centroid_Y_µm'

spe_all = list()
for (name in name_lst) {
  spe_all[name] = readRDS (file.path(out_path,'objects',paste0(name,'_clustered.rds')))
}
unique(spe_all[name][[name]]$X)
st_color_list = c('GPM'='#b62526','MTC'='#e97550','PPR'='#fb9e4f','B'='#439030','Lymphocyte'='#2f609f','NK'='#2f609f',
                  'Myeloid'='#57b894','Neutrophil'='#97d68f','Neuron'='#00a6ff', 'OPC'='#5dc6ff',
                  'Astrocyte' ='#efc3c5', 'Vascular'='#b0aad1')

#### plot celltype spatial for every sample
for (name in name_lst) {
  
  p = plotSpatial(spe_all[name][[name]],
                  img_id = "sample_id",
                  node_size_fix = 0.00005,coords=c(x_name,y_name), node_color_fix=st_color_list[spe_all[name][[name]]$X]) 
  ggsave(file.path('vis',paste0(name,"_overall_major.pdf")),dpi = 300, width = 45, height = 30,p)
}
##### get dataframe for clustering
all_df = data.frame()
for (name in name_lst) {
  # name = "179186" #165114
  spe = spe_all[name][[name]]
  prev_len = dim(all_df)[1]
  this_df = spe$aggregatedNeighbors
  this_df$X = spe$X
  this_df$name = name
  all_df = smartbind(all_df, this_df)
  later_len = dim(all_df)[1]
  this_df$group = spe$group
  
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
all_df_clustering = subset(all_df, select = -c(name,X))
colnames(all_df_clustering)
rds = 722
for (n_cl in c(3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) {
  set.seed(rds)
  if (!dir.exists(file.path('comb',paste0('cluster_',n_cl,'_',rds)))) {
    dir.create(file.path('comb',paste0('cluster_',n_cl,'_',rds)))
  }
  cn_1 <- kmeans(all_df_clustering, centers = n_cl)
  all_df$cluster = cn_1$cluster
  labels_cn = data.frame(cluster=cn_1$cluster)

  labels_cn$name = sapply(rownames(labels_cn), function(x) {
    split_element <- strsplit(x, "_")
    split_element[[1]][1]
  })
  all_df$cluster = as.factor(all_df$cluster)
  for (name in name_lst) {
    spe = spe_all[name][[name]]
    spe[[paste0('cluster_',n_cl,'_',rds)]] = as.factor(labels_cn[labels_cn$name == name,]$cluster)
    spe[[paste0('cluster_',n_cl,'_',rds)]] = as.factor(spe[[paste0('cluster_',n_cl,'_',rds)]])
    print(table(spe[[paste0('cluster_',n_cl,'_',rds)]]))
    spe_all[name][[name]] = spe
  }
  
}

for(name in name_lst){
  saveRDS(spe_all[name][[name]], file.path('objects',paste0(name,'_clustered.rds')))
}
