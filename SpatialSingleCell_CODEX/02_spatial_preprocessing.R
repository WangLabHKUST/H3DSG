library(imcRtools)
library(cytomapper)
library(lisaClust)
library(pheatmap)
library(viridis)


library(ggplot2)
library(viridis)

out_path = r"(D:\OneDrive - HKUST Connect\Work\SingleCell\Multiplex-IHC\spatial_analysis\0623)"
dir.create(out_path, showWarnings = FALSE)
setwd(out_path)
name_lst = c("X158665","X158785","X161865",'X164398',"X164921","X165114",'X179186','X179980')#
# all_df = read.csv(file.path('csv','all_df_0611.csv'))
celltype_col_name = 'celltype_agg'
for (name in name_lst) {
  if(file.exists(file.path('objects',sprintf('%s.rds',name))) ){
    print(paste("Skipping", name))
    next
  }
  
  print(name)
  cur_features = read.csv(file.path('csv',sprintf('%s_unique_tumor_renamed_0623.csv',name)))
  colnames(cur_features)
  dim(cur_features)
  cur_features = cur_features[cur_features[celltype_col_name] != 'Others',]
  table(cur_features[celltype_col_name])
counts <- cur_features[,!colnames(cur_features) %in% c(celltype_col_name)]

meta <- cur_features[,c(celltype_col_name)]
head(meta)
x_name = 'Centroid_X_µm'
y_name = 'Centroid_Y_µm'
coords <- cur_features[,c(x_name, y_name)]
spe <- SpatialExperiment(assays = list(counts = t(counts)),
                          colData = meta, 
                          sample_id = name,
                          spatialCoords = as.matrix(coords))

spe$coord_x = spatialCoords(spe)[,x_name]
spe$coord_y = spatialCoords(spe)[,y_name]
colnames(spe) = seq(1,ncol(spe))
spe_bak = spe
print("Start buildSpatialGraph")
spe <- buildSpatialGraph(spe, img_id = "sample_id", type = "expansion", threshold = 40,coords=c(x_name,y_name))

set.seed(722)
col_pair_name = colPairNames(spe)
print("Finished detectCommunity")
table(spe$X)
spe <- aggregateNeighbors(spe, 
                          colPairName = col_pair_name, 
                          aggregate_by = "metadata", 
                          count_by = "X")

saveRDS(spe, file.path(out_path,'objects',paste0(name,'.rds')))

}






