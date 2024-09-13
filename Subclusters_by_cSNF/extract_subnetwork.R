extract_subnetwork <- function(som_net,num_edg,net_type){
if (!require('igraph')) 
  install.packages('igraph')
som_net_for_visual <- rbind()
for(i in 1:dim(som_net)[1])
{
  for(j in 1:dim(som_net)[2])
  {
    if(i!=j){
      som_net_for_visual<- rbind(som_net_for_visual,
                                 c(rownames(som_net)[i],
                                   colnames(som_net)[j],som_net[i,j]))
    }
  }
}
dim(som_net_for_visual)

som_net_weight_sel <- rbind()
for(i in unique(som_net_for_visual[,1]))
{
  k1 <- which(som_net_for_visual[,1] %in% i)
  if(length(k1)>0)
  {
    used_net <- som_net_for_visual[k1,]
    used_net <- used_net[order(as.numeric(used_net[,3]),
                               decreasing = T),]
    som_net_weight_sel <- rbind(som_net_weight_sel,
                                used_net[c(1:num_edg),])
  }
}
dim(som_net_weight_sel)
library(igraph)
som_net_vis_<- graph_from_data_frame(som_net_weight_sel[,c(1,2)],directed = F)
som_net_vis_<- simplify(som_net_vis_,remove.multiple = T,remove.loops =T)
som_net_vis_edge <- get.edgelist(som_net_vis_)
dim(som_net_vis_edge)
som_net_weight <- rbind()
for(i in 1:dim(som_net_vis_edge)[1])
{
  k1 <- which(som_net_weight_sel[,1]%in% som_net_vis_edge[i,1])
  k2 <- which(som_net_weight_sel[,2]%in% som_net_vis_edge[i,2])
  k3 <- intersect(k1,k2)
  if(length(k3)>0)
  {
    som_net_weight <- rbind(som_net_weight,
                            c(som_net_vis_edge[i,1],som_net_vis_edge[i,2],
                              som_net_weight_sel[k3,3],net_type
                              #"meth"))#"rna"))#"cnv"))#"som"
                            ))
  }else{
    k1 <- which(som_net_weight_sel[,1]%in% som_net_vis_edge[i,2])
    k2 <- which(som_net_weight_sel[,2]%in% som_net_vis_edge[i,1])
    k3 <- intersect(k1,k2)
    if(length(k3)>0)
    {
      som_net_weight <- rbind(som_net_weight,
                              c(som_net_vis_edge[i,1],som_net_vis_edge[i,2],
                                som_net_weight_sel[k3,3],net_type
                                
                              ))
    }
  }
}
dim(som_net_weight)
som_net_weight <- as.data.frame(som_net_weight)
colnames(som_net_weight) <- c("Source","Target","Weight","Type")
som_net_weight <- as.data.frame(som_net_weight)
#som_net_weight$new_source <- rep("",length(som_net_weight$Source))
#som_net_weight$new_target <- rep("",length(som_net_weight$Target))
return(som_net_weight)
}