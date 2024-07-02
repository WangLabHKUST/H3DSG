rm(list=ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
source('./required_packages.R')
## Data load
mut_info <- (readRDS('./Four_omics_data/driver_SNV.RDS'))
cnv_info <- (readRDS('./Four_omics_data/cnv_info.RDS'))
meth_info <- (readRDS('./Four_omics_data/meth_info.RDS'))
RNA_info <- (readRDS('./Four_omics_data/rna_info.RDS'))
## Order by sample name
mut_info <- mut_info[order(rownames(mut_info)),]
cnv_info <- cnv_info[order(rownames(cnv_info)),]
RNA_info <- RNA_info[order(rownames(RNA_info)),]
meth_info <- meth_info[order(rownames(meth_info)),]
library(igraph)
library(SNFtool)
source('./internal.R')
## Distance calculation
Dist1 = dist2((as.matrix(mut_info)),(as.matrix(mut_info)))
Dist2 = dist2((as.matrix(cnv_info)),(as.matrix(cnv_info)))
Dist3 = dist2((as.matrix(RNA_info)),(as.matrix(RNA_info)))
Dist4 = dist2((as.matrix(meth_info)),(as.matrix(meth_info)))

## Parameter setting
K = 8 #Number of neighbors, usually (10~30)
alpha = 0.3 	# Hyperparameter, usually (0.3~0.8)
T1 = 15 # Number of Iterations, usually (10~20)

## Affinity matrix construction
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)
W4 = affinityMatrix(Dist4, K, alpha)
source('./cons_SNF.R')
### Label for each omics
lab <- c('D','D','R','E')
W_DSNF = contro_SNF(list(W1,W2,W3,W4),lab, K, T1)
### Number of cluster estimating based on integrated network 
cluster <- estimateNumberOfClustersGivenGraph(W_DSNF[[1]],NUMC=2:6)
print(cluster$`Eigen-gap best`)
Num_cluster <-cluster$`Eigen-gap best`
## Spectral clustering
group = spectralClustering(W_DSNF[[1]],Num_cluster)
displayClusters(W_DSNF[[1]], group)
length(rownames((as.matrix(Dist1)))[which(group==1)])
rownames((as.matrix(Dist1)))[which(group==1)]
length(rownames(t(as.matrix(Dist1)))[which(group==2)])
rownames(t(as.matrix(Dist1)))[which(group==2)]
source('./extract_subnetwork.R')
Split_net <- W_DSNF[[2]]
num_edg <- 1
som_net <- Split_net[[1]]
colnames(som_net) <- colnames(W1)
rownames(som_net) <- colnames(W1)
net_type <- "som"
net_som <- extract_subnetwork(som_net,num_edg,net_type)
cnv_net <- Split_net[[2]]
colnames(cnv_net) <- colnames(W1)
rownames(cnv_net) <- colnames(W1)
net_type <- "cnv"
net_cnv <- extract_subnetwork(cnv_net,num_edg,net_type)
rna_net <- Split_net[[3]]
colnames(rna_net) <- colnames(W1)
rownames(rna_net) <- colnames(W1)
net_type <- "rna"
net_rna <- extract_subnetwork(rna_net,num_edg,net_type)
meth_net <- Split_net[[4]]
colnames(meth_net) <- colnames(W1)
rownames(meth_net) <- colnames(W1)
net_type <- "meth"
net_meth <- extract_subnetwork(meth_net,num_edg,net_type)
## Combine_all_network_type
net_combine <- rbind(rbind(net_som,net_cnv),
                     rbind(net_rna,net_meth))
## Integrated_omics
net_combine_unique <- rbind()
for(i in unique(net_combine[,1]))
{
  k1 <- which(net_combine[,1] %in% i)
  net_sub <- net_combine[k1,]
  for(j in unique(net_sub[,2]))
  {
    k2 <- which(net_sub[,2] %in% j)
    net_combine_unique <- rbind(net_combine_unique,
                            c(i,j,length(k2)))
  }
}
colnames(net_combine_unique) <- c("Source","Targets","Weight")
### Store the combined network
#write.table(net_combine_unique,'net_combine_unique.txt',
#            sep="\t",row.names = F)