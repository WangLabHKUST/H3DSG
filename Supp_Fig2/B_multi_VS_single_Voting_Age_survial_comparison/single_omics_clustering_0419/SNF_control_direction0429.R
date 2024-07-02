rm(list=ls())
setwd('/Users/shanghaixia/Desktop/glioma_prognosis_method_multi-omics/Multi_omics_SP_related')
library(readxl)
integrated_label <- read_excel('Integrated_table.xlsx',sheet=1)
integrated_36 <- integrated_label[which(integrated_label$with_DNA == "Y"),]
mut_list <- read.table('Mutation_list_show.txt',header=F,sep="\t")
Driver_mutations <- read_excel('./SP_WES_results/savi_hg19_Snp_results/pDfilter/SP_dirver_mutation.xlsx',
                               sheet=1)
Driver_mutations <- Driver_mutations[which(as.numeric(Driver_mutations$primary_score)<10),]
library(stringr)

Mut_info <- matrix(0,nrow=length(unique(Driver_mutations$id_new)),
                   ncol=length(mut_list$V1))
rownames(Mut_info) <- paste('P',unique(Driver_mutations$id_new),sep="")
colnames(Mut_info) <- mut_list$V1
for(i in 1:dim(Mut_info)[1])
{
  for(j in 1:dim(Mut_info)[2])
  {
    k1 <- which(paste('P',Driver_mutations$id_new,sep="") %in% rownames(Mut_info)[i])
    k2 <- which(Driver_mutations$V54 %in% colnames(Mut_info)[j])
    if(length(intersect(k1,k2))>0)
    {
      Mut_info[i,j] <- 1
    }
  }
}
#apply(Mut_info,2,sum)

used_group <-cbind(integrated_36$ID2...2,
                   cbind(integrated_36$`TSNE cluster...3`,
                         integrated_36$final_dec))
used_group[,1] <- paste('P',used_group[,1],sep="")
#write.table(used_group,'./PPSN_multi_omics/used_group_info.txt',sep="\t",row.names = F,quote=F)

## CNV
sp_cnv <- readRDS('SP_36_CNV_cytoband_level.RDS')
cnv_var <- apply(sp_cnv,1,var)
cnv_select_var <- names(sort(cnv_var, decreasing=TRUE))[1:520]
used_cnv <- sp_cnv[which(rownames(sp_cnv) %in% cnv_select_var),]

## RNA
sp_RNA <- readRDS('SP_36_RNA.RDS')
RNA_var <- apply(sp_RNA,1,var)
RNA_select_var <- names(sort(RNA_var, decreasing=TRUE))[1:3000]
used_RNA <- sp_RNA[which(rownames(sp_RNA) %in% RNA_select_var),]

## methy
sp_Meth <- readRDS('SP_36_DNA_methylation.RDS')
colnames(sp_Meth) <- paste('P',colnames(sp_Meth),sep="")
Meth_var <- apply(sp_Meth,1,var)
Meth_select_var <- names(sort(Meth_var, decreasing=TRUE))[1:3000]#3000
used_Meth <- sp_Meth[which(rownames(sp_Meth) %in% Meth_select_var),]

library(SNFtool)
library(igraph)
library(Spectrum)
#library(kernlab)

##统一病人ID顺序
used_mut_with_mut <- Mut_info[order(rownames(Mut_info)),]
used_cnv_with_mut <- t(used_cnv)
used_cnv_with_mut <- used_cnv_with_mut[order(rownames(used_cnv_with_mut)),]
used_RNA_with_mut <- t(used_RNA)
used_RNA_with_mut <- used_RNA_with_mut[order(rownames(used_RNA_with_mut)),]
used_Meth_with_mut <- t(used_Meth)
used_Meth_with_mut <- used_Meth_with_mut[order(rownames(used_Meth_with_mut)),]

RNA_number <- 3000
DNA_number <- 250#50#100#200#300#400
cnv_number <- 520
brain1 <- t(used_mut_with_mut)
brain2 <- t(used_cnv_with_mut[,1:cnv_number])
brain3 <- t(used_RNA_with_mut[,1:RNA_number])#
brain4 <- t(used_Meth_with_mut[,1:DNA_number])#

Dist1 = dist2(t(as.matrix(brain1)),t(as.matrix(brain1)))
Dist2 = dist2(t(as.matrix(brain2)),t(as.matrix(brain2)))
Dist3 = dist2(t(as.matrix(brain3)),t(as.matrix(brain3)))
Dist4 = dist2(t(as.matrix(brain4)),t(as.matrix(brain4)))
### 自动决定K 
K = 12#12#10#15		# number of neighbors, usually (10~30)小于11不行
## 12～17
alpha = 0.5  	# hyperparameter, usually (0.3~0.8)
T1 = 15#10 	# Number of Iterations, usually (10~20)
## next, construct similarity graphs
W1 = affinityMatrix(Dist1, K, alpha)
W2 = affinityMatrix(Dist2, K, alpha)
W3 = affinityMatrix(Dist3, K, alpha)
W4 = affinityMatrix(Dist4, K, alpha)
## next, we fuse all the graphs
## then the overall matrix can be computed by similarity network fusion(SNF):
W_SNF = SNF(list(W1,W2,W3,W4), K, T1)#W2,
uu_SNF <- estimateNumberOfClustersGivenGraph(W_SNF,NUMC=2:6)
C_SNF = uu_SNF$`Eigen-gap best`#2


source('./SNFtool/R/internal.R')
#####
contro_SNF <- function(Wall,lab, K = 20, t = 20) 
{
  check_wall_names <- function(Wall) {
    name_match <- function(names_A, names_B) {
      return(identical(dimnames(names_A), dimnames(names_B)))
    }
    return(all(unlist(lapply(Wall, FUN = name_match, Wall[[1]]))))
  }
  wall.name.check <- check_wall_names(Wall)
  wall.names <- dimnames(Wall[[1]])
  if (!wall.name.check) {
    warning("Dim names not consistent across all matrices in Wall.\n            Returned matrix will have no dim names.")
  }
  LW <- length(Wall)
  normalize <- function(X) {
    row.sum.mdiag <- rowSums(X) - diag(X)
    row.sum.mdiag[row.sum.mdiag == 0] <- 1
    X <- X/(2 * (row.sum.mdiag))
    diag(X) <- 0.5
    return(X)
  }
  
  ###先正规化
  newW <- vector("list", LW)
  nextW <- vector("list", LW)
  for (i in 1:LW) {
    Wall[[i]] <- normalize(Wall[[i]])
    Wall[[i]] <- (Wall[[i]] + t(Wall[[i]]))/2
  }
  for (i in 1:LW) {
    newW[[i]] <- (.dominateset(Wall[[i]], K))
  }
  ###整合过程
  ###
  lab<- lab
  for (i in 1:t) {
    for (j in 1:LW) {
      sumWJ <- matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
      for (k in 1:LW) {
        if (k != j) {
          if( lab[k] != lab[j])
          {
          sumWJ <- sumWJ + Wall[[k]]
          }
        }
      }
      nextW[[j]] <- newW[[j]] %*% (sumWJ/(LW - 1)) %*% 
        t(newW[[j]])
    }
   ###再正规化再平均 
    for (j in 1:LW) {
      Wall[[j]] <- normalize(nextW[[j]])
      Wall[[j]] <- (Wall[[j]] + t(Wall[[j]]))/2
    }
  }
  W <- matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
  for (i in 1:LW) {
    W <- W + Wall[[i]]
  }
  W <- W/LW
  W <- normalize(W)
  W <- (W + t(W))/2
  if (wall.name.check) {
    dimnames(W) <- wall.names
  }
  return(W)
}
Wall <- list(W2,W3,W4)
K=K
t=alpha
lab <- c('D','D','R','E')
W_DSNF = contro_SNF(list(W1,W2,W3,W4),lab, K, T1)#W2,
uu_contr <- estimateNumberOfClustersGivenGraph(W_DSNF,NUMC=2:6)
C_contr = uu_contr$`Eigen-gap best`#2
group_DSNF = spectralClustering(W_DSNF,C_contr)	# the final subtypes information
displayClusters(W_DSNF, group_DSNF)
length(rownames(t(as.matrix(brain1)))[which(group_DSNF==1)])
length(intersect(rownames(t(as.matrix(brain1)))[which(group_DSNF==1)],
                 used_group[which(used_group[,3] %in% "YongG"),1]))
setdiff(rownames(t(as.matrix(brain1)))[which(group_DSNF==1)],
        used_group[which(used_group[,3] %in% "YongG"),1])

length(rownames(t(as.matrix(brain1)))[which(group_DSNF==2)])
length(intersect(rownames(t(as.matrix(brain1)))[which(group_DSNF==2)],
                 used_group[which(used_group[,3] %in% "OldG"),1]))

###转化成network
###有RNA 和甲基化的
###相似的流程，得到最后的W
###看它的邻居的label
data <- read_excel('Integrated_table.xlsx',sheet=1)
others <- data[which(!(data$with_DNA =="Y")),]
On_methy <- others[which(others$`RNA-seq`!= "NA"),]
RNA_methy <- others[which(!(others$`RNA-seq`!= "NA")),]
used <- data[which(data$with_DNA=="Y"),]
used_cons <- used
myNorm <- readRDS('DNA_methylation_45_samples_myNorm_0110.RDS')
myNorm_sel_RNA_meth <- myNorm[,which(colnames(myNorm)%in% as.character(RNA_methy$ID2...1))]
#dim()
used_myNorm_sel_RNA_meth <- myNorm_sel_RNA_meth[which(rownames(myNorm_sel_RNA_meth) %in% Meth_select_var),]
Only_RNA_meth <- matrix(0,nrow=dim(used_myNorm_sel_RNA_meth)[2],
                        ncol=dim(used_Meth)[2])
rownames(Only_RNA_meth) <- colnames(used_myNorm_sel_RNA_meth)
colnames(Only_RNA_meth) <- colnames(used_Meth)
for(i in 1:dim(Only_RNA_meth)[1])
{
  for(j in 1:dim(Only_RNA_meth)[2])
  {
    Only_RNA_meth[i,j] <-dist(rbind(used_myNorm_sel_RNA_meth[,i], used_Meth[,j])) 
  }
}
info_RNA_meythy <- data.frame(t(Only_RNA_meth))
lab_RNA_meythy <- rbind()
for(i in 1:dim(info_RNA_meythy)[1])
{
  k1 <- which(sel_pheno_info$ID== rownames(info_RNA_meythy)[i])
  lab_RNA_meythy <- rbind(lab_RNA_meythy,c(sel_pheno_info$ID[k1],sel_pheno_info$Group[k1]))
}
###只有甲基化
###计算distance
###聚类最近





