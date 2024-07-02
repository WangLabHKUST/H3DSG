rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
#clin_dir <- '/Users/shanghaixia/Desktop/glioma_prognosis_method_multi-omics/Multi_omics_SP_related/'
#library(readxl)
#clin <- read_excel(paste0(clin_dir,'Sp_45_subtype_ID_synology_230426.xlsx'),sheet=1)
info_master_Cpg <- readRDS('info_master_cpgs_DMP_info.RDS')
master_CpGs <- readRDS('Matser_cpgs_info_for_SP36.RDS')
meth_Group <- cbind(rownames(master_CpGs),master_CpGs[,c(which(colnames(master_CpGs)
                                                      %in% "Group"))])
colnames(meth_Group) <- c("ID","Group")
meth_Group <- as.data.frame(meth_Group)
meth_Group$Group <- paste0("Group",meth_Group$Group)
anno_info <- as.data.frame(meth_Group$Group)
colnames(anno_info) <- "Group"
rownames(anno_info) <- meth_Group$ID
anno_info$Group <- as.factor(anno_info$Group)
### clustering
library(pheatmap)
ann_colors = list(
  Group=c('Group1'='#f03b20','Group2'='#2b8cbe')#,
  #Gender = c(Female = 'sandybrown', Male = 'steelblue1'),
  #Loc = c(Upper = "lightblue4",Lower="lightsteelblue")
)

p1 <- pheatmap(t(as.matrix(master_CpGs[,which(colnames(master_CpGs)
                                              %in% rownames(info_master_Cpg))])),
               scale ="none",  cutree_col = 2,
               clustering_distance_cols = "correlation",
               color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
               clustering_method = "ward.D2",
               show_colnames = F,show_rownames = F,
               annotation_col=anno_info,
               annotation_colors = ann_colors,main="Cluster by master CpGs"
)
p1
pdf('A_unsupervised_clustering_bymaster_cpgs_SP36.pdf',width=6,height = 8)
p1
dev.off()

library(factoextra)
library(naivebayes)
meth.pca <- prcomp((as.matrix(master_CpGs[,which(colnames(master_CpGs)
                  %in% rownames(info_master_Cpg))])), scale = F)
res.ind <- get_pca_ind(meth.pca)
PCA_redut <- res.ind[["coord"]]
PCA_redut <- PCA_redut[,c(1:4)]

PCA_redut <- as.data.frame(PCA_redut)
PCA_redut$lab <- rep("0",length(PCA_redut$Dim.1))
for(i in 1:dim(PCA_redut)[1])
{
  k1 <- which(rownames(master_CpGs) %in% rownames(PCA_redut)[i])
  if(length(k1)>0)
  {
    PCA_redut$lab[i]<- master_CpGs$Group[k1]
  }
}
PCA_redut$lab <- as.factor(PCA_redut$lab)
model_fina_PCA <- naive_bayes(lab~., usekernel = T,kernel = "biweight",
                              data=PCA_redut[,c(1,5)])
model_fina_PCA

pdf('Master_cpG_PCA_dim1_model_learn_SP36.pdf',width=5,height =5)
plot(model_fina_PCA,
     prob="conditional",which=1,ask=F,
     arg.num = list(col = c("#FF3300","#2b8cbe"), lty = c(1,1),lwd=1.5,
                    main = "PCA of CpGs", legend.position = "topright",
                    legend.cex = 0.5))
dev.off()


feature_loc_age <- which(colnames(master_CpGs) %in%c("Age","Group"))
model_finaGauss_age <- naive_bayes(Group~., 
                    data=master_CpGs[,c(feature_loc_age)])
pdf('Age_Gauss_model_learn_SP36.pdf',width=5,height =5)
plot(model_finaGauss_age,
     prob="conditional",which=1,ask=F,
     arg.num = list(col = c("#FF3300","#2b8cbe"), lty = c(1,1),lwd=1.5,
                    main = "Age", legend.position = "topright",
                    legend.cex = 0.5))
dev.off()

feature_loc_KI67 <- which(colnames(master_CpGs) %in%c("KI67","Group"))
master_CpGs$KI67 <- as.numeric(master_CpGs$KI67)
model_fina_Gauss_KI67 <- naive_bayes(Group~., 
                    data=master_CpGs[,c(feature_loc_KI67)],
                    usekernel =F)
pdf('KI67_Gauss_model_learn_SP36.pdf',width=5,height =5)
plot(model_fina_Gauss_KI67,
     prob="conditional",which=1,ask=F,
     arg.num = list(col = c("#FF3300","#2b8cbe"), lty = c(1,1),lwd=1.5,
                    main = "KI67", legend.position = "topright",
                    legend.cex = 0.5))
dev.off()

feature_loc_TP53 <- which(colnames(master_CpGs) %in%c("TP53","Group"))
master_CpGs$TP53 <- as.character(master_CpGs$TP53)
model_fina_TP53 <- naive_bayes(Group~., 
                               data=master_CpGs[,c(feature_loc_TP53)])

pdf('TP53_model_learn_SP36.pdf',width=5,height =5)
plot(model_fina_TP53,
     prob="conditional",which=1,ask=F,
     arg.num = list(col = c("#FF3300","#2b8cbe"), lty = c(1,1),lwd=1.5,
                    main = "TP53", legend.position = "topright",
                    legend.cex = 0.5))
dev.off()

feature_loc_NF1 <- which(colnames(master_CpGs) %in%c("NF1","Group"))
master_CpGs$NF1 <- as.character(master_CpGs$NF1)
model_fina_NF1 <- naive_bayes(Group~., 
                               data=master_CpGs[,c(feature_loc_NF1)])
pdf('NF1_model_learn_SP36.pdf',width=5,height =5)
plot(model_fina_NF1,
     prob="conditional",which=1,ask=F,
     arg.num = list(col = c("#FF3300","#2b8cbe"), lty = c(1,1),lwd=1.5,
                    main = "NF1", legend.position = "topright",
                    legend.cex = 0.5))
dev.off()



