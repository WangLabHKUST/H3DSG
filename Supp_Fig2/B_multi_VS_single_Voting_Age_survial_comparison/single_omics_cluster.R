rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
info <- readxl::read_excel("single_omics_data_clin.xlsx",sheet=1)

SP_36_meth_selected <- t(readRDS('./Four_omics_data/meth_info.RDS'))
dim(SP_36_meth_selected)
colnames(SP_36_meth_selected)
renew_col_meth <- rbind()
for(i in 1:dim(SP_36_meth_selected)[2])
{
  k1 <- which(info$`Sample ID` %in% colnames(SP_36_meth_selected)[i])
  renew_col_meth <- rbind(renew_col_meth,
                          info$M_ID[k1])
}
colnames(SP_36_meth_selected ) <- renew_col_meth
anno_info <- as.data.frame(info[,c(4,5,6)])
rownames(anno_info) <- info$M_ID
colnames(anno_info)[c(2,3)] <- c("Gender","Loc")
#anno_info$Group <- as.factor(anno_info$Group)
anno_info$Gender <- as.factor(anno_info$Gender)
anno_info$Age <- as.numeric(anno_info$Age)
anno_info$Loc <- as.factor(anno_info$Loc)
library(pheatmap)
ann_colors = list(
  #Group=c(Group1='#f03b20',Group2='#2b8cbe'),
  Gender = c(Female = 'sandybrown', Male = 'steelblue1'),
  Loc = c(Upper = "lightblue4",Lower="lightsteelblue")
)

pdf('SFig1_Methy_cluster_add_annotation.pdf',width=8,height = 12)
pheatmap(as.matrix(SP_36_meth_selected),
         scale ="none",  cutree_col = 2,
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         clustering_method = "ward.D2",
         show_colnames = F,show_rownames = F,
         annotation_col=anno_info,
         annotation_colors = ann_colors,main="Cluster by methylation"
         )
dev.off()
meth_group <- pheatmap(as.matrix(SP_36_meth_selected),
                       scale ="none",  cutree_col = 2,
                       clustering_distance_cols = "correlation",
                       color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                       clustering_method = "ward.D2")
meth_col_cluster <- as.data.frame(cutree(meth_group$tree_col,k=2))
colnames(meth_col_cluster) <- "meth"

SP_36_cnv_selected_together <- t(readRDS('./Four_omics_data/cnv_info.RDS'))
renew_col_cnv_t <- rbind()
for(i in 1:dim(SP_36_cnv_selected_together)[2])
{
  k1 <- which(info$`Sample ID` %in% colnames(SP_36_cnv_selected_together)[i])
  renew_col_cnv_t <- rbind(renew_col_cnv_t,
                          info$M_ID[k1])
}
colnames(SP_36_cnv_selected_together) <- renew_col_cnv_t

pdf('SFig1_CNV_cluster_add_annotation.pdf',
    width=8,height = 10)
pheatmap(as.matrix(SP_36_cnv_selected_together),
         scale ="row",  cutree_col = 2,
         clustering_distance_cols = "correlation",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         clustering_method = "ward.D2",
         show_colnames = F,show_rownames = F,annotation_col=anno_info,
         annotation_colors = ann_colors,main="Cluster by CNV"
         )
dev.off()
cnv_group <- pheatmap(as.matrix(SP_36_cnv_selected_together),
                      scale ="row",  cutree_col = 2,
                      clustering_distance_cols = "correlation",
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                      clustering_method = "ward.D2",
                      show_colnames = F,show_rownames = F)
cnv_col_cluster <- as.data.frame(cutree(cnv_group$tree_col,k=2))
colnames(cnv_col_cluster) <- "cnv"

SP_36_RNA_seq_selected <- t(readRDS('./Four_omics_data/rna_info.RDS'))
renew_col_rna <- rbind()
for(i in 1:dim(SP_36_RNA_seq_selected)[2])
{
  k1 <- which(info$`Sample ID` %in% colnames(SP_36_RNA_seq_selected)[i])
  renew_col_rna  <- rbind(renew_col_rna ,
                           info$M_ID[k1])
}
colnames(SP_36_RNA_seq_selected) <- renew_col_rna 


pdf('SFig1_RNA_cluster_add_annotation.pdf',
    width=8,height = 10)
pheatmap(as.matrix(SP_36_RNA_seq_selected),
         cutree_col = 2,
         clustering_method = 'ward.D2',
         color = colorRampPalette(c("navy", "white", "firebrick3"))(15),
         scale="row",
         clustering_distance_cols ="correlation",
         show_rownames = F,show_colnames =F,annotation_col=anno_info,
         annotation_colors = ann_colors,main="Cluster by RNA"
         )
dev.off()
rna_group <- pheatmap(as.matrix(SP_36_RNA_seq_selected),
                      cutree_col = 2,
                      clustering_method = 'ward.D2',
                      color = colorRampPalette(c("navy", "white", "firebrick3"))(15),
                      scale="row",
                      clustering_distance_cols ="correlation",
                      show_rownames = F,show_colnames =F)
rna_col_cluster <- as.data.frame(cutree(rna_group$tree_col,k=2))
colnames(rna_col_cluster)<- "rna"

Mutat_matrix <- t(readRDS('./Four_omics_data/driver_SNV.RDS'))
#dim(Mutat_matrix)
renew_col_som <- rbind()
for(i in 1:dim(Mutat_matrix)[2])
{
  k1 <- which(info$`Sample ID` %in% colnames(Mutat_matrix)[i])
  renew_col_som  <- rbind(renew_col_som,
                          info$M_ID[k1])
}
colnames(Mutat_matrix) <- renew_col_som

pdf('SFig1_Mutation_cluster_add_annotation.pdf',
    width=8,height = 6)
pheatmap(as.matrix(Mutat_matrix),
         scale ="none",  cutree_col = 2,
         clustering_distance_cols = "correlation",
         color = c('white','black'),
         clustering_method = "ward.D2",
         show_colnames = F,show_rownames = T,
         annotation_col=anno_info,
         annotation_colors = ann_colors,main="Cluster by Somatic"
         )
dev.off()

mut_group <- pheatmap(as.matrix(Mutat_matrix),
                      scale ="none",  cutree_col = 2,
                      clustering_distance_cols = "correlation",
                      color = c('white','black'),
                      clustering_method = "ward.D2",
                      show_colnames = F,show_rownames = T)
mut_col_cluster <- as.data.frame(cutree(mut_group$tree_col,k=2))
colnames(mut_col_cluster) <- c("mut")
####generate the survival and age comparison
group <- merge(meth_col_cluster,cnv_col_cluster,by="row.names",all=T)
rownames(group)<- group$Row.names
group <- group[,-1]
group <- merge(group,rna_col_cluster,by="row.names",all=T)
rownames(group)<- group$Row.names
group <- group[,-1]
group <- merge(group,mut_col_cluster,by="row.names",all=T)
rownames(group)<- group$Row.names
group <- group[,-1]
used_info <- as.data.frame(info[,c(2:6,15:16)])
rownames(used_info) <- used_info$M_ID
used_info <- merge(used_info,group,by="row.names",all=T)
rownames(used_info) <- used_info$Row.names
used_info <- used_info[,-1]
###
median(used_info$Age[which(used_info$Group_multi_omics %in% 
                             "Group1")])
min(used_info$Age[which(used_info$Group_multi_omics %in% 
                          "Group1")])
max(used_info$Age[which(used_info$Group_multi_omics %in% 
                          "Group1")])
median(used_info$Age[which(used_info$Group_multi_omics %in% 
                             "Group2")])
min(used_info$Age[which(used_info$Group_multi_omics %in% 
                          "Group2")])
max(used_info$Age[which(used_info$Group_multi_omics %in% 
                          "Group2")])
wilcox.test(used_info$Age[which(used_info$Group_multi_omics %in% "Group1")],
          used_info$Age[which(used_info$Group_multi_omics %in% "Group2")],
          paired = F)
library("survival")
library("survminer")
multi_fit <- survfit(Surv(OS, sensor) ~ Group_multi_omics, data = used_info)
print(multi_fit)
ggsurvplot(multi_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#mth_clu_info <- used_info[,c(1,3,6:8)]#used_info[,c(1,3,6:7,9)]#
#x <- mth_clu_info
re_sta_index <- function(x)
{
  colnames(x)[5]<- 'sub'
  age_1 <- x$Age[which(x$sub==1)]
  age_2 <- x$Age[which(x$sub==2)]
  age_test <- wilcox.test(age_1,age_2 ,paired = F)
  age_info <- c(median(age_1),min(age_1),max(age_1),
                median(age_2),min(age_2),max(age_2),
                age_test$p.value
  )
  return(age_info)
}

table(used_info$meth)
table(used_info$cnv)
table(used_info$rna)
table(used_info$mut)
meth_age_sta <- re_sta_index(used_info[,c(1,3,6:8)])
meth_age_sta
meth_fit <- survfit(Surv(OS, sensor) ~ meth, data = used_info)
print(meth_fit)
ggsurvplot(meth_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

cnv_age_sta <- re_sta_index(used_info[,c(1,3,6:7,9)])
cnv_age_sta
cnv_fit <- survfit(Surv(OS, sensor) ~ cnv, data = used_info)
print(cnv_fit)
ggsurvplot(cnv_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


rna_age_sta <- re_sta_index(used_info[,c(1,3,6:7,10)])
rna_age_sta
rna_fit <- survfit(Surv(OS, sensor) ~ rna, data = used_info)
print(rna_fit)
ggsurvplot(rna_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

som_age_sta <- re_sta_index(used_info[,c(1,3,6:7,11)])
som_age_sta
som_fit <- survfit(Surv(OS, sensor) ~ mut, data = used_info)
print(som_fit)
ggsurvplot(som_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#write.table(used_info,'Clustering_by_single_omics.txt',
#            sep="\t",row.names = F,quote=F)

voting_SCM <- readxl::read_excel('Clustering_by_single_omics.xlsx',
                                 sheet=2)
table(voting_SCM$Voting)
re_sta_index_vo <- function(x)
{
  colnames(x)[5]<- 'sub'
  age_1 <- x$Age[which(x$sub%in% "Group1")]
  age_2 <- x$Age[which(x$sub%in% "Group2")]
  age_test <- wilcox.test(age_1,age_2 ,paired = F)
  age_info <- c(median(age_1),min(age_1),max(age_1),
                median(age_2),min(age_2),max(age_2),
                age_test$p.value
  )
  return(age_info)
}
voting_SCM_age_sta <- re_sta_index_vo(voting_SCM[,c(1,3,4:5,9)])
voting_SCM_age_sta
voting_SCM_fit <- survfit(Surv(OS, sensor) ~ Voting, data = voting_SCM)
print(voting_SCM_fit)
ggsurvplot(voting_SCM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
voting_SCR <- readxl::read_excel('Clustering_by_single_omics.xlsx',
                                 sheet=3)
table(voting_SCR$Voting)
voting_SCR_age_sta <- re_sta_index_vo(voting_SCR[,c(1,3,4:5,9)])
voting_SCR_age_sta
voting_SCR_fit <- survfit(Surv(OS, sensor) ~ Voting, data = voting_SCR)
print(voting_SCR_fit)
ggsurvplot(voting_SCR_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


voting_SMR <- readxl::read_excel('Clustering_by_single_omics.xlsx',
                                 sheet=4)
table(voting_SMR$Voting)
voting_SMR_age_sta <- re_sta_index_vo(voting_SMR[,c(1,3,4:5,9)])
voting_SMR_age_sta
voting_SMR_fit <- survfit(Surv(OS, sensor) ~ Voting, data = voting_SMR)
print(voting_SMR_fit)
ggsurvplot(voting_SMR_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
voting_CMR <- readxl::read_excel('Clustering_by_single_omics.xlsx',
                                 sheet=5)
table(voting_CMR$Voting)
voting_CMR_age_sta <- re_sta_index_vo(voting_CMR[,c(1,3,4:5,9)])
voting_CMR_age_sta
voting_CMR_fit <- survfit(Surv(OS, sensor) ~ Voting, data = voting_CMR)
print(voting_CMR_fit)
ggsurvplot(voting_CMR_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))