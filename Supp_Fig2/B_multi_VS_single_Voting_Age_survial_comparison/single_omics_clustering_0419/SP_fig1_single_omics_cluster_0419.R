rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
info <- readxl::read_excel("single_omics_data_integrated_0821.xlsx",sheet=1)
info$final_dec[which(info$final_dec %in% "YongG")] <- "Group1"
info$final_dec[which(info$final_dec %in% "OldG")] <- "Group2"

SP_36_meth_selected <- readRDS('DNA_methy_top20000_probes.RDS')
#write.table(rownames(SP_36_meth_selected),
#            "DNA_methy_top20000_probes_name.txt",sep="\t",
#            row.names = F,quote = F)
dim(SP_36_meth_selected)
colnames(SP_36_meth_selected)
renew_col_meth <- rbind()
for(i in 1:dim(SP_36_meth_selected)[2])
{
  k1 <- which(info$ID2 %in% colnames(SP_36_meth_selected)[i])
  renew_col_meth <- rbind(renew_col_meth,
                          info$M_ID[k1])
}
colnames(SP_36_meth_selected ) <- renew_col_meth
anno_info <- as.data.frame(info[,c(6,7,8)])#4,
rownames(anno_info) <- info$M_ID
colnames(anno_info)[c(2,3)] <- c("Gender","Loc")#1,"Group",
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

pdf('SFig1_Methy_cluster_add_annotation_0419.pdf',width=8,height = 12)
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


SP_36_cnv_selected_together <- readRDS('cnv_top550_cytobands_together.RDS')
renew_col_cnv_t <- rbind()
for(i in 1:dim(SP_36_cnv_selected_together)[2])
{
  k1 <- which(info$ID2 %in% colnames(SP_36_cnv_selected_together)[i])
  renew_col_cnv_t <- rbind(renew_col_cnv_t,
                          info$M_ID[k1])
}
colnames(SP_36_cnv_selected_together) <- renew_col_cnv_t

pdf('SFig1_CNV_cluster_add_annotation_0419.pdf',
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

SP_36_RNA_seq_selected <- readRDS('RNA_log2_tpm_1_top3000.RDS')
renew_col_rna <- rbind()
for(i in 1:dim(SP_36_RNA_seq_selected)[2])
{
  k1 <- which(info$ID2 %in% colnames(SP_36_RNA_seq_selected)[i])
  renew_col_rna  <- rbind(renew_col_rna ,
                           info$M_ID[k1])
}
colnames(SP_36_RNA_seq_selected) <- renew_col_rna 


pdf('SFig1_RNA_cluster_add_annotation_0419.pdf',
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

Mutat_matrix <- readRDS('mutation_for_clustering.RDS')
renew_col_som <- rbind()
for(i in 1:dim(Mutat_matrix)[2])
{
  k1 <- which(info$ID2 %in% colnames(Mutat_matrix)[i])
  renew_col_som  <- rbind(renew_col_som,
                          info$M_ID[k1])
}
colnames(Mutat_matrix) <- renew_col_som

pdf('SFig1_Mutation_cluster_add_annotation_0419.pdf',
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