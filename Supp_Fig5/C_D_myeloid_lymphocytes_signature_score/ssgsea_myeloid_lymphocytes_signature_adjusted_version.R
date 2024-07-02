#path="/Users/shanghaixia/Desktop/glioma_prognosis_method_multi-omics/Multi_omics_SP_related/Code_for_github/Fig2_volcano_plot" 
#setwd(path)
rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
col_order<- c("SPH3_37","SPH3_24","SPH3_14","SPH3_20",
          "SPH3_29","SPH3_07","SPH3_25","SPH3_12",
          "SPH3_40","SPH3_36","SPH3_06","SPH3_32",
          "SPH3_33","SPH3_44","SPH3_10","SPH3_21",
          "SPH3_13","SPH3_11","SPH3_16","SPH3_03",
          "SPH3_17","SPH3_18","SPH3_30","SPH3_09",
          "SPH3_38","SPH3_43","SPH3_45","SPH3_08",
          "SPH3_15","SPH3_23","SPH3_04","SPH3_19",
          "SPH3_27","SPH3_31","SPH3_05","SPH3_34")
mark_set <- vector('list',length=2)
## Lymphocytes signature genes
dat1_lym <- readxl::read_excel('cell_marker_lymphocyte.xlsx',sheet=1)
mark_set[[1]] <- as.matrix(dat1_lym$marker)
  #c("CD3","CD4","CD8","CD19","CD20","CD56","CD16")

## Myeloid signature genes
dat1_myel <- readxl::read_excel("myeloid_cell_marker_brain.xlsx",sheet=1)
mark_set[[2]] <- as.matrix(dat1_myel$marker)#c("CD15","CD16","CD14","CD68","CD11c","CD123",
                 #  "CD9","CD69","CD123","CD203c")
gene_exp <- readRDS('Raw_RAN_step1_tpm_log2_quan_for_36_sample1219.RDS')
group_info <- readRDS("label_for_twoGroups_SP36_1219.RDS")
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)    #need R 3.6
library(pheatmap)
gsva_matrix<- gsva(as.matrix(gene_exp), mark_set,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
new_name <- rbind()
for(i in 1:dim(gsva_matrix)[2])
{
  k1 <- which(group_info$ID %in% gsub("P","",colnames(gsva_matrix)[i]))
  new_name <- rbind(new_name,
                    group_info$name[k1])
}
colnames(gsva_matrix) <- new_name
gsav_matrix_with_order <- cbind()
for(i in 1:length(col_order))
{
  k1 <- which(colnames(gsva_matrix) %in% col_order[i])
  gsav_matrix_with_order <- cbind(gsav_matrix_with_order,
                                  gsva_matrix[,k1])
}
rownames(gsav_matrix_with_order) <- rownames(gsva_matrix)
colnames(gsav_matrix_with_order) <- col_order
gsav_matrix_with_order <- as.data.frame(t(gsav_matrix_with_order))
colnames(gsav_matrix_with_order) <- c("Lymphocytes_signature",
                                      "Myeloid_signature")
gsav_matrix_with_order$Group <- c(rep("SP_group1",16),
                                  rep("SP_group2",20))
gsav_matrix_with_order$Group <- as.factor(gsav_matrix_with_order$Group)
library(ggplot2)
library(ggpubr)
library(ggrepel)

library(ggbeeswarm)
NAN_plot <- ggplot(data=gsav_matrix_with_order,
                   aes(x=Group,y=Lymphocytes_signature)) + theme_classic() +
  geom_boxplot(data=gsav_matrix_with_order,
               aes(x=Group,y=Lymphocytes_signature),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=gsav_matrix_with_order,aes(x=Group,y=Lymphocytes_signature, color=Group, 
                                     fill=Group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  xlab("")+ylim(0,1.5)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",Lymphocytes_signature ~ Group,  
                   data = gsav_matrix_with_order,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.4)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Lymphocytes_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=gsav_matrix_with_order,
                   aes(x=Group,y=Myeloid_signature)) + theme_classic() +
  geom_boxplot(data=gsav_matrix_with_order,
               aes(x=Group,y=Myeloid_signature),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=gsav_matrix_with_order,aes(x=Group,y=Myeloid_signature, color=Group, 
                                                   fill=Group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  xlab("")+ylim(0.6,1.3)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",Myeloid_signature ~ Group,  
                   data = gsav_matrix_with_order,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.22)
)
figure_21<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Myeloid_signature_box_comparison.pdf",
       plot=figure_21,bg = 'white', width =6, height = 8, 
       units = 'cm', dpi = 600)






if(F){
plt_GBP2 <- df_GBP2_3C

plt_GBP2$group <- ifelse(plt_GBP2$MES_P == "Yes", "Mut", "WT")
NAN_plot <- ggplot(data=plt_GBP2,aes(x=group,y=GBP2_scale)) + theme_classic() 
NAN_plot <- NAN_plot + 
  geom_boxplot(data=plt_GBP2,aes(x=group,y=GBP2_scale),width = 0.4,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=plt_GBP2,aes(x=group,y=GBP2_scale, color=group, fill=group),width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)
NAN_plot <-NAN_plot + theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
                            plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
                            legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
                            axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
                            axis.title.x=element_text(size=14,face='plain',color='black'),axis.title.y=element_text(size=14,face='plain',color='black'))
NAN_plot
NAN_plot <- NAN_plot + scale_y_continuous(expand=c(0,0),limits=c(-2,4),breaks = seq(-10,20,1)) 
NAN_plot<- NAN_plot +ylab("GBP2 Relative Abundance") +xlab(NULL)
NAN_plot <- NAN_plot + scale_fill_manual(name=NULL,values = c(Mut="#ff7f00",WT="#984ea3"))
NAN_plot <- NAN_plot + scale_color_manual(name=NULL,values =  c(Mut="#ff7f00",WT="#984ea3"))
library(ggpubr)
the<-compare_means(method = "wilcox.test",GBP2_scale ~ group,  data = plt_GBP2,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(4)
)
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="./figures/Fig6/1213_2D_GBP2_3C_validatio.pdf", plot=figure_2,bg = 'white', width =9.5, height = 9, units = 'cm', dpi = 600)
}

