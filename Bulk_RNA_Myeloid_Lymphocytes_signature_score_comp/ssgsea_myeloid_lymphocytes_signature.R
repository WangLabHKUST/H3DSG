#path="/Users/shanghaixia/Desktop/glioma_prognosis_method_multi-omics/Multi_omics_SP_related/Code_for_github/Fig2_volcano_plot" 
#setwd(path)
rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
gsav_matrix_with_order <- readRDS('gsav_matrix_with_order.RDS')

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








