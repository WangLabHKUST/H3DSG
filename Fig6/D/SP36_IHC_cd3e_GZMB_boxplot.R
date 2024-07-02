rm(list=ls())
library(rstudioapi)
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
IHC_info_36 <- readRDS('IHC_info_36.RDS')
used_exp_CD33_36 <- as.data.frame(IHC_info_36[,c(1,2,3,6)])
used_exp_CD33_36$CD3  <- as.numeric(used_exp_CD33_36$CD3...16)
used_exp_CD33_36  <- used_exp_CD33_36[which(is.na(used_exp_CD33_36$CD3)==F),]

NAN_plot <- ggplot(data=used_exp_CD33_36,
                   aes(x=final_dec,y=CD3)) + theme_classic() +
  geom_boxplot(data=used_exp_CD33_36,
               aes(x=final_dec,y=CD3),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_CD33_36,
                   aes(x=final_dec,y=CD3, color=final_dec, 
                       fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#FF3300',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#FF3300',Group2='#2b8cbe'))+
  xlab("")+ylim(0,1.7*100)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",CD3 ~ final_dec,  
                   data = used_exp_CD33_36,paired = F )

NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.65*100)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig6_CD3E_36_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)

colnames(IHC_info_36)[which(colnames(IHC_info_36)%in% "GZMA")] <- "GZMB"
used_exp_GZMB_36 <- as.data.frame(IHC_info_36[,c(1,2,3,12)])
used_exp_GZMB_36$GZMB  <- as.numeric(used_exp_GZMB_36$GZMB)
used_exp_GZMB_36  <- used_exp_GZMB_36[which(is.na(used_exp_GZMB_36$GZMB)==F),]

NAN_plot <- ggplot(data=used_exp_GZMB_36,
                   aes(x=final_dec,y=GZMB)) + theme_classic() +
  geom_boxplot(data=used_exp_GZMB_36,
               aes(x=final_dec,y=GZMB),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_GZMB_36,
                   aes(x=final_dec,y=GZMB, color=final_dec, 
                       fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#FF3300',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#FF3300',Group2='#2b8cbe'))+
  xlab("")+ylim(0,22*10)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",GZMB ~ final_dec,  
                   data = used_exp_GZMB_36,paired = F )

NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(19*10)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig6_GZMB_36_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)

