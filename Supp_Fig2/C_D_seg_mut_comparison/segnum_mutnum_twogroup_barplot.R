rm(list = ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
seg_TMB <- readxl::read_excel('seg_TMB_cm_info.xlsx',sheet = 1)
seg_TMB$log10_transMutcount
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
NAN_plot <- ggplot(data=seg_TMB,
            aes(x=final_dec,y=segments)) + theme_classic() +
  geom_boxplot(data=seg_TMB,
               aes(x=final_dec,y=segments),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=seg_TMB,aes(x=final_dec,y=segments, 
                                        color=final_dec, 
                                                   fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  xlab("")+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",segments ~ final_dec,  
                   data = seg_TMB,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(8)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="segments_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, 
       units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=seg_TMB,
                   aes(x=final_dec,y=log10_transMutcount)) + theme_classic() +
  geom_boxplot(data=seg_TMB,
               aes(x=final_dec,y=log10_transMutcount),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=seg_TMB,aes(x=final_dec,y=log10_transMutcount, 
                                    color=final_dec, 
                                    fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(SP_group1='#f03b20',SP_group2='#2b8cbe'))+
  xlab("")+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",log10_transMutcount ~ final_dec,  
                   data = seg_TMB,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(2.5)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="TMB_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, 
       units = 'cm', dpi = 600)


