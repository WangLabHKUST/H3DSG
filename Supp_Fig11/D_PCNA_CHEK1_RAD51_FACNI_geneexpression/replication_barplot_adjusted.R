rm(list = ls())
library(rstudioapi)
library(readxl)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())

used_exp <- readRDS('DNA_replication_gene.RDS')
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=PCNA)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=PCNA),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=PCNA, color=group, 
                                        fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#990000',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#990000',Group2='#2b8cbe'))+
  xlab("")+ylim(0,12)+
  theme(axis.text.x=element_blank())+##f03b20
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",PCNA ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(11)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="PCNA_RNA_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, 
       height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=RAD51)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=RAD51),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=RAD51, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#990000',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#990000',Group2='#2b8cbe'))+
  xlab("")+ylim(0,7)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",RAD51 ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(7)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="RAD51_RNA_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=FANCI)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=FANCI),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=FANCI, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#990000',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#990000',Group2='#2b8cbe'))+
  xlab("")+ylim(0,7)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",FANCI ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(7)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="FANCI_RNA_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=CHEK1)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=CHEK1),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=CHEK1, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#990000',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#990000',Group2='#2b8cbe'))+
  xlab("")+ylim(0,7)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",CHEK1 ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(7)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="CHEK1_RNA_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

