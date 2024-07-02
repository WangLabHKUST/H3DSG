rm(list = ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
if (!require('ggplot2')) 
  install.packages('ggplot2')
if (!require('ggpubr')) 
  install.packages('ggpubr')
if (!require('ggrepel')) 
  install.packages('ggrepel')
if (!require('ggbeeswarm')) 
  install.packages('ggbeeswarm')

library(rstudioapi)
library(readxl)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
used_exp <- readRDS("E2F1_MYC_MEX3A_gene_exp_tpm_log2.RDS")
unique(used_exp$group)
used_exp$group <- as.character(used_exp$group)
used_exp$group[which(used_exp$group %in% "YongG")] <- "Group1"
used_exp$group[which(used_exp$group %in% "OldG")] <- "Group2"
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=E2F1)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=E2F1),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=E2F1, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,9)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",E2F1 ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(8.8)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_E2F1_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=MYC)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=MYC),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=MYC, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,11)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",MYC ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(10.8)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_MYC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

NAN_plot <- ggplot(data=used_exp,
                   aes(x=group,y=MEX3A)) + theme_classic() +
  geom_boxplot(data=used_exp,
               aes(x=group,y=MEX3A),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp,aes(x=group,y=MEX3A, color=group, 
                                     fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,8)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",MEX3A ~ group,  
                   data = used_exp,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(7.8)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_MEX3A_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)

#### IHC
used_exp_IHC <- readRDS("E2F1_MYC_MEX3A_IHC.RDS")
colnames(used_exp_IHC)[3] <- "MYC"
used_exp_IHC <- used_exp_IHC[which(is.na(used_exp_IHC$E2F1)==F),]
NAN_plot <- ggplot(data=used_exp_IHC,
                   aes(x=final_dec,y=E2F1)) + theme_classic() +
  geom_boxplot(data=used_exp_IHC,
               aes(x=final_dec,y=E2F1),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_IHC,aes(x=final_dec,y=E2F1, color=final_dec, 
                  fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,110)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",E2F1 ~ final_dec,  
                   data = used_exp_IHC,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(108)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_E2F1_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)


used_exp_IHC <- readRDS("E2F1_MYC_MEX3A_IHC.RDS")
colnames(used_exp_IHC)[3] <- "MYC"
used_exp_IHC <- used_exp_IHC[which(is.na(used_exp_IHC$MYC)==F),]

NAN_plot <- ggplot(data=used_exp_IHC,
                   aes(x=final_dec,y=MYC)) + theme_classic() +
  geom_boxplot(data=used_exp_IHC,
               aes(x=final_dec,y=MYC),
               width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_IHC,aes(x=final_dec,y=MYC, 
                color=final_dec, fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,100)+ylab("c-MYC")+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",MYC ~ final_dec,  
                   data = used_exp_IHC,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(96)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_MYC_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)


used_exp_IHC <- readRDS("E2F1_MYC_MEX3A_IHC.RDS")
colnames(used_exp_IHC)[3] <- "MYC"
used_exp_IHC <- used_exp_IHC[which(is.na(used_exp_IHC$MEX3A)==F),]
NAN_plot <- ggplot(data=used_exp_IHC,
                   aes(x=final_dec,y=MEX3A)) + theme_classic() +
  geom_boxplot(data=used_exp_IHC,
               aes(x=final_dec,y=MEX3A),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_IHC,aes(x=final_dec,y=MEX3A, color=final_dec, 
                                          fill=final_dec),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,140)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",MEX3A ~ final_dec,  
                   data = used_exp_IHC,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(134)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig5_MEX3A_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height = 8, units = 'cm', dpi = 600)


