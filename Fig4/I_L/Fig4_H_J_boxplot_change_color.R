rm(list = ls())
library(rstudioapi)
library(readxl)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)

used_exp_IHC<- readRDS("PCDNA_RAD51_PKCs_H2AX_FANCI_IHC.RDS")

used_exp_FANCI <- as.data.frame(used_exp_IHC[,c(1,2,7)])
used_exp_FANCI$FANCI  <- as.numeric(used_exp_FANCI$FANCI)
used_exp_FANCI  <- used_exp_FANCI[which(is.na(used_exp_FANCI$FANCI)==F),]

NAN_plot <- ggplot(data=used_exp_FANCI ,
                   aes(x=final_dec,y=FANCI)) + theme_classic() +
  geom_boxplot(data=used_exp_FANCI ,
               aes(x=final_dec,y=FANCI),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_FANCI ,aes(x=final_dec,y=FANCI, color=final_dec, 
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
  xlab("")+ylim(0,18*100)+
  theme(axis.text.y=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",FANCI ~ final_dec,  
                   data = used_exp_FANCI,paired = F )

NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(17.5*100)
)
NAN_plot<- NAN_plot+coord_flip()
NAN_plot
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig6_FANCI_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =8, height = 6, units = 'cm', dpi = 600)

used_exp_PCNA <- as.data.frame(used_exp_IHC[,c(1,2,5)])
used_exp_PCNA$PCNA  <- as.numeric(used_exp_PCNA$pCNA)
used_exp_PCNA  <- used_exp_PCNA[which(is.na(used_exp_PCNA$PCNA)==F),]

NAN_plot <- ggplot(data=used_exp_PCNA,
                   aes(x=final_dec,y=PCNA)) + theme_classic() +
  geom_boxplot(data=used_exp_PCNA,
               aes(x=final_dec,y=PCNA),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_PCNA,aes(x=final_dec,y=PCNA, color=final_dec, 
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
  xlab("")+ylim(0,1.5*100)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",PCNA ~ final_dec,  
                   data = used_exp_PCNA,paired = F )
NAN_plot
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.45*100)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig4_PCNA_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)

used_exp_RAD51 <- as.data.frame(used_exp_IHC[,c(1,2,6)])
used_exp_RAD51$RAD51  <- as.numeric(used_exp_RAD51$RAD51)
used_exp_RAD51  <- used_exp_RAD51[which(is.na(used_exp_RAD51$RAD51)==F),]

NAN_plot <- ggplot(data=used_exp_RAD51,
                   aes(x=final_dec,y=RAD51)) + theme_classic() +
  geom_boxplot(data=used_exp_RAD51,
               aes(x=final_dec,y=RAD51),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_RAD51,aes(x=final_dec,y=RAD51, color=final_dec, 
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
  xlab("")+ylim(0,1.3*100)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",RAD51 ~ final_dec,  
                   data = used_exp_RAD51,paired = F )
NAN_plot
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.25*100)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig4_RAD51_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)

used_exp_DNApk <- as.data.frame(used_exp_IHC[,c(1,2,4)])
used_exp_DNApk$DNApk  <- as.numeric(used_exp_DNApk$DNApk)
used_exp_DNApk  <- used_exp_DNApk[which(is.na(used_exp_DNApk$DNApk)==F),]

NAN_plot <- ggplot(data=used_exp_DNApk,
                   aes(x=final_dec,y=DNApk)) + theme_classic() +
  geom_boxplot(data=used_exp_DNApk,
               aes(x=final_dec,y=DNApk),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_DNApk,aes(x=final_dec,y=DNApk, color=final_dec, 
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
  xlab("")+ylim(0,1.5*100)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",DNApk ~ final_dec,  
                   data = used_exp_DNApk,paired = F )
NAN_plot
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.45*100)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig4_DNApk_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)


used_exp_rH2AX <- as.data.frame(used_exp_IHC[,c(1,2,3)])
used_exp_rH2AX$rH2AX  <- as.numeric(used_exp_rH2AX$rH2AX)
used_exp_rH2AX  <- used_exp_rH2AX[which(is.na(used_exp_rH2AX$rH2AX)==F),]

NAN_plot <- ggplot(data=used_exp_rH2AX,
                   aes(x=final_dec,y=rH2AX)) + theme_classic() +
  geom_boxplot(data=used_exp_rH2AX,
               aes(x=final_dec,y=rH2AX),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+#,color=final_dec
  geom_quasirandom(data=used_exp_rH2AX,aes(x=final_dec,y=rH2AX, color=final_dec, 
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
  xlab("")+ylim(0,1.3*100)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",rH2AX ~ final_dec,  
                   data = used_exp_rH2AX,paired = F )
NAN_plot
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.25*100)
)
#NAN_plot<- NAN_plot+coord_flip()
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="Fig4_rH2AX_IHC_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, height =8, units = 'cm', dpi = 600)
