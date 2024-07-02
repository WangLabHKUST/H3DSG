rm(list = ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
CXCL14_related <- readRDS('CXCL14_related.RDS')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)

#### boxplot for CNV
p1 <- ggplot(CXCL14_related,aes(x=group,y=CXCL14_cnv,color=group))+
  geom_boxplot(data=CXCL14_related,
               aes(x=group,y=CXCL14_cnv),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=CXCL14_related,aes(x=group,y=CXCL14_cnv, color=group,fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme_classic()+stat_compare_means()+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("")+ylab('CXCL14 (.segmean value)')+
  scale_fill_manual(name=NULL,values = c(Pons='#f03b20',Thalamus='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Pons='#f03b20',Thalamus='#2b8cbe'))+
  ylim(-1,0.25)+#'#f03b20'
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  theme(legend.position = "none")#+ylab("")
pdf('CXCL14_CNV_segmean_boxplot.pdf',width=4,height = 4)
p1
dev.off()
  
#### scatter
p2 <- ggplot(CXCL14_related,aes(x=CXCL14_cnv,
                          y=log2(tpm+0.1)))+geom_point(aes(color=group))+
  theme_classic()+ #geom_text_repel(aes(label=X))+
  scale_fill_manual(name=NULL,values = c(Pons='#f03b20',Thalamus='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Pons='#f03b20',Thalamus='#2b8cbe'))+
  geom_smooth(method='lm',se = FALSE,color="grey60")+ylab('CXCL14 expression log2(TPM+0.1)')+
  xlab('CXCL14 CNV (.segmean)')+
  stat_cor(label.y = 9)+theme(legend.position = "none") #'#f03b20'#990000

pdf('CXCL14_CNV_segmean_correlation_gene_expr.pdf',width=4,height = 4)
p2
dev.off()


