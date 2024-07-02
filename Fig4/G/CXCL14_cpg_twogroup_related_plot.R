rm(list = ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
CXCL14_related <- readRDS('only_thalamus_for_reg.RDS')
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)
only_thalamus <- readRDS('only_thalamus_for_reg.RDS')

used_only_thalamus <- only_thalamus[,which(colnames(only_thalamus)
                      %in% c("cg01821923","cg10493994",
                      "cg03859527","CXCL14_tpm","CXCL14_cnv"))]
used_only_thalamus$CXCL14_tpm_log2 <- log2(used_only_thalamus$CXCL14_tpm+0.1)



thres2 <- quantile(used_only_thalamus$cg10493994,
                   seq(0.1,1,0.1))
used_only_thalamus$lab_cpg2 <- rep("l",length(used_only_thalamus$cg10493994))
used_only_thalamus$lab_cpg2[which(used_only_thalamus$cg10493994>
                                    thres2[9])] <- "g"
ggplot(used_only_thalamus,
       aes(lab_cpg2,CXCL14_tpm_log2,color=lab_cpg2))+
  geom_point()+geom_boxplot()+theme_classic()+
  stat_compare_means()+ggtitle('cg10493994')+xlab("")
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)
p1 <- ggplot(used_only_thalamus,aes(x=lab_cpg2,y=CXCL14_tpm_log2,color=lab_cpg2))+
  geom_boxplot(data=used_only_thalamus,
               aes(x=lab_cpg2,y=CXCL14_tpm_log2),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=used_only_thalamus,aes(x=lab_cpg2,y=CXCL14_tpm_log2, 
                                               color=lab_cpg2,fill=lab_cpg2),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme_classic()+stat_compare_means()+
  #geom_jitter(shape=16, position=position_jitter(0.2))+
  xlab("")+ylab('CXCL14 log2(tpm+0.1)')+
  scale_fill_manual(name=NULL,values = c(g='#FF0000',l='#9999FF'))+
  scale_color_manual(name=NULL,values =  c(g='#FF0000',l='#9999FF'))+
  ylim(0,11)+#'#f03b20'
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  theme(legend.position = "none")#+ylab("")
p1
pdf('Thalamus_like_cpg_two_group_boxplot.pdf',width=4,height = 4)
p1
dev.off()




