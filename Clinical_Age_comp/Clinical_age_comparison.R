rm(list=ls())
if (!require('rstudioapi'))
  install.packages('rstudioapi')
if (!require('ggplot2'))
  install.packages('ggplot2')
if (!require('ggpubr'))
  install.packages('ggpubr')
library(rstudioapi)

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
info_plot <- read.table('Age_comp.txt',sep="\t",header=T)
info_plot <- as.data.frame(info_plot)
info_plot$Age <- as.numeric(info_plot$Age)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(ggrepel)
my_comparisons <- list( c("H3-Medulla", "H3-Pons"),
                        c("SP_group1", "SP_group2"),
                        c("H3-Medulla", "SP_group2"),
                        c("H3-Pons","SP_group1"))
p1 <- ggplot(data=info_plot,aes(x=reorder(final_dec,c(1:length(info_plot$M_ID))),
                           y=Age,fill=(final_dec)))+
  geom_boxplot(data=info_plot,width=0.5,fill="transparent")+
  geom_dotplot(binaxis='y', stackdir='center',dotsize=0.6)  +
  stat_compare_means(comparisons = my_comparisons)+
  theme_classic()+xlab("")+
  scale_fill_manual(values=c('#a6bddb','#feb24c','#990000',#'#f03b20',
                             '#2b8cbe'))+
  theme(axis.text=element_text(angle =90,size=10,color="black"),
        axis.title=element_text(size=10))+
  theme(legend.position="none")+ylab("Age at diagnosis (years)")
p1
pdf('Clinical_age_comparison.pdf',width=5,height = 5)
p1
dev.off()



