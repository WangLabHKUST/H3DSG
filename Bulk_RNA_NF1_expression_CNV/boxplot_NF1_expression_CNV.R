rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
library(ggplot2)
library(ggpubr)
#library("tidyverse")
#library(gapminder)
data <- readxl::read_excel("NF1_seg_exp_hg38.xlsx",sheet=1)
#ggboxplot(data, x = "subgroup", y = "NF1_exp",
#          color = "NF1_del_loss_group", palette = "jco",
#          add = "jitter")+
#  stat_compare_means()+xlab("")+theme_classic()
p1 <- ggplot(data, aes(x=subgroup, y=NF1_exp, fill=NF1_del_loss_group)) + 
  geom_boxplot()+theme_classic()+
  stat_compare_means()+xlab("")#method = "wilcox.test"


pdf('NF1_exp_cnv_boxplot.pdf',
    width=6,height = 6)
p1
dev.off()
  

