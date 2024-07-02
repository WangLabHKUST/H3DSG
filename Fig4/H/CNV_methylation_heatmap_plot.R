rm(list = ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
library(ComplexHeatmap)
CNV <- readRDS('CNV_heatmap_plot.RDS')
library(copynumber)
library(GenomicRanges)
library(circlize)
chr_df = read.chromInfo()$df
chr_df = chr_df[chr_df$chr %in% c(paste0("chr", 1:22)),]#,'chrX','chrY'), ]
chr_gr = GRanges(seqnames = chr_df[, 1], ranges = IRanges(chr_df[, 2] + 1, chr_df[, 3]))
chr_gr
figure_width <- 14.8
figure_height <- 0.28
size_numer_font_size <- 9#9.5#9.5#10.5
lwd_number <- 0.1
split <- c(rep(1,16),rep(2,20))
library(EnrichedHeatmap)
chr_window = makeWindows(chr_gr, w = 3e6)
chr_window
chr = as.vector(seqnames(chr_window))
chr_level = paste0("chr", 1:22)
chr = factor(chr, levels = chr_level)
p1 <- Heatmap(as.matrix(CNV), name = "CNV", 
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        row_split =chr, cluster_rows = FALSE, show_column_dend = F,
        cluster_columns=F,
        column_split = split,cluster_column_slices = FALSE,
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10,
            fontface="italic", col="black"), border = TRUE,
        row_gap = unit(0, "points"))
getwd()
pdf('CNV_heatmap36_sample.pdf',width=5,height = 7)
p1
dev.off()

methyZscore <- readRDS('Methylation_nummat_data_z_score.RDS')
methySplit <- readRDS('DNA_methy_split_label.RDS')
p2 <- Heatmap(t(as.matrix(methyZscore)), name = "z_meth", 
        #col = colorRamp2(c(0, 0.2,0.4,0.6, 0.8,1), 
        #      c("#4575b4","#91bfdb","#e0f3f8","#fee090",
        #        "#fc8d59","#d73027")),
        row_split =chr, cluster_rows = FALSE, show_column_dend = F,
        cluster_columns=F,
        column_split = t(as.matrix(methySplit)),
        cluster_column_slices = FALSE,
        column_gap = unit(1, "points"),
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10,
                                               fontface="italic", col="black"), border = TRUE,
        border_gp = gpar(col = "grey30", lty = 1,lwd=0.4),
        row_gap = unit(0, "points"),
        show_column_names = F)

pdf('meth_heatmap36_sample.pdf',width=5,height = 7)
p2
dev.off()

all_index_info<- readRDS('CIN_methylation_correlation0124.RDS')
library(ggside)
library(tidyverse)
library(tidyquant)
library(circlize)
cytoband.file = system.file(package = "circlize",
                            "extdata", "cytoBand.txt")
cytoband.df = read.table(cytoband.file, colClasses = c("character", "numeric",
                                                       "numeric", "character", "character"),sep = "\t")
dim(cytoband.df[which(!(cytoband.df$V1 %in% c("chrX","chrY"))),])
cytoband.df_sel <- cytoband.df[which(!(cytoband.df$V1 %in% c("chrX","chrY"))),]
region_cor <- rbind() 
for(kkk in 1:dim(cytoband.df_sel)[1])
{
  only_1st_region <- all_index_info[intersect(intersect(which(all_index_info[,1] %in% cytoband.df_sel$V1[kkk]),
                                                        which(as.numeric(all_index_info[,2]) == cytoband.df_sel$V2[kkk])),
                                              which(as.numeric(all_index_info[,3]) == cytoband.df_sel$V3[kkk])),]
  if(dim(only_1st_region)[1]>0){
    only_1st_region <- as.data.frame(only_1st_region)
    colnames(only_1st_region)[c(6,7)] <- c("CIN",'Methy')
    uu <- cor.test(as.numeric(only_1st_region$Methy),
                   as.numeric(only_1st_region$CIN),method = "spearman")
    region_cor <- rbind(region_cor,
                        c(cytoband.df_sel$V1[kkk],cytoband.df_sel$V2[kkk],
                          cytoband.df_sel$V3[kkk],cytoband.df_sel$V4[kkk],
                          cytoband.df_sel$V5[kkk],uu$p.value,
                          uu[["estimate"]][["rho"]]))
  }
  
}
dim(region_cor)
region_cor <- as.data.frame(region_cor)
region_cor$V7 <- as.numeric(region_cor$V7)
region_cor <- region_cor[order(region_cor$V7,decreasing = F),]
region_cor$lab <- rep("No_sig",length(region_cor$V1))
region_cor$lab[which(as.numeric(region_cor$V6)<0.05)] <- "sig"
region_cor$lab <- as.factor(region_cor$lab)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
require(MASS)

sig_1st_region <- all_index_info[intersect(intersect(which(all_index_info[,1] %in% region_cor$V1[1]),
                                                     which(as.numeric(all_index_info[,2]) == as.numeric(region_cor$V2[1]))),
                                           which(as.numeric(all_index_info[,3]) == as.numeric(region_cor$V3[1]))),]

sig_1st_region <- as.data.frame(sig_1st_region)
colnames(sig_1st_region)[c(6,7)] <- c("CIN",'Methy')
sig_1st_region$Methy <- as.numeric(sig_1st_region$Methy)
sig_1st_region$CIN <- as.numeric(sig_1st_region$CIN)
p11 <- ggplot(region_cor, aes(x=V7)) + 
  #geom_density()+
  geom_histogram(aes(fill=lab),binwidth=0.02, colour="white")+#fill="white"#black#y=..density..,
  #geom_density(aes(fill=lab),lwd = 0.6,alpha=.2,colour ="magenta1")+#fill="#FF6666"2
  xlab("Correlation coefficient by spearman between CIN and methylation level")+
  theme_classic()+xlim(-0.7,0.3)+
  theme(axis.text=element_text(angle =0,size=10,color="black"),#,face="bold",face="bold"
        axis.title=element_text(size=10))+scale_fill_manual(values=c('skyblue',
                                                                     'blue3'),#steel
                                                            labels=c("ns","p<0.05"))+
  ylab("Freq")+#geom_vline(xintercept = -0.33, linetype="dashed", 
  #               color = "black", size=0.5)+ labs(fill=NULL)+xlab( expression(CM_cor))+
  theme(legend.position = c(0.8,0.6))+
  labs(fill=NULL)+xlab( expression(CM_cor))+
  theme(legend.position = c(0.8,0.6))+
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey", size=0.5)
#Density#grey
p22 <- ggplot(sig_1st_region,aes(x = Methy, y = CIN)) +
  geom_point(aes(color = Group),
             shape = 16)+theme_classic()+# alpha = .5,
  geom_smooth(method = lm, se=FALSE,color="purple")+#,linetype="dashed""grey50"
  stat_cor(label.y = 0.95,method = "spearman")+
  xlim(0.2,0.6)+ylim(0,1)+
  scale_color_manual(values = c("#FF3300","#2b8cbe"))+#"grey",
  #stat_regline_equation(label.y = 0.7)+#"#990000"
  theme(axis.text=element_text(angle =0,size=10,color="black"),#,face="bold",face="bold"
        axis.title=element_text(size=10))+xlab('Beta value of Region:chr4q25')+
  ylab('CIN of Region:chr4q25')+theme(legend.position = "none")
pdf('Cor_between_CIN_and_methylation.pdf',width=5.6,height = 5.6)
ggarrange(p22, p11, p11+theme(legend.position = "none"),
          #labels = c("A", "B", "C"),
          ncol = 2, nrow =2)
dev.off()
