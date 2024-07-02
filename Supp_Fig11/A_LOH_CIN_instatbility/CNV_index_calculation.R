#https://cran.rstudio.com/web/packages/CINmetrics/vignettes/CINmetrics.html
rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
library(CINmetrics)
getwd()

GAB_seg <- read.table('Group_AB_seg.seg',sep="\t",header=T)
GAB_seg <- GAB_seg[,c(2:4,6,1)]
colnames(GAB_seg) <- c("Chromosome","Start","End","Segment_Mean","Sample")

cinmetrics.test_GAB <- CINmetrics(cnvData = GAB_seg)


#fraction.genome.test <- fga(cnvData = GAB_seg,segmentMean = 0.1)
clin_info <- read.csv('SP_36_clinic_with_purity_0504.csv',header=T)
lab <- rbind()
for(i in 1:dim(cinmetrics.test_GAB)[1])
{
  k1 <- which(paste0('P',clin_info$ID2...2) %in% cinmetrics.test_GAB$sample_id[i])
  if(length(k1)>0)
  {
    lab <- rbind(lab,c(cinmetrics.test_GAB$sample_id[i],
                       clin_info$final_dec[k1])) 
  }
}
cinmetrics.test_GAB$group <- lab[,2]
cinmetrics.test_GAB$group[which(cinmetrics.test_GAB$group %in% "YongG")] <- "Group1"
cinmetrics.test_GAB$group[which(cinmetrics.test_GAB$group %in% "OldG")] <- "Group2"
lab <- as.data.frame(lab)
rownames(lab) <- lab[,1]
modified.tai.test <- taiModified(cnvData = GAB_seg)
rownames(modified.tai.test) <- modified.tai.test$sample_id
modified.tai.test <- merge(modified.tai.test,lab,by="row.names",all=T)
modified.tai.test$modified_tai <- as.numeric(modified.tai.test$modified_tai)
modified.tai.test$V2[which(modified.tai.test$V2 %in% "YongG")] <- "Group1"
modified.tai.test$V2[which(modified.tai.test$V2 %in% "OldG")] <- "Group2"
colnames(modified.tai.test)[5]<- "group"
modified.tai.test$group <- as.factor(modified.tai.test$group)


cinmetrics.test_GAB$tai <- as.numeric(cinmetrics.test_GAB$tai)
cinmetrics.test_GAB$cna <- as.numeric(cinmetrics.test_GAB$cna)
cinmetrics.test_GAB$base_segments <- as.numeric(cinmetrics.test_GAB$base_segments)
cinmetrics.test_GAB$break_points <- as.numeric(cinmetrics.test_GAB$break_points)
cinmetrics.test_GAB$fga <- as.numeric(cinmetrics.test_GAB$fga)
cinmetrics.test_GAB$group <- as.factor(cinmetrics.test_GAB$group)
library(ggplot2)
library(ggrepel)
library(ggpubr)
p2 <- ggplot(cinmetrics.test_GAB, aes(x=group, y=cna,
       fill=group)) + #label=sample_id,
       geom_boxplot(width=0.35)+
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
       #geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
       stat_compare_means()+#ggtitle('Copy Number Aberration')+
       scale_fill_manual(values=c('#f03b20','#2b8cbe'))+
       theme_classic()+xlab("")+
       theme(axis.text=element_text(angle =90,size=10,color="black"),#,face="bold",face="bold"
        axis.text.x = element_text(angle =0,size=10,color="black"),
        axis.title=element_text(size=10))+ylab('Copy Number Aberration')+
  theme(legend.position = "none")+ylim(0,100)

p2_no_label<- p2+xlab("")+ylab("")+
  theme(axis.text.y = element_blank(),axis.text.x = element_blank())
pdf('CIN_CNA.pdf',width=4.8,
    height = 3)
ggarrange(p2, #p2_no_label, 
          labels = c("",""),
          ncol = 2, nrow = 1)
dev.off() 

p3 <- ggplot(cinmetrics.test_GAB, aes(x=group, y=base_segments,
       fill=group)) + #label=sample_id,
       geom_boxplot(width=0.35)+
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
       #geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
       stat_compare_means()+#ggtitle('Altered Base segments')+
       scale_fill_manual(values=c('#f03b20','#2b8cbe'))+
        theme_classic()+xlab("")+
        theme(axis.text=element_text(angle =90,size=10,color="black"),#,face="bold",face="bold"
        axis.text.x = element_text(angle =0,size=10,color="black"),
        axis.title=element_text(size=10))+theme(legend.position = "none")+
  ylim(0,2.1e+09)

p3_no_label<- p3+xlab("")+ylab("")+
  theme(axis.text.y = element_blank(),axis.text.x = element_blank())
pdf('CIN_base_segment.pdf',width=4.8,
    height = 3)
ggarrange(p3, #p3_no_label, 
          labels = c("",""),
          ncol = 2, nrow = 1)
dev.off() 

p4 <- ggplot(cinmetrics.test_GAB, aes(x=group, y=break_points,
      fill=group)) + # label=sample_id,
       geom_boxplot(width=0.35)+
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
       #geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
       stat_compare_means()+#ggtitle('Number of Break Points')+
       scale_fill_manual(values=c('#f03b20','#2b8cbe'))+
        theme_classic()+xlab("")+
        theme(axis.text=element_text(angle =90,size=10,color="black"),#,face="bold",face="bold"
        axis.text.x = element_text(angle =0,size=10,color="black"),
        axis.title=element_text(size=10))+
  theme(legend.position = "none")+ylim(0,260)
p4_no_label<- p4+xlab("")+ylab("")+
  theme(axis.text.y = element_blank(),axis.text.x = element_blank())
pdf('CIN_breakpoint.pdf',width=5.4,
    height = 3.2)
ggarrange(p4, #p4_no_label, 
          labels = c("",""),
          ncol = 2, nrow = 1)
dev.off() 

p5 <- ggplot(cinmetrics.test_GAB, aes(x=group, y=fga,
       fill=group)) + #label=sample_id,
       geom_boxplot(width=0.35)+
       geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)+
       #geom_text_repel(box.padding = 0.5, max.overlaps = Inf)+
       stat_compare_means()+#ggtitle('Fraction of Genome Altered')+
       scale_fill_manual(values=c('#f03b20','#2b8cbe'))+
       theme_classic()+xlab("")+
       theme(axis.text=element_text(angle =90,size=10,color="black"),#,face="bold",face="bold"
        axis.text.x = element_text(angle =0,size=10,color="black"),
        axis.title=element_text(size=10))+theme(legend.position = "none")+ylim(0,1)+
     ylab('Fraction of Genome Altered')

p5_no_label<- p5+xlab("")+ylab("")+
  theme(axis.text.y = element_blank(),axis.text.x = element_blank())
pdf('CIN_FGA.pdf',width=4.8,
    height = 3)
ggarrange(p5, #p5_no_label, 
          labels = c("",""),
          ncol = 2, nrow = 1)
dev.off() 


