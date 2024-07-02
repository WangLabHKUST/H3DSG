rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
only_thalamus<- readRDS('only_thalamus_for_reg.RDS')
library(ggplot2)
library(ggpubr)


p2 <- ggplot(only_thalamus,aes(x=cg10493994,
            y=log2(CXCL14_tpm+0.1)))+geom_point()+
  theme_classic()+ 
  geom_smooth(method='lm',se = FALSE,color="grey60")+
  ylab('CXCL14 expression log2(TPM+0.1)')+
  stat_cor(label.y = 10.5)+ylim(0,11)+
  theme(legend.position = "none")
pdf('Thalamus_cg10493994_CXCL14_neg.pdf',width=6,height = 6)
p2
dev.off()





meth_exp_combine_heat<- readRDS('meth_exp_combine_heat.RDS')
#meth_exp_combine_heat <-as.data.frame((meth_exp_combine_heat))
meth_exp_combine_heat$cg01821923 <- as.numeric(meth_exp_combine_heat$cg01821923 )
meth_exp_combine_heat$cg10493994 <- as.numeric(meth_exp_combine_heat$cg10493994 )
meth_exp_combine_heat$cg03859527 <- as.numeric(meth_exp_combine_heat$cg03859527 )
meth_exp_combine_heat$CXCL14_tpm <- (as.numeric(meth_exp_combine_heat$CXCL14_tpm ))
#meth_exp_combine_heat$CXCL14_tpm <- log2(as.numeric(meth_exp_combine_heat$CXCL14_tpm )+0.1)
#meth_exp_combine_heat$EGFR_tpm <- as.numeric(meth_exp_combine_heat$EGFR_tpm)
#meth_exp_combine_heat$EGFR_tpm <- log2(as.numeric(meth_exp_combine_heat$EGFR_tpm )+0.1)
meth_exp_combine_heat$CXCL14_cnv <- as.numeric(meth_exp_combine_heat$CXCL14_cnv)
split <- meth_exp_combine_heat$group.y
split[which(split%in% "Pons")] <- 1
split[which(split%in% "Thalamus")] <- 2
PP1 <- ggplot(meth_exp_combine_heat,
       aes(CXCL14_cnv,CXCL14_tpm))+
  geom_point(aes(color=group.y))+
  geom_smooth(method = lm, se = FALSE,
              color="grey60")+theme_classic()+
  scale_color_manual(name=NULL,
        values = c(Pons='#f03b20',Thalamus='#2b8cbe'))+
  xlab('CXCL14 CNV (.segmean)')+ylab('CXCL14 (Log2(TPM+0.1))')+
  theme(legend.position = "none")+stat_cor(label.y = 9)
pdf('CXCL14_CNV_exp_regression.pdf',
    width = 4,height = 4)
PP1
dev.off()

if(F)
{
table(split)
library(ComplexHeatmap)
p1_heat <- Heatmap(as.matrix(t(meth_exp_combine_heat[,c(1:3)])), name = "meth", 
        cluster_rows = FALSE, show_column_dend = F,
        cluster_columns=F,
        column_split = (as.matrix(split)),
        cluster_column_slices = FALSE,
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10,
                                               fontface="italic", col="black"), border = TRUE,
        row_gap = unit(0, "points"))

p2_heat<- Heatmap(as.matrix(t(meth_exp_combine_heat[,c(4)])), name = "exp", 
        cluster_rows = FALSE, show_column_dend = F,
        cluster_columns=F,row_labels =  "CXCL14_exp",row_title_side = 'right',
        column_split = (as.matrix(split)),
        cluster_column_slices = FALSE,
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10,
        fontface="italic", col="black"), border = TRUE,
        row_gap = unit(0, "points"))

p3_heat <- Heatmap(as.matrix(t(meth_exp_combine_heat[,c(7)])), name = "cnv", 
        cluster_rows = FALSE, show_column_dend = F,
        
        cluster_columns=F,row_labels =  "CXCL14_CNV",row_title_side = 'right',
        column_split = (as.matrix(split)),
        cluster_column_slices = FALSE,
        row_title_rot = 0, row_title_gp = gpar(fontsize = 10,
        fontface="italic", col="black"), border = TRUE,
        row_gap = unit(0, "points"))

p3_heat%v%p1_heat%v%p2_heat
pdf('CXCL14_cpg_CNV_heatmap.pdf',width = 6,height = 5)
p3_heat%v%p1_heat%v%p2_heat
dev.off()
p3 <- ggplot(only_thalamus,aes(x=cg03859527,
                               y=log2(CXCL14_tpm+0.1)))+geom_point()+
  theme_classic()+ 
  geom_smooth(method='lm',se = FALSE,color="grey60")+
  ylab('CXCL14 expression log2(TPM+0.1)')+
  stat_cor(label.y = 10.5)+ylim(0,11)+
  theme(legend.position = "none")
pdf('Thalamus_cg03859527_CXCL14_neg.pdf',width=6,height = 6)
p3
dev.off()
p1 <- ggplot(only_thalamus,aes(x=cg01821923,
                               y=log2(CXCL14_tpm+0.1)))+geom_point()+
  theme_classic()+ 
  geom_smooth(method='lm',se = FALSE,color="grey60")+
  ylab('CXCL14 expression log2(TPM+0.1)')+
  stat_cor(label.y = 10.5)+ylim(0,11)+
  theme(legend.position = "none")
pdf('Thalamus_cg01821923_CXCL14_neg.pdf',width=6,height = 6)
p1
dev.off()
p1
}

