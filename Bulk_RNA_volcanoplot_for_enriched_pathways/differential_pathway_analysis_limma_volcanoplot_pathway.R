rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(getwd())
###getting from single sample GSEA
### similar to differential gene analysis 

c2.cp.kegg.medicus <- read.csv('gsave_matrix_ssgsea_info_c2.cp.kegg_medicus.csv',
                               header=T,row.names = 1)
dim(c2.cp.kegg.medicus)
c2.cp.kegg.legacy <- read.csv('gsave_matrix_ssgsea_info_c2.cp.kegg_legacy.csv',
                              header=T,row.names = 1)
dim(c2.cp.kegg.legacy)
c2.cp.reactome <- read.csv('gsave_matrix_ssgsea_info_c2.cp.reactome.csv',
                           header=T,row.names = 1)
dim(c2.cp.reactome)
c5.gp.bp <- read.csv('gsave_matrix_ssgsea_info_c5.gp.bp.csv',
                     header=T,row.names = 1)
c2.hallmark <- read.csv('gsave_matrix_ssgsea_info.csv',
                        header=T,row.names = 1)
dim(c2.hallmark)
all_score <- rbind(c2.cp.kegg.medicus,
             rbind(rbind(c2.cp.kegg.legacy,c2.cp.reactome),
                  c2.hallmark))#)#rbind(c5.gp.bp,
#### label load
data_label <- readRDS('label_for_twoGroups_SP36_1219.RDS')
Pons_loc <- data_label$RNA[which(data_label$group %in% "YongG")]
Thalamus_loc <- data_label$RNA[which(data_label$group %in% "OldG")]

pons_score <- all_score[,which(colnames(all_score)
                            %in% Pons_loc)]
thalamus_score <-  all_score[,which(colnames(all_score)
                              %in% Thalamus_loc)]
score_combine <- as.matrix(cbind(pons_score,
                       thalamus_score))
group_info <- c(rep("pons",16),
                rep("thalamus",20))
library(limma)
design <- model.matrix(~0+factor(group_info))
colnames(design) <- levels(factor(group_info))
rownames(design) <- colnames(score_combine)
contrast.matrix <- makeContrasts(pons-thalamus,levels = design)
fit <- lmFit(score_combine,design) #非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DEG <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DEG <- na.omit(DEG)
DEG$regulate <- rep("unchanged",length(DEG$logFC))
DEG$regulate[intersect(which(DEG$logFC>0.05),
                       which(DEG$adj.P.Val<0.05))] <- "up-regulated"

DEG$regulate[intersect(which(DEG$logFC < (-0.05)),
                       which(DEG$adj.P.Val<0.05))] <- "down-regulated"

table(DEG$regulate)
rownames(DEG)[which(DEG$regulate %in% "down-regulated")]
rownames(DEG)[which(DEG$regulate  %in% "up-regulated")]
#write.table(DEG,"differential_pathway_analysis_limma.txt",
#            sep="\t",row.names = T,quote=F)
#library(edgeR)
DEG$label <- rep("no",length(DEG$logFC))
up_pathway <- c("KEGG_MEDICUS_REFERENCE_PRE_IC_FORMATION",
                #"KEGG_MEDICUS_REFERENCE_TRAIP_DEPENDENT_REPLISOME_DISASSEMBLY",
                #"REACTOME_UNWINDING_OF_DNA",
                #"KEGG_MEDICUS_REFERENCE_ORIGIN_UNWINDING_AND_ELONGATION",
                "KEGG_MEDICUS_REFERENCE_DNA_REPLICATION_TERMINATION",
                #"KEGG_MEDICUS_REFERENCE_BREAK_INDUCED_REPLICATION",
                "KEGG_MEDICUS_REFERENCE_DNA_REPLICATION_LICENSING"
                )
down_pathway <- c(#"KEGG_MEDICUS_REFERENCE_REGULATION_OF_COMPLEMENT_CASCADE_MAC_INHIBITION",
                  "HALLMARK_BILE_ACID_METABOLISM",
                  #"KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION",
                  #"KEGG_ASTHMA","KEGG_MEDICUS_REFERENCE_NLRC4_INFLAMMASOME_SIGNALING_PATHWAY",
                  #"REACTOME_INTERLEUKIN_1_PROCESSING",
                   "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES",
                    "KEGG_MEDICUS_REFERENCE_NLRP1_INFLAMMASOME_SIGNALING_PATHWAY")
DEG$label[which(rownames(DEG) %in% c(up_pathway,down_pathway))] <- "yes"

library(ggplot2)
library(ggrepel)
p_thre=0.05
b_thre=0.05
####other color
p1 <- ggplot(DEG,aes(logFC,
                -1*log10(as.numeric(DEG$adj.P.Val)))) +    
  geom_point(aes(color = DEG$regulate),size=1)+theme_bw()+
  theme(panel.grid = element_blank())+ylim(0,6)+xlim(-0.2,0.2)+
  theme(legend.position = "none")+xlab("log2 fold change")+
  ylab("-log10 adjusted p-value")+
  scale_color_manual(values = c("#2b8cbe","grey","#FF3300"))+#"#FF3300","#2b8cbe"
  geom_vline(xintercept =b_thre,linetype="dashed")+
  #geom_text(aes(b_thre,4,label =as.character(b_thre), hjust = - 1))+
  geom_vline(xintercept =-b_thre,linetype="dashed")+
  geom_hline(yintercept=-1*log10(p_thre),
             linetype="dashed")+
  geom_text_repel(label=ifelse(DEG$label=="yes",
                               rownames(DEG),""),
                  max.overlaps = 1000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
                  size=2,                                  # 字体大小
                  box.padding=unit(0.35,'lines'),           # 标记的边距
                  point.padding=unit(0.3, 'lines'), 
                  #segment.color='black',                   # 标记线条的颜色
                  show.legend=FALSE)
pdf('Differential_pathway_analysis_volcanoplot.pdf',
    width=6,height = 6)
p1
dev.off()
  

