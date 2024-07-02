rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
AUP_ROC <- as.data.frame(read_excel("AUC_ROC_LOOCV_plot_barPlot_index.xlsx",sheet=1))
AUP_ROC <- AUP_ROC[-c(5,6,7),]
AUP_ROC$AUC <- as.numeric(AUP_ROC$AUC)
AUP_ROC$AUPRC <- as.numeric(AUP_ROC$AUPRC)
AUP_ROC_auc <- as.data.frame(AUP_ROC)
AUP_ROC_auc <- AUP_ROC_auc[order(AUP_ROC_auc$AUC),]
AUP_ROC_auc <- rbind(AUP_ROC_auc,AUP_ROC_auc[,c(1,2,4,3)])
AUP_ROC_auc$type <- c(rep("AUC",12),
                      rep("AUPRC",12))
AUP_ROC_auc$value <- c(AUP_ROC_auc$AUC[1:12],
                       AUP_ROC_auc$AUPRC[1:12])
#AUP_ROC_auprc <- as.data.frame(AUP_ROC)
#AUP_ROC_auprc <- AUP_ROC_auprc[order(AUP_ROC_auprc$AUPRC),]
library(ggplot2)
dim(AUP_ROC_auc)
# Basic barplot
p1<-ggplot(data=AUP_ROC_auc, aes(x=reorder(Label,c(1:24)), 
                                 y=value,fill=type)) +
  geom_bar(stat="identity",position="dodge",width=0.75,color="black")+theme_classic2()+
  scale_x_discrete(guide = guide_axis(angle = 90))+xlab("")+
  theme(axis.text=element_text(angle =0,size=10,color="black"),#,face="bold",face="bold"
        axis.title=element_text(size=10))+#ylim(0.8,1)
  scale_fill_manual(name="Index",values=c("orange2", "steelblue"))+ylab("Value")+
  xlab("")
p1#"#999999", "#E69F00"

pdf('AUC_AUPRC_plot.pdf',width = 5,height = 5)
p1
dev.off()