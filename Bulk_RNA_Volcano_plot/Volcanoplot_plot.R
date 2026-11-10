#setwd(path)
rm(list=ls())
print( getwd() )
info <- read.table('DE_genes_hg38_foldchange_plot.txt',sep="\t",header=T)
#info$lab<- "ns"
#info$lab[intersect(which(info$log2FoldChange< (-1)),
#                    which(info$wilcox_log2_tpm_01<0.05))] <- "down"
#info$lab[intersect(which(info$log2FoldChange> (1)),
#                    which(info$wilcox_log2_tpm_01<0.05))] <- "up"


library(ggplot2)
library(ggrepel)
p_thre=0.05
b_thre=1
p1 <- ggplot(info,aes((info$log2foldchange),
          -1*log10(as.numeric(info$wilcox_p)))) +    
  geom_point(aes(color = lab)) +           
  labs(title="",                                
       x=expression(paste(log[2],("fold change"),sep="")),#"Log2 (Fold change)", 
       y=expression(paste(-log[10],('P value'),sep="")))+theme_bw()+
  theme(panel.grid = element_blank())+
  geom_hline(yintercept=-1*log10(p_thre),
             linetype="dashed")+
  #geom_text(aes(0, -1*log10(p_thre), label = as.character(round(-1*log10(p_thre),3)), vjust = - 1))+
  geom_vline(xintercept =b_thre,linetype="dashed")+
  #geom_text(aes(b_thre,4,label =as.character(b_thre), hjust = - 1))+
  geom_vline(xintercept =-b_thre,linetype="dashed")+
  #geom_text(aes(-b_thre,4,label =as.character(-b_thre) , hjust = - 1))+#2.5
  geom_text_repel(label=ifelse(info$show=="yes" ,
                               info$sy ,""),
                  max.overlaps = 1000,                    # max overlaps; with many points some labels are hidden, raise this to keep them
                  size=3,                                  # font size
                  box.padding=unit(0.5,'lines'),           # padding around each label
                  point.padding=unit(0.1, 'lines'),
                  segment.color='black',                   # colour of the label segments
                  show.legend=FALSE)+
  scale_color_manual(values = c("#2b8cbe","grey","#FF3300"))+#"#ece7f2"
  labs(colour="")+theme(legend.position = "none")
pdf('Volcanoplot.pdf',width=8,height = 8)
p1
dev.off()