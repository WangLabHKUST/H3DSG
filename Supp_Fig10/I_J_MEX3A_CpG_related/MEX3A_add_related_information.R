rm(list = ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
library(MASS)
### correlation test between gene expression, methylation,
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(readxl)
cur_gene <- "MEX3A"
cur_exp_cnv_cpg_TF <- readRDS('MEX3A_TF_CpG_plot.RDS')
ptf_regress1<- ggplot(cur_exp_cnv_cpg_TF,aes(cg16549043,exp,fill="white"))+
  geom_point(aes(color=lab))+
  scale_color_manual(values=c('#f03b20','#2b8cbe'))+
  stat_smooth(method = "rlm", 
              geom = "smooth",se=FALSE,color="grey40")+
  xlab("Beta value of cg16549043")+ylab(paste0(cur_gene," expression"))+
  theme_classic()+
  theme(axis.text=element_text(angle =0,size=10,color="black"),#,face="bold",face="bold"
        axis.title=element_text(size=10))+ylim(0,1.1)+
  theme(legend.position = "none")+stat_cor(label.y = 1.04)
pdf(paste0(cur_gene,'_cg16549043_scatter_related_info.pdf'),
    width=3,height = 3)
ggarrange(ptf_regress1, #ptf_regress1, 
          labels = c(""),
          ncol = 1, nrow = 1)
dev.off()

library(ggbeeswarm)
NAN_plot <- ggplot(data=cur_exp_cnv_cpg_TF,
                   aes(x=lab,y=cg16549043)) + theme_classic() +
  geom_boxplot(data=cur_exp_cnv_cpg_TF,
               aes(x=lab,y=cg16549043),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=cur_exp_cnv_cpg_TF,aes(x=lab,y=cg16549043, 
                                                   color=lab, 
                                        fill=lab),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black'),axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,1.1)+
  theme(axis.text.x=element_blank())+
  theme(legend.position = "none")
the<-compare_means(method = "wilcox.test",cg16549043 ~ lab,  
                   data = cur_exp_cnv_cpg_TF,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(1.05)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="MEX3A_cg16549043_box_comparison.pdf",
       plot=figure_2,bg = 'white', width =6,
       height = 8, units = 'cm', dpi = 600)

 



