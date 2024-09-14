
rm(list = ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
#LOH_index <- read.csv('SP_36_LOH_statistics.csv',header=T)
#LOH_index$group[which(LOH_index$group %in% "YongG")] <- "Group1"
#LOH_index$group[which(LOH_index$group %in% "OldG")] <- "Group2"
#LOH_index$group <- as.factor(LOH_index$group)
#LOH_index <- LOH_index[,-c(1)]
#write.table(LOH_index,'LOH_index_plot.txt',sep="\t",
#            row.names = F,quote=F)
LOH_index <- read.table('LOH_index_plot.txt',sep="\t",
                         header=T)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
NAN_plot <- ggplot(data=LOH_index,
                   aes(x=group,y=fra)) + theme_classic() +
  geom_boxplot(data=LOH_index,
               aes(x=group,y=fra),width = 0.6,size=0.5,fill="transparent",
               outlier.colour = "white")+
  geom_quasirandom(data=LOH_index,aes(x=group,y=fra, color=group, 
                                                   fill=group),
                   width = 0.25,size=1.5,alpha=0.7,stroke=0.8, varwidth = T)+
  theme(panel.background=element_rect(fill='transparent',color='transparent'),
        plot.margin=unit(c(2,2,2,2),'lines'),
        plot.title=element_text(size=34,vjust=0.5,hjust=0.5,face='bold.italic',color='transparent'),text=element_text(size=14,face='bold'),
        legend.key.width=unit(0.6,'cm'),legend.key.height=unit(0.6,'cm'),legend.position='right',legend.text=element_text(size=14,hjust=0,face='bold'),
        axis.text.x=element_text(size=12,face='bold',color='black',
                                 angle = 90, vjust = 0.5),
        axis.text.y=element_text(size=12,face='bold',color='black'),
        axis.title.x=element_text(size=14,face='plain',color='black'),
        axis.title.y=element_text(size=14,face='plain',color='black'))+
  scale_fill_manual(name=NULL,values = c(Group1='#f03b20',Group2='#2b8cbe'))+
  scale_color_manual(name=NULL,values =  c(Group1='#f03b20',Group2='#2b8cbe'))+
  xlab("")+ylim(0,0.6)+ylab('CNVs/LOH')+
  #theme(axis.text.x=element_blank())+
  theme(legend.position = "none")
#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
the<-compare_means(method = "wilcox.test",fra ~ group,  
                   data = LOH_index,paired = F )
NAN_plot <- NAN_plot +stat_pvalue_manual(
  the, label = NULL, size = 3,linetype = 1,
  y.position = c(0.58)
)
figure_2<-rbind(ggplotGrob(NAN_plot),size="last")
ggsave(file="LOH_index_comparison.pdf",
       plot=figure_2,bg = 'white', width =6, 
       height = 8, units = 'cm', dpi = 600)
