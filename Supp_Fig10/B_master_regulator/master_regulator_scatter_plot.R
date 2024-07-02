
rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )
library(ggplot2)
library(ggrepel)
RegNet_set_info_plot <- readRDS("RegNet_set_info_plot.RDS")
Set_by_sel_info<- readRDS("Set_by_sel_info_for_our_method.RDS")
fig1 <- ggplot()+ theme_classic()+
  xlab(expression("TF activity difference in two groups"))+
  ylab(expression("TF expression difference in two groups"))
fig1
fig1<- fig1+geom_point(data=RegNet_set_info_plot,aes(x=MS_A_minus_B,y=Exp_A_minus_B),color="grey")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())+
  geom_point(data=Set_by_sel_info,aes(x=MS_A_minus_B,y=Exp_A_minus_B,fill =Group)
             ,size=round(abs(Set_by_sel_info$NES)*2.5),shape = 21,
             colour = "#696969")+
  scale_fill_manual(values=c("#FF3300","#2b8cbe"))+
  geom_text_repel(aes(x=Set_by_sel_info$MS_A_minus_B[union(which(Set_by_sel_info$MS_A_minus_B < (-2)),
                                                           which(Set_by_sel_info$MS_A_minus_B > 2))],
                      y=Set_by_sel_info$Exp_A_minus_B[union(which(Set_by_sel_info$MS_A_minus_B < (-2)),
                                                    which(Set_by_sel_info$MS_A_minus_B > 2))],
                      label=Set_by_sel_info$TF[union(which(Set_by_sel_info$MS_A_minus_B<(-2)),
                                                     which(Set_by_sel_info$MS_A_minus_B>2))]),
                  #color=exp_diff$lab[union(which(exp_diff$NES< (-1.53)),
                  #                         which(exp_diff$NES>1.8))]),
                  max.overlaps = 1000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
                  size=6,                                  # 字体大小
                  box.padding=unit(0.5,'lines'),           # 标记的边距
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',                   # 标记线条的颜色
                  show.legend=FALSE)+theme(legend.position="no")+
  labs(color='TF types') +
  theme(axis.text.x = element_blank(),#element_text(angle =0, vjust = 0.5, hjust=1,color="black"),
        axis.text.y = element_blank(),#element_text(angle = 0, vjust = 0.5, hjust=1,color="black"),
        axis.ticks = element_blank(),
        text = element_text(size = 16,color="black"))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "grey")

fig1
pdf('master_regulator_scatter_plot.pdf',width=7,height = 7)
fig1
dev.off()


