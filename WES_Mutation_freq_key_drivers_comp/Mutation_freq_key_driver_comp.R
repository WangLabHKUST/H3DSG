rm(list=ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
if (!require('ggplot2')) 
  install.packages('ggplot2')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
used_line_plot_somatic <- read.table('Genomic_frequency_comprison_H3.3_alteration_0523.txt',
            sep="\t",header=T)
used_line_plot_somatic$Type <- gsub("_freq","",
                              used_line_plot_somatic$Type)

library(ggplot2)
pdf('Mutation_freq_key_driver_comp.pdf',
    width=5,height = 3.5)
ggplot(data=used_line_plot_somatic, aes(x=reorder(Gene,
              c(1:dim(used_line_plot_somatic)[1])), y=freq, group=Type,
              color=Type)) +
  geom_line(linetype = "solid",size=0.6)+#"dashed" 0.8
  geom_point(aes(shape=shape))+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1#,face="bold"
                                   ))+
  xlab("")+ylab("Somatic mutation freq")+
  scale_shape_manual(values=c(3,17))+
  theme(axis.text=element_text(size=10),#,face="bold",face="bold",face="bold"
        axis.title=element_text(size=10))+#,face="bold"
  #scale_fill_discrete(name = "New Legend Title")+
  guides(color = guide_legend(title = "Type"))+
  scale_color_manual(values=c('#a6bddb','#feb24c','#990000',
                              '#2b8cbe'))#'#f03b20',
dev.off()



