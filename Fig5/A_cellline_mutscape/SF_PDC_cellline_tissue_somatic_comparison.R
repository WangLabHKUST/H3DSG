rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
data <- readxl::read_excel('PDC_cellline_5_mutations.xlsx',sheet=1)
som <- matrix(0,nrow=length(data$mutation),ncol=length(colnames(data))-1)
rownames(som) <- data$mutation
colnames(som) <- colnames(data)[-1]
for(i in 2:dim(data)[2])
{
  for(j in 1:dim(data)[1])
  {
    if(is.na(data[j,i])==F)
    {
      if(data[j,i]=="missense")
      {
        som[j,i-1] <-1
      }else{
        if(data[j,i]=="frameshift")
        {
          som[j,i-1] <-2
        }else{
          if(data[j,i]=="stop_gained")
          {
            som[j,i-1] <-3
          }
        }
      }
    }
  }
}
som <- som[-3,]
library(ComplexHeatmap)
library(circlize)

figure_width <- 2.6
figure_height <- 0.28
size_numer_font_size <- 8#9.5#9.5#10.5
lwd_number <- 0.1




#row_ha <- rowAnnotation(freq = anno_barplot(apply(Mutat_matrix,1,function(x)length(which(x>0)))))
combine_mut <- Heatmap(som,
                       name="Som3",col=c('white', 'black','hotpink1',"#996633"),
                       column_split = c(rep(1,2),rep(2,2),rep(3,2),rep(4,3),rep(5,2)),
                       show_column_names = T,
                       row_split=c(rep(1,3),rep(2,3),rep(3,6)),
                       rect_gp = gpar(col = "grey20", lwd = lwd_number),
                       cluster_rows = FALSE,cluster_columns = FALSE,
                       width = unit(figure_width, "cm"), height = unit(0.25*14, "cm"),
                       row_names_side='left',
                       row_title = NULL,column_title = NULL,
                       row_names_gp = gpar(fontsize =size_numer_font_size, fontface="italic", col="black"),
                       column_names_gp = gpar(fontsize = size_numer_font_size,
                                              fontface="italic",col = c(rep(c("blue","red"),3),
                                               rep("blue",2),"red",rep("blue",2))),
                       heatmap_legend_param = list(title = "Mutation", at = c(1,2,3), 
                                                   labels = c("Missense","Frameshift","Stop_gain"
                                                              ),
                                                    border = "black"
                       )
)

combine_mut
pdf(paste('SP_PDC_5_cell_line_som_tissue.pdf',sep="/"),
    width=6,height = 6)
combine_mut
dev.off()