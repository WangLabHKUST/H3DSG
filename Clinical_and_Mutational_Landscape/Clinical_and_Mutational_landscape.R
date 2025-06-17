rm(list=ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
source('./required_packages.R')
library(ComplexHeatmap)
library(circlize)
figure_width <- 6
figure_height <- 0.28
size_numer_font_size <- 9
lwd_number <- 0.2


data_info = read.csv('Data_clinical_info.csv')
rownames(data_info) <- data_info[,1]
data_info <- data_info[,-1]

############# Fig 1. b ##############
###### Sequencing data
wes_matrix <- (as.matrix(data_info['WES',]))
wes_matrix <- t(apply(wes_matrix,1,function(x) as.numeric(x)))
colnames(wes_matrix) <- colnames(data_info)

WES_plot <- Heatmap(wes_matrix,
                    rect_gp = gpar(col = "black", lwd =lwd_number),#'goldenrod1','darkorchid',
                    col = c('#a7cae2'),#'white',"plum3",
                    row_labels = "WES",
                    cluster_rows = FALSE,cluster_columns = FALSE,
                    width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                    row_names_side='left',
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                    column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                    heatmap_legend_param = list(title = "WES",
                                                at = c(1), labels = c("Available"),border = "black",ncol=4),
                    show_heatmap_legend =T,border = T)
WES_plot
### 
methy_matrix <- (as.matrix(data_info['DNA_methylation',]))
methy_matrix <- t(apply(methy_matrix,1,function(x) as.numeric(x)))
colnames(methy_matrix) <- colnames(data_info)

methy_plot <- Heatmap(methy_matrix,
                      rect_gp = gpar(col = "black", lwd =lwd_number),
                      col = c("#f0c9c3"),
                      row_labels = "Methylation",
                      cluster_rows = FALSE,cluster_columns = FALSE,
                      width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                      row_names_side='left',
                      show_column_names = F,
                      row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                      column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                      heatmap_legend_param = list(title = "Methylation", at = c(1), 
                                                  labels = c("Available"),border = "black",ncol=3),
                      show_heatmap_legend =T,border = T)
methy_plot

RNA_matrix <- (as.matrix(data_info['RNA',]))
RNA_matrix <- t(apply(RNA_matrix,1,function(x) as.numeric(x)))
colnames(RNA_matrix) <- colnames(data_info)
RNA_plot <- Heatmap(RNA_matrix,
                    rect_gp = gpar(col = "black", lwd =lwd_number),#'goldenrod1','darkorchid',
                    col = c('#b1e18d'),#'white',
                    row_labels = "RNA-seq",
                    cluster_rows = FALSE,cluster_columns = FALSE,
                    width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                    row_names_side='left',
                    show_column_names = F,
                    row_names_gp = gpar(fontsize = size_numer_font_size, col="black"),#fontface="italic", 
                    column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),#fontface="italic",
        
                    heatmap_legend_param = list(title = "RNA",
                                                at = c(1), labels = c("Available"),border = "black",ncol=4),
                    show_heatmap_legend =T,border = T)

RNA_plot


mut_load_matrix <-(as.matrix(data_info['Count',]))
mut_load_matrix <- t(apply(mut_load_matrix,1,function(x) as.numeric(x)))
colnames(mut_load_matrix) <- colnames(data_info)
rownames(mut_load_matrix) <- "Mut_load"

col_ha <- HeatmapAnnotation(bar1 = anno_barplot(as.numeric(mut_load_matrix), 
                                                gp = gpar(fill = "grey"),border = F),
                            show_annotation_name = F)


location_matrix <- (as.matrix(data_info[10:29,]))
location_matrix <- apply(as.matrix(location_matrix),2,
                    function(x) as.numeric(x))
rownames(location_matrix) <- rownames(data_info)[10:29]
Row_sum <- apply(as.matrix(location_matrix),1,function(x) sum(as.numeric(x)))
loc_sel <- which(Row_sum>0)
location_matrix <- location_matrix[loc_sel,]

foc_loc <- apply(location_matrix,2,function(x) sum(as.numeric(x)))
col_ha_loc <- HeatmapAnnotation(freq=anno_barplot(foc_loc,axis = T,
              gp = gpar(fill = "blueviolet"),annotation_height = 1,
              bar_width = 0.35,border = F,axis_param=list(gp=gpar(fontsize = 4))),
              show_annotation_name = FALSE) 
foo_loc= apply(location_matrix,1,function(x) sum(as.numeric(x)))
row_ha_loc <- rowAnnotation(freq = anno_barplot(foo_loc,axis = T,
                bar_width = 0.35,gp = gpar(fill = "blueviolet"),
                border = F,axis_param=list(gp=gpar(fontsize = 4))),
                show_annotation_name = FALSE)
location <- Heatmap(location_matrix,
            name="Location",col=c('white', '#fc9272'),
            rect_gp = gpar(col = "black", lwd =lwd_number),
            cluster_rows = FALSE,cluster_columns = FALSE,
            width = unit(figure_width, "cm"), 
            height = unit(figure_height*10, "cm"),
            row_names_side='left',show_column_names = F,
            row_names_gp = gpar(fontsize = 3, col="black"),
            column_names_gp = gpar(fontsize = 3,col="black"),
            show_heatmap_legend = F,top_annotation=col_ha_loc,
            right_annotation = row_ha_loc,border = T)

loc_sim_matrix <- (as.matrix(data_info['loc',]))
loc_sim_matrix <- t(apply(loc_sim_matrix,1,function(x) as.numeric(x)))
colnames(loc_sim_matrix) <- colnames(data_info)
rownames(loc_sim_matrix) <- "loc"
colors_loc_sim <- c("lightblue4","lightsteelblue")
loc_sim_plot <- Heatmap(loc_sim_matrix,
                        name="Location_sim",col= colors_loc_sim,
                        rect_gp = gpar(col ="black", lwd =lwd_number),
                        row_labels = "Location",
                        show_column_names = F,
                        cluster_rows = FALSE,cluster_columns = FALSE,
                        width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                        row_names_side='left',
                        row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                        column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                        heatmap_legend_param = list(title = "Location", at = c(2,1), 
                                                    labels = c("Upper", "Lower"),ncol=2,
                                                    border = "black"),
                        show_heatmap_legend = T,border = T)
fig1b = WES_plot%v%methy_plot%v%RNA_plot%v%col_ha%v%location
fig1b


############# Fig 1. c ##############
data_info_ordered = data_info[,order(colnames(data_info))]
age_matrix <- (as.matrix(data_info_ordered['Age',]))
age_matrix <- t(apply(age_matrix,1,function(x) as.numeric(x)))
colnames(age_matrix) <- colnames(data_info_ordered)

Age_plot <- Heatmap(age_matrix,
                    rect_gp = gpar(col = "black", lwd =lwd_number),
                    col = colorRamp2(c(0, 20,40, 60), 
                                     c('#ffffcc','#c2e699','#78c679',"#238443")),
                    row_labels = "Age",
                    cluster_rows = FALSE,cluster_columns = FALSE,
                    width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                    row_names_side='left',show_column_names = T,
                    row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                    #column_names_gp = gpar(size=5),
                    column_names_gp = gpar(fontsize = 5,col="black"),
                    heatmap_legend_param = list(title = "Age", border = "black"),border = T
)
Age_plot
gender_matrix <- (as.matrix(data_info_ordered['Gender',]))
gender_matrix <- t(apply(gender_matrix,1,function(x) as.numeric(x)))
colnames(gender_matrix) <- colnames(data_info_ordered)

Gender_plot <- Heatmap(gender_matrix,
                       rect_gp = gpar(col = "black", lwd =lwd_number),#'goldenrod1','darkorchid',
                       col = c('coral','cyan4'),#,'steelblue1','sandybrown'
                       row_labels = "Sex",#purple3
                       cluster_rows = FALSE,cluster_columns = FALSE,
                       width = unit(figure_width, "cm"),
                       height = unit(figure_height, "cm"),
                       row_names_side='left',
                       show_column_names = F,# fontface="italic",fontface="italic",
                       row_names_gp = gpar(fontsize = size_numer_font_size, col="black"),
                       column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                       heatmap_legend_param = list(title = "Sex", at = c(2,3), labels = c("Female","Male"),
                                                   border = "black",ncol=2,
                                                   row_gap = unit(5, "mm")),
                       show_heatmap_legend = T,border = T
)

states_matrix <- (as.matrix(data_info_ordered['Status',]))
states_matrix <- t(apply(states_matrix,1,function(x) as.numeric(x)))
colnames(states_matrix) <- colnames(data_info_ordered)
rownames(states_matrix) <- "Sensor"

surv_states <- Heatmap((states_matrix),
                       rect_gp = gpar(col = "black", lwd =lwd_number),#'goldenrod1','darkorchid',
                       col = c('yellow','purple4'),#'white',
                       row_labels = "Vital status",
                       cluster_rows = FALSE,cluster_columns = FALSE,
                       width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                       row_names_side='left',
                       show_column_names = F,
                       row_names_gp = gpar(fontsize = size_numer_font_size, col="black"),#fontface="italic", 
                       column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),#fontface="italic",
                       heatmap_legend_param = list(title = "Vital status", at = c(0,1), labels = c("Living",
                                                                                                   "Decreased"),
                                                   border = "black",ncol=2),#top_annotation=row_ha_surv,
                       show_heatmap_legend =T,border = T
)
surv_states

Time_matrix <- (as.matrix(data_info_ordered['Months',]))
Time_matrix <- t(apply(Time_matrix,1,function(x) as.numeric(x)))
colnames(Time_matrix) <- colnames(data_info_ordered)
rownames(Time_matrix) <- "Time"

surv_time <- Heatmap((Time_matrix),
                     rect_gp = gpar(col = "black", lwd =lwd_number),
                     col = colorRamp2(quantile(Time_matrix),
                                      c("#f1eef6","#bdc9e1","#74a9cf","#2b8cbe","#045a8d")),
                     row_labels = "Survival months",
                     cluster_rows = FALSE,cluster_columns = FALSE,
                     width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                     row_names_side='left',
                     show_column_names = F,
                     row_names_gp = gpar(fontsize = size_numer_font_size, col="black"),#fontface="italic", 
                     column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),#fontface="italic",
                     heatmap_legend_param = list(title = "Survival months",border = "black"),
                     show_heatmap_legend =T,border = T
)
surv_time


ki67_matrix <- (as.matrix(data_info_ordered['KI67',]))
ki67_matrix <- t(apply(ki67_matrix,1,function(x) as.numeric(x)))
colnames(ki67_matrix) <- colnames(data_info_ordered)
rownames(ki67_matrix) <- "KI67"

ki67_plot <-  Heatmap(ki67_matrix,
                 rect_gp = gpar(col = "black", lwd =lwd_number),
                 col = colorRamp2(c(0, 25,50, 75), 
                                  c('#fee8c3','#ffcfa1','#faada5',"#f66b5d")),
                 row_labels = "Ki67",
                 cluster_rows = FALSE,cluster_columns = FALSE,
                 width = unit(figure_width, "cm"), height = unit(figure_height, "cm"),
                 row_names_side='left',show_column_names = T,
                 row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                 column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
                 heatmap_legend_param = list(title = "Ki67", border = "black"),border = T
)
ki67_plot
idx1 = which(rownames(data_info_ordered)=='H3-3A')
idx2 = which(rownames(data_info_ordered) == 'FGFR1')
cnvmut_matrix <- (as.matrix(data_info[idx1:idx2,]))
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "white", col ="grey30")),
  NAN1= function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "white", col = NA)),
  CNV_Amp = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#d84e58", col = NA)),
  CNV_Gain=function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#f4c3c8", col = NA)),
  CNV_Del = function(x, y, w, h)
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#474dec", col = NA)),
  CNV_Loss = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#98cce8", col = NA)),
  CNV_loss = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#98cce8", col = NA)),
  
  #Funct = function(x, y, w, h) 
  #  grid.rect(x, y, w*0.6, h*0.6, gp = gpar(fill = "seagreen", col = NA)),
  Focal_Del = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#474dec", col = NA)),
  Focal_Amp = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#d84e58", col = NA)),
  #red rectangles
  Mut_Missense= function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'darkturquoise', col = NA)),
  Mut_Frameshift= function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'sienna4', col = NA)),
  Mut_Splice= function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'hotpink1', col = NA)),
  Mut_Stop_gain= function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'black', col = NA)),
  Mut_Inframe= function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = 'darkturquoise', col = NA))#,
  
)
pCNV_mut <- oncoPrint(cnvmut_matrix,
          row_order =c(1:dim(cnvmut_matrix)[1]),
          column_order = c(1:dim(cnvmut_matrix)[2]),
          alter_fun = alter_fun,
          width = unit(figure_width, "cm"),
          show_pct = F,show_row_names = T,
          show_column_names = T,
          row_names_side='left',
          row_names_gp = gpar(fontsize = size_numer_font_size,
                              col="black"),
          column_names_gp = gpar(fontsize = 5),
          height = unit(figure_height*19, "cm"),
          heatmap_legend_param = list(border=T,title="Alterations"),
          border = T)
pCNV_mut

group_matrix <- (as.matrix(data_info_ordered['final_dec',]))
group_matrix <- t(apply(group_matrix,1,function(x) as.numeric(x)))
colnames(group_matrix) <- colnames(data_info_ordered)
rownames(group_matrix) <- "group"
#colors_group_sim <- c("#FF3300","#2b8cbe")
colors_group_sim <- c("#990000","#2b8cbe")

group_sim_plot <- Heatmap(group_matrix,
              name="Group",col= colors_group_sim,
              rect_gp = gpar(col ="black", lwd =lwd_number),
              row_labels = "Group",show_column_names = F,
              cluster_rows = FALSE,cluster_columns = FALSE,
              width = unit(figure_width, "cm"),
              height = unit(figure_height, "cm"),
              row_names_side='left',row_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
              column_names_gp = gpar(fontsize = size_numer_font_size,col="black"),
              show_heatmap_legend = F,
              border = T)

# pCNV_mut%v% group_sim_plot

fig1c <- Age_plot%v%Gender_plot%v%surv_states%v%surv_time%v%ki67_plot%v%pCNV_mut%v%group_sim_plot

fig1c 

fig1bc = (fig1b)%v%fig1c
pdf('Clinical_and_muational_landscape_focal_CNV.pdf',width=15,height = 18)
draw(fig1bc, 
 ht_gap = unit(0.07, "cm"))
dev.off()


