rm(list=ls())

# run from RStudio, or with Rscript from anywhere
source('required_packages.R')
library(ComplexHeatmap)
library(circlize)
figure_width <- 6
figure_height <- 0.28
size_numer_font_size <- 9
lwd_number <- 0.2

data_path = '.'
out_path = '.'
# TP53 = Mut_Missense was filled in for P09 and P18. Both cells had been empty
#   since the 2026-01 version, but the hg38 columns of Table S2 list SP09 TP53 R209Q
#   and SP18 TP53 R209Q + S202A, and Table S1 and every other source call them Mut.
data_info = read.csv(file.path(data_path, 'Data_clinical_info.csv'))
# data_info = data.frame(readxl::read_excel(file.path(data_path, 'Data_clinical_info - 0615_archive.xlsx'),sheet=1))
rownames(data_info) <- data_info[,1]
data_info <- data_info[,-1]
# Columns (samples) are sorted P01..P36 for every block. Fig 1b (tumour location,
#   mutation burden bars) used to take the csv order (P09, P13, P07, ...) while Fig 1c
#   used the sorted order; stacking them makes ComplexHeatmap align by position rather
#   than by sample name, so the column labels of the top half and the bottom half were
#   not the same sample (only 1 of 36 columns agreed, and draw() warned "Column names
#   of heatmap N are not consistent with the main heatmap").
data_info <- data_info[, order(colnames(data_info))]
data_info[1:3,1:3]
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



mut_load_matrix <-(as.matrix(data_info['Count',]))
mut_load_matrix <- t(apply(mut_load_matrix,1,function(x) as.numeric(x)))
colnames(mut_load_matrix) <- colnames(data_info)
rownames(mut_load_matrix) <- "Mut_load"

col_ha <- HeatmapAnnotation(bar1 = anno_barplot(as.numeric(mut_load_matrix), 
                                                gp = gpar(fill = "grey"),border = F),
                            show_annotation_name = F)

loc_name1 = 'C1'; loc_name2 = 'L1'#'S5'
loc_idx1 = which(rownames(data_info) == loc_name1);loc_idx2 = which(rownames(data_info) == loc_name2);
location_matrix <- (as.matrix(data_info[loc_idx1:loc_idx2,]))
location_matrix <- apply(as.matrix(location_matrix),2,
                    function(x) as.numeric(x))


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
# pdf(file.path(out_path, 'fig1b.pdf'),width=15,height = 18)
# draw(fig1b, 
#  ht_gap = unit(0.07, "cm"))
# dev.off()


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
idx1 = which(rownames(data_info_ordered)=='H3-3A')
idx2 = which(rownames(data_info_ordered) == 'FGFR1')
cnvmut_matrix <- (as.matrix(data_info_ordered[idx1:idx2,]))
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
## replace "NA" with NA
cnvmut_matrix[cnvmut_matrix == "NA"] <- NA

# Where one sample carries two mutation types in the same gene (from Table S2, see
#   multi_hit_from_TableS2.R), the second type is drawn as a right triangle over the
#   lower left half of the cell, on top of the full-cell rectangle. The order of the
#   alter_fun list is the drawing order, so the triangles come last.
mut_fill <- c(Mut_Missense = 'darkturquoise', Mut_Frameshift = 'sienna4', Mut_Splice = 'hotpink1',
              Mut_Stop_gain = 'black', Mut_Inframe = 'darkturquoise')
tri_lower_left <- function(fill) function(x, y, w, h)
  grid.polygon(unit.c(x - w*0.45, x - w*0.45, x + w*0.45),
               unit.c(y + h*0.45, y - h*0.45, y - h*0.45),
               gp = gpar(fill = fill, col = NA))
alter_fun_tri <- c(alter_fun, setNames(lapply(mut_fill, tri_lower_left), paste0(names(mut_fill), '_2nd')))

multi_hit <- read.csv(file.path(data_path, 'landscape_multi_hit.csv'), check.names = FALSE)
multi_hit <- multi_hit[!is.na(multi_hit$triangle_type), ]
cnvmut_matrix_tri <- cnvmut_matrix
for (i in seq_len(nrow(multi_hit))) {
  g <- multi_hit$gene[i]; p <- multi_hit$sample[i]
  stopifnot(identical(cnvmut_matrix[g, p], multi_hit$drawn_type[i]))    # the type already drawn must match
  cnvmut_matrix_tri[g, p] <- paste0(cnvmut_matrix[g, p], ';', multi_hit$triangle_type[i], '_2nd')
}

make_oncoprint <- function(mat, af) {
  types <- unique(unlist(strsplit(na.omit(as.vector(mat)), ';')))
  legend_types <- intersect(names(af), setdiff(types, grep('_2nd$', types, value = TRUE)))   # _2nd is not listed in the legend
  oncoPrint(mat,
          row_order =c(1:dim(mat)[1]),
          column_order = c(1:dim(mat)[2]),
          alter_fun = af, alter_fun_is_vectorized = FALSE,
          width = unit(figure_width, "cm"),
          show_pct = F,show_row_names = T,
          show_column_names = T,
          row_names_side='left',
          row_names_gp = gpar(fontsize = size_numer_font_size,
                              col="black"),
          column_names_gp = gpar(fontsize = 5),
          height = unit(figure_height*19, "cm"),
          heatmap_legend_param = list(border=T,title="Alterations", at = legend_types, labels = legend_types),
          border = T)
}
pCNV_mut_tri <- make_oncoprint(cnvmut_matrix_tri, alter_fun_tri)

# percent (count) of mutated patients per gene, out of all samples; a patient with two
# mutation types in one cell is counted once
n_mut_patients <- rowSums(!is.na(cnvmut_matrix) & grepl('Mut_', cnvmut_matrix))
mut_pct_labels <- setNames(sprintf('%.2f%%(%d)', 100 * n_mut_patients / ncol(cnvmut_matrix), n_mut_patients),
                           rownames(cnvmut_matrix))

group_matrix <- (as.matrix(data_info_ordered['final_dec',]))
group_matrix <- t(apply(group_matrix,1,function(x) as.numeric(x)))
colnames(group_matrix) <- colnames(data_info_ordered)
rownames(group_matrix) <- "group"
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

fig1bc_tri = (fig1b)%v%(group_sim_plot%v%Age_plot%v%Gender_plot%v%surv_states%v%surv_time%v%ki67_plot%v%pCNV_mut_tri)

# Only the final version is written: thick vertical lines between the samples (columns) of
# every heatmap body (all other line widths unchanged), the second-mutation triangles and
# the per-gene frequency labels. The earlier plain / no-triangle versions were removed on
# 2026-11-05.
col_sep_lwd <- 1.2          # column separator width (ordinary cell border is lwd_number = 0.2)
# row_text: named vector (row name -> text), drawn to the right of the heatmap body whose
# row names match it exactly; added after the plot is drawn, so it changes no layout
draw_fig1bc <- function(file, col_sep = NULL, fig = fig1bc_tri, row_text = NULL) {
  pdf(file, width = 15, height = 18)
  if (!is.null(row_text)) {                          # make room for the text: shift the legend right, heatmap size unchanged
    old_pad <- ht_opt$HEATMAP_LEGEND_PADDING
    ht_opt$HEATMAP_LEGEND_PADDING <- unit(2.6, "cm")
    on.exit(ht_opt$HEATMAP_LEGEND_PADDING <- old_pad, add = TRUE)
  }
  ht <- draw(fig, ht_gap = unit(0.07, "cm"))
  if (!is.null(row_text)) {
    target <- Filter(function(nm) inherits(ht@ht_list[[nm]], 'Heatmap') &&
                       identical(rownames(ht@ht_list[[nm]]@matrix), names(row_text)), names(ht@ht_list))
    stopifnot(length(target) == 1)
    n <- length(row_text)
    decorate_heatmap_body(target, grid.text(unname(row_text), x = unit(1, "npc") + unit(1.5, "mm"),
                                            y = unit(1 - (seq_len(n) - 0.5) / n, "npc"), just = "left",
                                            gp = gpar(fontsize = size_numer_font_size, col = "black")))
  }
  # every block must have identical column names (aligned by sample)
  cn <- lapply(Filter(function(h) inherits(h, 'Heatmap'), ht@ht_list), function(h) colnames(h@matrix))
  cn <- Filter(Negate(is.null), cn)
  stopifnot(all(vapply(cn, identical, logical(1), cn[[1]])))
  if (!is.null(col_sep)) {
    for (nm in names(ht@ht_list)) {
      if (!inherits(ht@ht_list[[nm]], 'Heatmap')) next
      n <- ncol(ht@ht_list[[nm]]@matrix)
      decorate_heatmap_body(nm, {
        x <- (seq_len(n - 1)) / n
        grid.segments(x0 = x, x1 = x, y0 = 0, y1 = 1, gp = gpar(col = "black", lwd = col_sep))
      })
    }
  }
  dev.off()
}
# thick separators plus the second-mutation triangles (replaces pasting triangles in PowerPoint)
draw_fig1bc(file.path(out_path, 'Clinical_and_muational_landscape_focal_CNV_thick_col_sep_multi_hit.pdf'),
            col_sep = col_sep_lwd, fig = fig1bc_tri, row_text = mut_pct_labels)

