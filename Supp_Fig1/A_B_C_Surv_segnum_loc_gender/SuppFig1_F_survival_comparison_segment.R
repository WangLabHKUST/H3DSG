rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
info_plot <- read.table('Survival_comp_segments_0521.txt',
                        sep="\t",header=T)
median(info_plot$segments)
median(info_plot$segments)
max(info_plot$segments)
min(info_plot$segments)
info_plot$final_dec <- as.factor(info_plot$final_dec)
info_plot$OS_1 <- as.numeric(info_plot$OS_1)
info_plot$sensor_1 <- as.numeric(info_plot$sensor_1)
median_seg <- median(info_plot$segments)
mean_seg <- mean(info_plot$segments)
info_plot$lab <- rep('Lower',length(info_plot$segments))
info_plot$lab[which(info_plot$segments>median_seg)] <- 'Higher'
library(ggplot2)
library(ggpubr)
library(ggrepel)



library("survival")
library("survminer")
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggbeeswarm)
fit_4 <- survfit(Surv(OS_1, sensor_1) ~ lab, data = info_plot)
p1 <- ggsurvplot(fit_4,
                 pval = TRUE, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('#8856a7','#9ebcda'),
                 xlab = "Time in months",
                 legend.title=""#,
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
pdf('SuppFig1_F_survival_comparison_seg_number.pdf',width = 4.6,height = 4,
    onefile=F)
p1
dev.off()
p1 <- ggsurvplot(fit_4,
                 pval = F, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('#8856a7','#9ebcda'),
                 xlab = "Time in months",
                 legend.title=""#,
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
p1
#'#f03b20',
#'#2b8cbe'
#rm(list=ls())
#library(rstudioapi)
#current_path = rstudioapi::getActiveDocumentContext()$path
#setwd(dirname(current_path))
#data <- readxl::read_excel('Survival_comp_segments.xlsx',sheet=2)
#library("survival")
#library("survminer")
#library(ggplot2)
#library(ggpubr)
#library(ggrepel)
#info_plot <- data[,c(3,4,9,12)]
#info_plot$
fit2 <- survfit(Surv(OS_1, sensor_1) ~ Sex, data = info_plot)
p2 <- ggsurvplot(fit2,
                 pval = TRUE, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('coral','cyan4'),
                 xlab = "Time in months",
                 legend.title=""#,'#8856a7','#9ebcda'
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
pdf('SuppFig1_F_survival_comparison_sex.pdf',width = 4.6,height = 4,
    onefile=F)
p2
dev.off()
p11 <- ggsurvplot(fit2,
                 pval = F, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('coral','cyan4'),
                 xlab = "Time in months",
                 legend.title=""#,'#8856a7','#9ebcda'
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)

fit3 <- survfit(Surv(OS_1, sensor_1) ~ type, data = info_plot)
p3 <- ggsurvplot(fit3,
                 pval = TRUE, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c("lightblue4","lightsteelblue"),
                 xlab = "Time in months",
                 legend.title=""#,'#8856a7','#9ebcda'
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
pdf('SuppFig1_F_survival_comparison_location_sim.pdf',width = 4.6,height = 4,
    onefile=F)
p3
dev.off()
p22 <- ggsurvplot(fit3,
                 pval = F, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c("lightblue4","lightsteelblue"),
                 xlab = "Time in months",
                 legend.title=""#,'#8856a7','#9ebcda'
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
p22


