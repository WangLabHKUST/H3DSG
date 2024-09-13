rm(list=ls())
if (!require('rstudioapi')) 
  install.packages('rstudioapi')
if (!require('survival')) 
  install.packages('survival')
if (!require('survminer')) 
  install.packages('survminer')

library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
info_plot <- read.table('Survival_comp_0521.txt',
                        sep="\t",header=T)
info_plot$final_dec <- as.factor(info_plot$final_dec)
info_plot$OS <- as.numeric(info_plot$OS)
info_plot$sensor <- as.numeric(info_plot$sensor)
library("survival")
library("survminer")
fit_4 <- survfit(Surv(OS, sensor) ~ final_dec, data = info_plot)
p1 <- ggsurvplot(fit_4,
                 pval = TRUE, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('#a6bddb','#feb24c','#990000',#'#f03b20',
                             '#2b8cbe'),
                 xlab = "Time in months",
                 legend.title="",
                 legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
pdf('Clinical_survival_comparison.pdf',width = 4.6,height = 4,
    onefile=F)
p1
dev.off()



p2 <- ggsurvplot(fit_4,
                 pval = F, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('#a6bddb','#feb24c','#990000',#'#f03b20',
                             '#2b8cbe'),
                 xlab = "Time in months",
                 legend.title="",
                 legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
)
p2

info_plot_SP <- info_plot[which(info_plot$final_dec %in% 
                                  c("SP_group1","SP_group2")),]
fit_2_sp <- survfit(Surv(OS, sensor) ~ final_dec, data = info_plot_SP)
p2_sp <- ggsurvplot(fit_2_sp,
                 pval = T, 
                 conf.int =F,
                 risk.table = F, 
                 risk.table.col = "",
                 #linetype = "strata",
                 font.tickslab = c(14),
                 surv.median.line = "hv",
                 ggtheme = theme_classic(),
                 palette = c('#f03b20',
                             '#2b8cbe'),#'#a6bddb','#feb24c',
                 xlab = "Time in months",
                 legend.title="",
                 legend.labs =c("SP_group1","SP_group2")#,
)#"H3-Medulla","H3-Pons",
p2_sp

info_plot_H3 <- info_plot[which(info_plot$final_dec %in% 
                                  c("H3-Pons","H3-Medulla")),]
fit_2_H3 <- survfit(Surv(OS, sensor) ~ final_dec, data = info_plot_H3)
p2_H3 <- ggsurvplot(fit_2_H3,
                    pval = T, 
                    conf.int =F,
                    risk.table = F, 
                    risk.table.col = "",
                    #linetype = "strata",
                    font.tickslab = c(14),
                    surv.median.line = "hv",
                    ggtheme = theme_classic(),
                    palette = c('#f03b20',
                                '#2b8cbe'),#'#a6bddb','#feb24c',
                    xlab = "Time in months",
                    legend.title="",
                    legend.labs =c("H3-Medulla","H3-Pons")#,
)#
p2_H3

info_plot_pons <- info_plot[which(info_plot$final_dec %in% 
                                  c("H3-Pons","SP_group1")),]
fit_2_Pons <- survfit(Surv(OS, sensor) ~ final_dec, data = info_plot_pons)
p2_pons <- ggsurvplot(fit_2_Pons,
                    pval = T, 
                    conf.int =F,
                    risk.table = F, 
                    risk.table.col = "",
                    #linetype = "strata",
                    font.tickslab = c(14),
                    surv.median.line = "hv",
                    ggtheme = theme_classic(),
                    palette = c('#f03b20',
                                '#feb24c'),#'#a6bddb','#feb24c',
                    xlab = "Time in months",
                    legend.title="",
                    legend.labs =c("H3-Pons","SP_roup1")#,
)#
p2_pons

info_plot_thalamus <- info_plot[which(info_plot$final_dec %in% 
                                  c("SP_group2","H3-Medulla")),]
fit_2_thalamus <- survfit(Surv(OS, sensor) ~ final_dec, data = info_plot_thalamus)
p2_thalamus <- ggsurvplot(fit_2_thalamus,
                    pval = T, 
                    conf.int =F,
                    risk.table = F, 
                    risk.table.col = "",
                    #linetype = "strata",
                    font.tickslab = c(14),
                    surv.median.line = "hv",
                    ggtheme = theme_classic(),
                    palette = c('#a6bddb',
                                '#2b8cbe'),#'#a6bddb','#feb24c',
                    xlab = "Time in months",
                    legend.title="",
                    legend.labs =c("SP_group2","H3-Pons")#,
)#
p2_thalamus



if(F){
data0520 <- readxl::read_excel("metadata_original_update0520.xlsx",sheet=1)
data0520 <- data0520[which(!(data0520$WES %in% "NA")),]

fit_2_sp_new <- survfit(Surv(data0520$`OS-1`, data0520$`sensor-1`) ~ final_dec, 
                    data = data0520)
p2_sp_new <- ggsurvplot(fit_2_sp_new,
                    pval = T, 
                    conf.int =F,
                    risk.table = F, 
                    risk.table.col = "",
                    #linetype = "strata",
                    font.tickslab = c(14),
                    surv.median.line = "hv",
                    ggtheme = theme_classic(),
                    palette = c('#f03b20',
                                '#2b8cbe'),#'#a6bddb','#feb24c',
                    xlab = "Time in months",
                    legend.title="",
                    legend.labs =c("SP_group1","SP_group2")#,
)#"H3-Medulla","H3-Pons",
p2_sp_new
}
