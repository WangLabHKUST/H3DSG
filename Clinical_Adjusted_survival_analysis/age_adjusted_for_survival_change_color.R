
rm(list=ls())
### package load
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
#### data import
###age read
#age_dir <-'/Users/shx/Downloads/H3DSG-main/Clinical_Age_comp'
age_info <- read.table('Age_comp.txt',sep="\t",header=T)
rownames(age_info) <- age_info$M_ID
rownames(age_info)[which(rownames(age_info) %in% "SPH3_7")] <- "SPH3_07"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_6")] <- "SPH3_06"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_4")] <- "SPH3_04"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_3")] <- "SPH3_03"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_5")] <- "SPH3_05"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_9")] <- "SPH3_09"
rownames(age_info)[which(rownames(age_info) %in% "SPH3_8")] <- "SPH3_08"
## survival read

surv_info <- read.table('Survival_comp_0521.txt',sep="\t",header=T)
rownames(surv_info) <- surv_info$M_ID
merge_age_surv <- merge(age_info,
                        surv_info,by="row.names",all=T)
rownames(merge_age_surv) <- merge_age_surv$Row.names
merge_age_surv <- merge_age_surv[,-1]
#
### survival analysis
merge_age_surv_only <- merge_age_surv[which(merge_age_surv$final_dec.x
                                    %in% c("SP_group2","SP_group1")),]
library("survival")
library("survminer")
#fit_4 <- survfit(Surv(OS, sensor) ~ final_dec.x,
#                 data = merge_age_surv_only)
#p1 <- ggsurvplot(fit_4,pval = TRUE, conf.int =F,
#                 risk.table = F, risk.table.col = "",
                 #linetype = "strata",
#                 font.tickslab = c(14),surv.median.line = "hv",
#                 ggtheme = theme_classic(),
                 #palette = c('#a6bddb','#feb24c','#990000',#'#f03b20',
                 #            '#2b8cbe'),
 #                xlab = "Time in months",legend.title=""#,
                 #legend.labs =c("H3-Medulla","H3-Pons","SP_group1","SP_group2")#,
#)
#p1
##https://robindenz1.github.io/adjustedCurves/
library(adjustedCurves)
set.seed(31)
##
SPH3group <- merge_age_surv[which(merge_age_surv$final_dec.x
                %in% c("SP_group2","SP_group1")),]
####estimate a cox-regression for the outcome
SPH3group$final_dec.x <- as.factor(SPH3group$final_dec.x)
cox_mod <- coxph(Surv(OS, sensor)~Age+final_dec.x,
                 data=SPH3group,x=T)
#### use it to estimate adjusted survival curves
adjsurv <- adjustedsurv(data=SPH3group,variable="final_dec.x",
                        ev_time="OS",event="sensor",
                        method="direct",bootstrap=T,
            n_boot = 1000,outcome_model=cox_mod,conf_int=TRUE)

p11 <- plot(adjsurv, conf_int=F, risk_table=T, risk_table_stratify=T,
     risk_table_digits=0, x_n_breaks=10,
     ggtheme = theme_classic(),
     censoring_ind="points",
     median_surv_lines=T,#pval = TRUE,
    custom_colors = c('#990000',#'#f03b20',#
                 '#2b8cbe'),#'#a6bddb','#feb24c',
    censoring_ind_size=2)
pdf('adjusted_survival_SP.pdf',width=6,height = 6)
p11
dev.off()

adj_test <- adjusted_curve_test(adjsurv, from=0, to=50)
summary(adj_test)

