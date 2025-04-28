rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
### Load Train data
Train <- readRDS('Train_data_for_classifier.RDS')
source('LOOCV_related_function.R')
LOOCV_TP53 <- LOOCV_TP53_results(Train)

TP53_loocv_index_05 <- confusionMatrix(as.factor(LOOCV_TP53$Pre_Group),
                                       as.factor(LOOCV_TP53$Group))
#### index calculation
print("index for using TP53:")
TP53_loocv_index_05$byClass
TP53_loocv_index_05$overall

nb2.auprc_Train_results_TP53 <- evalmod(labels= as.numeric(LOOCV_TP53$Group),
                  scores = as.numeric(LOOCV_TP53$Prob_Group2))

#### ROC/PRC plot
print("PRC_plot for using TP53:")
plot(nb2.auprc_Train_results_TP53, "PRC")
print("ROC_plot for using TP53:")
plot(nb2.auprc_Train_results_TP53,"ROC")
#### auc calculation
print("AUC for using TP53:")
auc(nb2.auprc_Train_results_TP53)

LOOCV_NF1 <- LOOCV_NF1_results(Train)

nb2.auprc_Train_results_NF1 <- evalmod(labels= as.numeric(LOOCV_NF1$Group),
                                       scores = as.numeric(LOOCV_NF1$Prob_Group2))
NF1_loocv_index_05 <- confusionMatrix(as.factor(LOOCV_NF1$Pre_Group),
                                      as.factor(LOOCV_NF1$Group))
#### index calculation
print("index for using NF1:")
NF1_loocv_index_05$byClass
NF1_loocv_index_05$overall

#### ROC/PRC plot
print("PRC_plot for using NF1:")
plot(nb2.auprc_Train_results_NF1, "PRC")
print("ROC_plot for using NF1:")
plot(nb2.auprc_Train_results_NF1,"ROC")
#### auc calculation
print("AUC for using NF1:")
auc(nb2.auprc_Train_results_NF1)

LOOCV_age <- LOOCV_age_results(Train)
age_loocv_index_05 <- confusionMatrix(as.factor(LOOCV_age$Pre_Group),
                                      as.factor(LOOCV_age$Group))
#### index calculation
print("index for using age:")
age_loocv_index_05$byClass
age_loocv_index_05$overall
nb2.auprc_Train_results_age<- evalmod(labels= as.numeric(LOOCV_age$Group),
                                      scores = as.numeric(LOOCV_age$Prob_Group2))
### ROC/PRC plot
print("PRC_plot for using age:")
plot(nb2.auprc_Train_results_age, "PRC")
print("ROC_plot for using age:")
plot(nb2.auprc_Train_results_age,"ROC")
### auc calculation
print("AUC for using age:")
auc(nb2.auprc_Train_results_age)

LOOCV_KI67 <- LOOCV_KI67_results(Train)

nb2.auprc_Train_results_KI67<- evalmod(labels= as.numeric(LOOCV_KI67$Group),
                                       scores = as.numeric(LOOCV_KI67$Prob_Group2))#,
KI67_loocv_index_05 <- confusionMatrix(as.factor(LOOCV_KI67$Pre_Group),
                                       as.factor(LOOCV_KI67$Group))
#### index calculation
print("index for using KI67:")
KI67_loocv_index_05$byClass
KI67_loocv_index_05$overall

print("PRC_plot for using KI67:")
plot(nb2.auprc_Train_results_KI67, "PRC")
print("ROC_plot for using KI67:")
plot(nb2.auprc_Train_results_KI67,"ROC")
### auc calculation
print("AUC for using KI67:")
auc(nb2.auprc_Train_results_KI67)

LOOCV_cpg <- LOOCV_cpg_results(Train)
cpg_loocv_index_05 <- confusionMatrix(as.factor(LOOCV_cpg$Pre_Group),
                                      as.factor(LOOCV_cpg$Group))
### index calculation
print("index for using cpg:")
cpg_loocv_index_05$byClass
cpg_loocv_index_05$overall
nb2.auprc_Train_results_cpg <- evalmod(labels= as.numeric(LOOCV_cpg$Group),
                            scores = as.numeric(LOOCV_cpg$Prob_Group2))

### ROC/PRC plot
print("PRC_plot for using cpg:")
plot(nb2.auprc_Train_results_cpg, "PRC")
print("ROC_plot for using cpg:")
plot(nb2.auprc_Train_results_cpg,"ROC")
### auc calculation
print("AUC for using cpg:")
auc(nb2.auprc_Train_results_cpg)

