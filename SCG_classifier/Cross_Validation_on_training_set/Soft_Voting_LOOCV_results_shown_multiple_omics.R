rm(list=ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
### Load Train data
Train <- readRDS('Train_data_for_classifier.RDS')
source('LOOCV_related_function.R')
LOOCV_TP53 <- LOOCV_TP53_results(Train)
LOOCV_NF1 <- LOOCV_NF1_results(Train)
LOOCV_age <- LOOCV_age_results(Train)
LOOCV_KI67 <- LOOCV_KI67_results(Train)
LOOCV_cpg <- LOOCV_cpg_results(Train)

### Voting by methy, age, KI67
MAK_voting <- cbind(LOOCV_age,cbind(LOOCV_cpg,LOOCV_KI67))
MAK_voting$votingProb_Group1 <- rep(0,length(MAK_voting$ID))
MAK_voting$votingProb_Group2 <- rep(0,length(MAK_voting$ID))
MAK_voting$voting_Group <- rep(0,length(MAK_voting$ID))
for(i in 1:dim(MAK_voting)[1])
{
  MAK_voting$votingProb_Group1[i] <- mean(c(as.numeric(MAK_voting[i,3]),
                                            as.numeric(MAK_voting[i,8]),
                                            as.numeric(MAK_voting[i,13])))
  MAK_voting$votingProb_Group2[i] <- mean(c(as.numeric(MAK_voting[i,4]),
                                            as.numeric(MAK_voting[i,9]),
                                            as.numeric(MAK_voting[i,14])))
  if(MAK_voting$votingProb_Group1[i]>MAK_voting$votingProb_Group2[i])
  {
    MAK_voting$voting_Group[i] <- 1
  }else{
    MAK_voting$voting_Group[i] <- 2
  }
}
cpg_age_ki67_loocv_index_05 <- confusionMatrix(as.factor(MAK_voting$voting_Group),
                                               as.factor(MAK_voting$Group))

print("Index cal for Voting by using cpg, age, ki67:")
cpg_age_ki67_loocv_index_05$overall
cpg_age_ki67_loocv_index_05$byClass


### Voting by methy, age, TP53
MAT_voting <- cbind(LOOCV_age,cbind(LOOCV_cpg,LOOCV_TP53))
MAT_voting$votingProb_Group1 <- rep(0,length(MAT_voting$ID))
MAT_voting$votingProb_Group2 <- rep(0,length(MAT_voting$ID))
MAT_voting$voting_Group <- rep(0,length(MAT_voting$ID))
for(i in 1:dim(MAT_voting)[1])
{
  MAT_voting$votingProb_Group1[i] <- mean(c(as.numeric(MAT_voting[i,3]),
                                            as.numeric(MAT_voting[i,8]),
                                            as.numeric(MAT_voting[i,13])))
  MAT_voting$votingProb_Group2[i] <- mean(c(as.numeric(MAT_voting[i,4]),
                                            as.numeric(MAT_voting[i,9]),
                                            as.numeric(MAT_voting[i,14])))
  if(MAT_voting$votingProb_Group1[i]>MAT_voting$votingProb_Group2[i])
  {
    MAT_voting$voting_Group[i] <- 1
  }else{
    MAT_voting$voting_Group[i] <- 2
  }
}
cpg_age_TP53_loocv_index_05 <- confusionMatrix(as.factor(MAT_voting$voting_Group),
                                as.factor(MAT_voting$Group))
print("Index cal for Voting by using cpg, age, TP53:")
cpg_age_TP53_loocv_index_05$overall
cpg_age_TP53_loocv_index_05$byClass



### Voting by methy, age, NF1
MAN_voting <- cbind(LOOCV_age,cbind(LOOCV_cpg,LOOCV_NF1))
MAN_voting$votingProb_Group1 <- rep(0,length(MAN_voting$ID))
MAN_voting$votingProb_Group2 <- rep(0,length(MAN_voting$ID))
MAN_voting$voting_Group <- rep(0,length(MAN_voting$ID))
for(i in 1:dim(MAN_voting)[1])
{
  MAN_voting$votingProb_Group1[i] <- mean(c(as.numeric(MAN_voting[i,3]),
                                            as.numeric(MAN_voting[i,8]),
                                            as.numeric(MAN_voting[i,13])))
  MAN_voting$votingProb_Group2[i] <- mean(c(as.numeric(MAN_voting[i,4]),
                                            as.numeric(MAN_voting[i,9]),
                                            as.numeric(MAN_voting[i,14])))
  if(MAN_voting$votingProb_Group1[i]>MAN_voting$votingProb_Group2[i])
  {
    MAN_voting$voting_Group[i] <- 1
  }else{
    MAN_voting$voting_Group[i] <- 2
  }
}
cpg_age_NF1_loocv_index_05 <- confusionMatrix(as.factor(MAN_voting$voting_Group),
                                              as.factor(MAN_voting$Group))
print("Index cal for Voting by using cpg, age, NF1:")
cpg_age_NF1_loocv_index_05$overall
cpg_age_NF1_loocv_index_05$byClass

### Voting by methy, KI67, TP53
MKP_voting <- cbind(LOOCV_KI67,cbind(LOOCV_cpg,LOOCV_TP53))
MKP_voting$votingProb_Group1 <- rep(0,length(MKP_voting$ID))
MKP_voting$votingProb_Group2 <- rep(0,length(MKP_voting$ID))
MKP_voting$voting_Group <- rep(0,length(MKP_voting$ID))
for(i in 1:dim(MKP_voting)[1])
{
  MKP_voting$votingProb_Group1[i] <- mean(c(as.numeric(MKP_voting[i,3]),
                                            as.numeric(MKP_voting[i,8]),
                                            as.numeric(MKP_voting[i,13])))
  MKP_voting$votingProb_Group2[i] <- mean(c(as.numeric(MKP_voting[i,4]),
                                            as.numeric(MKP_voting[i,9]),
                                            as.numeric(MKP_voting[i,14])))
  if(MKP_voting$votingProb_Group1[i]>MKP_voting$votingProb_Group2[i])
  {
    MKP_voting$voting_Group[i] <- 1
  }else{
    MKP_voting$voting_Group[i] <- 2
  }
}
cpg_KI67_TP53_loocv_index_05 <- confusionMatrix(as.factor(MKP_voting$voting_Group),
                                                as.factor(MKP_voting$Group))
print("Index cal for Voting by using cpg, ki67, TP53:")
cpg_KI67_TP53_loocv_index_05$overall
cpg_KI67_TP53_loocv_index_05$byClass

### Voting by methy, KI67, NF1
MKN_voting <- cbind(LOOCV_KI67,cbind(LOOCV_cpg,LOOCV_NF1))
MKN_voting$votingProb_Group1 <- rep(0,length(MKN_voting$ID))
MKN_voting$votingProb_Group2 <- rep(0,length(MKN_voting$ID))
MKN_voting$voting_Group <- rep(0,length(MKN_voting$ID))
for(i in 1:dim(MKN_voting)[1])
{
  MKN_voting$votingProb_Group1[i] <- mean(c(as.numeric(MKN_voting[i,3]),
                                            as.numeric(MKN_voting[i,8]),
                                            as.numeric(MKN_voting[i,13])))
  MKN_voting$votingProb_Group2[i] <- mean(c(as.numeric(MKN_voting[i,4]),
                                            as.numeric(MKN_voting[i,9]),
                                            as.numeric(MKN_voting[i,14])))
  if(MKN_voting$votingProb_Group1[i]>MKN_voting$votingProb_Group2[i])
  {
    MKN_voting$voting_Group[i] <- 1
  }else{
    MKN_voting$voting_Group[i] <- 2
  }
}
cpg_KI67_NF1_loocv_index_05 <- confusionMatrix(as.factor(MKN_voting$voting_Group),
                                               as.factor(MKN_voting$Group))
print("Index cal for Voting by using cpg, ki67, NF1:")
cpg_KI67_NF1_loocv_index_05$overall
cpg_KI67_NF1_loocv_index_05$byClass

### Voting by methy, TP53, NF1
MTN_voting <- cbind(LOOCV_TP53,cbind(LOOCV_cpg,LOOCV_NF1))
MTN_voting$votingProb_Group1 <- rep(0,length(MTN_voting$ID))
MTN_voting$votingProb_Group2 <- rep(0,length(MTN_voting$ID))
MTN_voting$voting_Group <- rep(0,length(MTN_voting$ID))
for(i in 1:dim(MTN_voting)[1])
{
  MTN_voting$votingProb_Group1[i] <- mean(c(as.numeric(MTN_voting[i,3]),
                                            as.numeric(MTN_voting[i,8]),
                                            as.numeric(MTN_voting[i,13])))
  MTN_voting$votingProb_Group2[i] <- mean(c(as.numeric(MTN_voting[i,4]),
                                            as.numeric(MTN_voting[i,9]),
                                            as.numeric(MTN_voting[i,14])))
  if(MTN_voting$votingProb_Group1[i]>MTN_voting$votingProb_Group2[i])
  {
    MTN_voting$voting_Group[i] <- 1
  }else{
    MTN_voting$voting_Group[i] <- 2
  }
}
cpg_TP53_NF1_loocv_index_05 <- confusionMatrix(as.factor(MTN_voting$voting_Group),
                                               as.factor(MTN_voting$Group))
print("Index cal for Voting by using cpg,TP53, NF1:")
cpg_TP53_NF1_loocv_index_05$overall
cpg_TP53_NF1_loocv_index_05$byClass


### Voting by ki67, age, NF1
KAN_voting <- cbind(LOOCV_KI67,cbind(LOOCV_age,LOOCV_NF1))
KAN_voting$votingProb_Group1 <- rep(0,length(KAN_voting$ID))
KAN_voting$votingProb_Group2 <- rep(0,length(KAN_voting$ID))
KAN_voting$voting_Group <- rep(0,length(KAN_voting$ID))
for(i in 1:dim(KAN_voting)[1])
{
  KAN_voting$votingProb_Group1[i] <- mean(c(as.numeric(KAN_voting[i,3]),
                                            as.numeric(KAN_voting[i,8]),
                                            as.numeric(KAN_voting[i,13])))
  KAN_voting$votingProb_Group2[i] <- mean(c(as.numeric(KAN_voting[i,4]),
                                            as.numeric(KAN_voting[i,9]),
                                            as.numeric(KAN_voting[i,14])))
  if(KAN_voting$votingProb_Group1[i]>KAN_voting$votingProb_Group2[i])
  {
    KAN_voting$voting_Group[i] <- 1
  }else{
    KAN_voting$voting_Group[i] <- 2
  }
}
age_KI67_NF1_loocv_index_05 <- confusionMatrix(as.factor(KAN_voting$voting_Group),
                                               as.factor(KAN_voting$Group))
print("Index cal for Voting by using age,ki67, NF1:")
age_KI67_NF1_loocv_index_05$overall
age_KI67_NF1_loocv_index_05$byClass

### Voting by ki67, age, TP53
KAT_voting <- cbind(LOOCV_KI67,cbind(LOOCV_age,LOOCV_TP53))
KAT_voting$votingProb_Group1 <- rep(0,length(KAT_voting$ID))
KAT_voting$votingProb_Group2 <- rep(0,length(KAT_voting$ID))
KAT_voting$voting_Group <- rep(0,length(KAT_voting$ID))
for(i in 1:dim(KAT_voting)[1])
{
  KAT_voting$votingProb_Group1[i] <- mean(c(as.numeric(KAT_voting[i,3]),
                                            as.numeric(KAT_voting[i,8]),
                                            as.numeric(KAT_voting[i,13])))
  KAT_voting$votingProb_Group2[i] <- mean(c(as.numeric(KAT_voting[i,4]),
                                            as.numeric(KAT_voting[i,9]),
                                            as.numeric(KAT_voting[i,14])))
  if(KAT_voting$votingProb_Group1[i]>KAT_voting$votingProb_Group2[i])
  {
    KAT_voting$voting_Group[i] <- 1
  }else{
    KAT_voting$voting_Group[i] <- 2
  }
}
age_KI67_TP53_loocv_index_05 <- confusionMatrix(as.factor(KAT_voting$voting_Group),
                                                as.factor(KAT_voting$Group))
print("Index cal for Voting by using age, ki67, TP53:")
age_KI67_TP53_loocv_index_05$overall
age_KI67_TP53_loocv_index_05$byClass

### Voting by ki67, NF1, TP53
KNT_voting <- cbind(LOOCV_KI67,cbind(LOOCV_NF1,LOOCV_TP53))
KNT_voting$votingProb_Group1 <- rep(0,length(KNT_voting$ID))
KNT_voting$votingProb_Group2 <- rep(0,length(KNT_voting$ID))
KNT_voting$voting_Group <- rep(0,length(KNT_voting$ID))
for(i in 1:dim(KAT_voting)[1])
{
  KNT_voting$votingProb_Group1[i] <- mean(c(as.numeric(KNT_voting[i,3]),
                                            as.numeric(KNT_voting[i,8]),
                                            as.numeric(KNT_voting[i,13])))
  KNT_voting$votingProb_Group2[i] <- mean(c(as.numeric(KNT_voting[i,4]),
                                            as.numeric(KNT_voting[i,9]),
                                            as.numeric(KNT_voting[i,14])))
  if(KNT_voting$votingProb_Group1[i]>KNT_voting$votingProb_Group2[i])
  {
    KNT_voting$voting_Group[i] <- 1
  }else{
    KNT_voting$voting_Group[i] <- 2
  }
}
NF1_KI67_TP53_loocv_index_05 <- confusionMatrix(as.factor(KNT_voting$voting_Group),
                                                as.factor(KNT_voting$Group))
print("Index cal for Voting by using NF1, ki67, TP53:")
NF1_KI67_TP53_loocv_index_05$overall
NF1_KI67_TP53_loocv_index_05$byClass

### Voting by methy, age, TP53,KI67, NF1
MAKTN_voting <- cbind(cbind(cbind(LOOCV_age,LOOCV_KI67),cbind(LOOCV_cpg,LOOCV_TP53)),
                      LOOCV_NF1)
MAKTN_voting$votingProb_Group1 <- rep(0,length(MAKTN_voting$ID))
MAKTN_voting$votingProb_Group2 <- rep(0,length(MAKTN_voting$ID))
MAKTN_voting$voting_Group <- rep(0,length(MAKTN_voting$ID))
for(i in 1:dim(MAKTN_voting)[1])
{
  MAKTN_voting$votingProb_Group1[i] <- mean(c(as.numeric(MAKTN_voting[i,3]),
                                              as.numeric(MAKTN_voting[i,8]),
                                              as.numeric(MAKTN_voting[i,13]),
                                              as.numeric(MAKTN_voting[i,18]),
                                              as.numeric(MAKTN_voting[i,23])
  ))
  MAKTN_voting$votingProb_Group2[i] <- mean(c(as.numeric(MAKTN_voting[i,4]),
                                              as.numeric(MAKTN_voting[i,9]),
                                              as.numeric(MAKTN_voting[i,14]),
                                              as.numeric(MAKTN_voting[i,19]),
                                              as.numeric(MAKTN_voting[i,24])
  ))
  if(MAKTN_voting$votingProb_Group1[i]>MAKTN_voting$votingProb_Group2[i])
  {
    MAKTN_voting$voting_Group[i] <- 1
  }else{
    MAKTN_voting$voting_Group[i] <- 2
  }
}
cpg_age_TP53_ki67_NF1_loocv_index_05 <- confusionMatrix(as.factor(MAKTN_voting$voting_Group),
                                        as.factor(MAKTN_voting$Group))
print("Index cal for Voting by using age, cpg, NF1, ki67, TP53:")
cpg_age_TP53_ki67_NF1_loocv_index_05$overall
cpg_age_TP53_ki67_NF1_loocv_index_05$byClass



