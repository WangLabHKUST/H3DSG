
library(naivebayes)
library(caret)
library(precrec)
LOOCV_TP53_results <- function(Train)
{
  lab_loc <- which(colnames(Train) %in% "Group")
  feature_loc_TP53 <- which(colnames(Train) %in%c("TP53"))
  LOOCV_TP53 <- rbind()
  for(i in 1:dim(Train)[1])
  {
    model_fina_TP53 <- naive_bayes(Group~., data=Train[-i,c(feature_loc_TP53,lab_loc)])
    train_TP53_test <- as.data.frame(Train[i,feature_loc_TP53])
    colnames(train_TP53_test) <- "TP53"
    Test_results_TP53 <- predict(model_fina_TP53, train_TP53_test,type="prob")
    Test_results_TP53_group <- predict(model_fina_TP53, train_TP53_test)
    LOOCV_TP53 <- rbind(LOOCV_TP53,c(rownames(Train)[i],
                                     Train$Group[i],Test_results_TP53,Test_results_TP53_group))
  }
  colnames(LOOCV_TP53) <- c("ID","Group","Prob_Group1","Prob_Group2","Pre_Group")
  LOOCV_TP53 <- as.data.frame(LOOCV_TP53)
  return(LOOCV_TP53)
}

LOOCV_NF1_results <- function(Train)
{
  lab_loc <- which(colnames(Train) %in% "Group")
feature_loc_NF1 <- which(colnames(Train) %in%c("NF1"))
LOOCV_NF1 <- rbind()
for(i in 1:dim(Train)[1])
{
  model_fina_NF1 <- naive_bayes(Group~., data=Train[-i,c(feature_loc_NF1,lab_loc)],
                                laplace = T)
  train_NF1_test <- as.data.frame(Train[i,feature_loc_NF1])
  colnames(train_NF1_test) <- "NF1"
  train_results_NF1 <- predict(model_fina_NF1, train_NF1_test,type="prob")
  train_results_NF1_group <- predict(model_fina_NF1, train_NF1_test)
  LOOCV_NF1 <- rbind(LOOCV_NF1,c(rownames(Train)[i],
                                 Train$Group[i],train_results_NF1,train_results_NF1_group))
}
colnames(LOOCV_NF1) <- c("ID","Group","Prob_Group1","Prob_Group2","Pre_Group")
LOOCV_NF1 <- as.data.frame(LOOCV_NF1)
return(LOOCV_NF1)
}

LOOCV_age_results <- function(Train)
{
  lab_loc <- which(colnames(Train) %in% "Group")
  feature_loc_age <- which(colnames(Train) %in%c("Age"))
  LOOCV_age <- rbind()
  for(i in 1:dim(Train)[1])
  {
    model_fina_age <- naive_bayes(Group~., data=Train[-i,c(feature_loc_age,lab_loc)])
    train_age_test <- as.data.frame(Train[i,feature_loc_age])
    colnames(train_age_test) <- "Age"
    train_results_age <- predict(model_fina_age, train_age_test,type="prob")
    train_results_group <- predict(model_fina_age, train_age_test)
    LOOCV_age <- rbind(LOOCV_age,c(rownames(Train)[i],
                                   Train$Group[i],train_results_age,train_results_group))
  }
  colnames(LOOCV_age) <- c("ID","Group","Prob_Group1","Prob_Group2","Pre_Group")
  LOOCV_age <- as.data.frame(LOOCV_age)
  return(LOOCV_age)
}

LOOCV_KI67_results <- function(Train)
{
  lab_loc <- which(colnames(Train) %in% "Group")
  feature_loc_KI67<- which(colnames(Train) %in%c("KI67"))
  LOOCV_KI67 <- rbind()
  for(i in 1:dim(Train)[1])
  {
    model_fina_KI67 <- naive_bayes(Group~., data=Train[-i,c(feature_loc_KI67,lab_loc)],
                                   usekernel = T
    )
    train_KI67_test <- as.data.frame(Train[i,feature_loc_KI67])
    colnames(train_KI67_test) <- "KI67"
    train_results_KI67 <- predict(model_fina_KI67, train_KI67_test,type="prob")
    train_results_KI67_group <- predict(model_fina_KI67, train_KI67_test)
    LOOCV_KI67 <- rbind(LOOCV_KI67,c(rownames(Train)[i],
                                     Train$Group[i],train_results_KI67,train_results_KI67_group))
  }
  colnames(LOOCV_KI67) <- c("ID","Group","Prob_Group1","Prob_Group2","Pre_Group")
  LOOCV_KI67 <- as.data.frame(LOOCV_KI67)
  return(LOOCV_KI67)
}

LOOCV_cpg_results <- function(Train)
{
  lab_loc <- which(colnames(Train) %in% "Group")
  feature_loc_cpg <- which(colnames(Train) %in%c("cg16549043"))
  LOOCV_cpg <- rbind()
  for(i in 1:dim(Train)[1])
  {
    model_fina_cpg <- naive_bayes(Group~., data=Train[-i,c(feature_loc_cpg,lab_loc)],
                                  usekernel = T,kernel = "biweight"
    )
    train_cpg_test <- as.data.frame(Train[i,feature_loc_cpg])
    colnames(train_cpg_test) <- "cg16549043"
    train_results_cpg <- predict(model_fina_cpg, train_cpg_test,type="prob")
    train_results_cpg_group <- predict(model_fina_cpg, train_cpg_test)
    LOOCV_cpg <- rbind(LOOCV_cpg,c(rownames(Train)[i],
                                   Train$Group[i],train_results_cpg,train_results_cpg_group))
  }
  colnames(LOOCV_cpg) <- c("ID","Group","Prob_Group1","Prob_Group2","Pre_Group")
  LOOCV_cpg <- as.data.frame(LOOCV_cpg)
}



