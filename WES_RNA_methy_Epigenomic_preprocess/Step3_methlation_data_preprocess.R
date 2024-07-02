rm(list=ls())
setwd('/Users/shanghaixia/Desktop/glioma_prognosis_method_multi-omics/Multi_omics_SP_related')
library(readxl)
clin <- read_excel('Sp_45_subtype_ID_synology_230426.xlsx',sheet=1)
baseDir<- '/Volumes/spinal_multiomics/DNA_methylation_H3_integrated_0206'
length(list.files(baseDir))
library("ChAMP")
## load data
myLoad <- champ.import(baseDir,arraytype="EPIC")
myfiler <- champ.filter(beta=myLoad$beta,
                        pd=myLoad$pd,detP=myLoad$detP,
                        beadcount=myLoad$beadcount,arraytype ="EPIC")
myNorm <- champ.norm(beta=myfiler$beta,arraytype="EPIC")
SP_36_meth <- myNorm[,which(colnames(myNorm) %in% 
                              as.character(clin$ID2...2)[which(clin$WES != "NA")])]

###Batch detect
champ.SVD(beta = myNorm,
          rgSet=NULL,
          pd=myLoad$pd,
          RGEffect=FALSE,
          PDFplot=TRUE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages/")
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

meth_Group <- rbind()
for(i in 1:dim(SP_36_meth)[2])
{
  k1 <- which(as.character(clin$ID2...2) %in% colnames(SP_36_meth)[i])
  if(length(k1)>0)
  {
    meth_Group<- rbind(meth_Group,
                       c(colnames(SP_36_meth)[i],clin$final_dec[k1]))
  }
}
colnames(meth_Group) <- c("ID","Group")
meth_Group <- as.data.frame(meth_Group)

myDMP <- champ.DMP(beta=SP_36_meth,pheno=meth_Group$Group,
                   compare.group=c("OldG","YongG"),
                   #adjPVal = 0.05,
                   arraytype = "EPIC")
DMP_info <- myDMP[["OldG_to_YongG"]]
beta_change <- DMP_info$deltaBeta

DMP_info$beta_change <-beta_change
lab <- rep('no',length(beta_change))
p_thre <- 0.05
b_thre <- 0.1#0.2
lab[intersect(which(DMP_info$adj.P.Val<p_thre),
              which((beta_change)>b_thre))] <-'hyper'
lab[intersect(which(DMP_info$adj.P.Val<p_thre),
              which((beta_change)< (-b_thre)))] <- 'hypo'
DMP_info$lab <- lab
length(which(DMP_info$lab %in% "hyper"))#13987
length(which(DMP_info$lab %in% "hypo"))#127837

