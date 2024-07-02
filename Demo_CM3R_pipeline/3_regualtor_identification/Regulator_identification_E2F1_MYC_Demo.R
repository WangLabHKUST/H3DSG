
rm(list=ls())
library(rstudioapi)
library(readxl)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
Regress_path <- gsub("3_regualtor_identification","",dirname(current_path))
Regress_path <- paste(Regress_path,"/2_constrained_Ridge_regression/Regression_results/",sep="")
Info_path <- gsub("3_regualtor_identification","",dirname(current_path))
Info_path <- paste(Info_path,"/1_multi_omics_data_preparation/Prepared_data/",sep="")
file_dir <- list.files(Regress_path)
library(stringr)
imode <- function(x)
{
  uniq <- unique(x)
  freq <- tabulate(match(x,uniq))
  mf <- max(freq)
  return(uniq[freq==mf])
}

for(i in 1:length(file_dir))
{
  
    file_cur <- file_dir[i]
    gene_cur <- gsub("_min_max_neg_cpg_cnv_tf_sigmoidtrans_Ridge_regression.csv",
                     "",file_cur)
    dd <- read.csv(paste(Regress_path,file_cur[1],sep=""),header=T)
    prepare_info <- read.csv(paste(Info_path,
                    paste0(gene_cur,'__info_min_max_score.csv'),sep="/"),header=T,row.names = 1)
   
  dd$ecludiean.distance <- round(dd$ecludiean.distance,3)
  ed_imode <- max(imode(as.numeric(dd$ecludiean.distance)))
  num_ed_imode <- length(which(dd$ecludiean.distance == ed_imode))
  used_dd <- dd[(num_ed_imode),]
  intercep_loc <- which(colnames(used_dd) == "intercept")
  rr <- as.data.frame(t(used_dd[,-c(1:2,intercep_loc)]))
  colnames(rr) <- "gene"
  
  PCC_all_cur_gene <- as.data.frame(as.matrix(prepare_info))
  #rownames(PCC_all_cur_gene) <- prepare_info[,1]
  PCC_all_cur_gene <- cor(PCC_all_cur_gene)[-1,1:2]
  PCC_all_cur_gene <- PCC_all_cur_gene[order(abs(as.numeric(PCC_all_cur_gene[,1])),
                                             decreasing = T),]
  store_PCC_name <- paste0(gene_cur,'_PCC_order.csv')
  write.csv(PCC_all_cur_gene[,1],
                       store_PCC_name,row.names = T,quote=F)
  
  
  if(length(unique(dd$ecludiean.distance))>=1 & sum(as.numeric(rr$gene))>0){
    rr <- cbind(as.matrix(rownames(rr)[order(rr$gene,decreasing = T)]),
                rr$gene[order(rr$gene,decreasing = T)])
    label <- rbind()
    for(nn in rr[,1])
    {
      if(str_detect(nn,'_cnv')==T)
      {
        label <- rbind(label,'cnv')
      }else{
        if(str_detect(nn,'cg')==T)
        {
          label <- rbind(label,'cpg')
        }else{
          label <- rbind(label,'tf')
        }
      }
    }
    rr <- cbind(rr,label)
    rr <- cbind(rr,c(1:length(label)))
    rownames(rr) <- rr[,1]
    cnv_rank <- as.numeric(rr[which(rr[,3] %in% 'cnv'),4])
    cnv_value <- as.numeric(rr[which(rr[,3] %in% 'cnv'),2])
    pp_exp_cnv<- prepare_info[,c(2,3)]
    colnames(pp_exp_cnv)<- c('exp','cnv')
    pp_exp_cnv <- as.data.frame(pp_exp_cnv)
    pp_exp_cnv$exp <- as.numeric(pp_exp_cnv$exp)
    pp_exp_cnv$cnv <- as.numeric(pp_exp_cnv$cnv)
    
    pcc_cnv <- round(cor(prepare_info[,2],
                         prepare_info[,3]),3)
    
    
   
    min_cpg_rank <- which(as.numeric(rr[which(rr[,3] %in% 'cpg'),4])
                          == min(as.numeric(rr[which(rr[,3] %in% 'cpg'),4])))
    cpg_rank <- rr[which(rr[,3] %in% 'cpg')[min_cpg_rank],4]
    cpg_value <-rr[which(rr[,3] %in% 'cpg')[min_cpg_rank],2]
    cpg_name <- rr[which(rr[,3] %in% 'cpg')[min_cpg_rank],1]
    cpg_pcc <- round(cor(as.numeric(prepare_info[,2]),
                         as.numeric(prepare_info[,which(colnames(prepare_info) %in% cpg_name)])),3)
    
    min_tf_rank <- which(as.numeric(rr[which(rr[,3] %in% 'tf'),4])
                         == min(as.numeric(rr[which(rr[,3] %in% 'tf'),4])))
    tf_rank <- rr[which(rr[,3] %in% 'tf')[min_tf_rank],4]
    tf_value <- rr[which(rr[,3] %in% 'tf')[min_tf_rank],2]
    tf_name <- rr[which(rr[,3] %in% 'tf')[min_tf_rank],1]
    tf_pcc <- round(cor(as.numeric(prepare_info[,2]),
                        as.numeric(prepare_info[,which(colnames(prepare_info) %in% tf_name)])),3)
    
    oo <- rbind(c('cnv',cnv_rank,cnv_value,'cnv',pcc_cnv),
                c('cpg',cpg_rank,cpg_value,cpg_name,cpg_pcc),
                c('tf',tf_rank,tf_value,tf_name,tf_pcc))
    k1 <- which(as.numeric(oo[,2])==min(as.numeric(oo[,2])))
    k2 <- which(as.numeric(oo[,2])==max(as.numeric(oo[,2])))
    k3 <- setdiff(c(1:length(oo[,1])),c(k1,k2))
    
    rr <- merge(rr,PCC_all_cur_gene,by="row.names",all=T)
    rr <- rr[,-c(1,dim(rr)[2])]
    rr  <- cbind(rep(gene_cur,dim(rr)[1]),rr)
    colnames(rr) <- c("Targets","Regulator","Coefficient","type","Rank","Correlation_EXP")
    write.csv(rr,paste(paste("Regulation_network_of_",gene_cur,sep=""),".csv",sep=""),row.names=F,quote=F)
    master_loc <- which(as.numeric(oo[,2])==min(as.numeric(oo[,2])))
    master <- paste(oo[master_loc,1],oo[master_loc,4],sep=":")
    print(paste(paste(paste("Master regulator of ",gene_cur,sep=""),"is: ",sep=" "),
                master))

  }

}



