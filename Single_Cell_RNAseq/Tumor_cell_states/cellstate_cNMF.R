library(ggplot2)
library(dendextend)
library(purrr)
library(dplyr)

setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/scRNA/rerun_cNMF/filteredmalignant/')
# scRNAseq: SP04,SP05,SP11,SP13,SP15,SP22,SP24,SP41; SP20, SP21;


custom_magma <- c(colorRampPalette(c("white", rev(viridis::magma(323, begin = 0.15))[1]))(10), rev(viridis::magma(323, begin = 0.18)))

#correlate qc metric and program usage
qcmetric_programweight_corr <-function(sid = 'SP04'){
  load(paste0('../../rerun_preprocessing/',sid,'/seurat.Rda'))
  qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
  cc_cols <- c("S.Score", "G2M.Score")
  
  mqc = seurat@meta.data[,c(cc_cols, qc_cols)]
  
  rm(seurat)
  
  fl = list.files(path = paste0('../filteredmalignant/',sid),pattern = 'dt_0_02.consensus.txt',full.names = T)
  fl=fl[grep('usages',fl)]
  mpw = read.delim(file = fl[1], row.names = 1, check.names = F)
  names(mpw) = paste0(sid,'_P',names(mpw))
  ovlapvells = intersect(rownames(mqc), rownames(mpw))
  
  corm <-cor(mqc[ovlapvells,],mpw[ovlapvells,]) %>% t() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var = "Program") %>% 
    tibble::add_column(.before = 1, "Sample" = sid) 
  return(corm)
}

sc_samples = c('SP04','SP05','SP11','SP13','SP15','SP22','SP24','SP41','SP20T5','SP21T5')
for (i in 1:length(sc_samples)){
  tmp = qcmetric_programweight_corr(sid = sc_samples[i])
  if (i == 1){
    cor_qcpw <-tmp
  }else{
    cor_qcpw <-rbind(cor_qcpw ,tmp)
  }
}

#
extract_program <-function(sid = 'SP04'){
  gs_p689 = read.delim(paste0(sid,'.programs.genes100.csv'), sep = ",", row.names = 1, check.names = F)
  names(gs_p689) = paste0(sid,'_P', names(gs_p689))
  selK = ncol(gs_p689)
  pw_p689 = read.delim(paste0(sid,'/',sid,'.usages.k_',selK,'.dt_0_02.consensus.txt'))
  
  pw_p689$group = apply(pw_p689[,-1],1,which.max)
  pw_p689$group =paste0(sid,'_P',pw_p689$group )
  pw.ncl = as.data.frame(table(pw_p689$group))
  pw.prop = as.data.frame(prop.table(table(pw_p689$group)))
  
  p2xkd = unique(c(pw.ncl$Var1[pw.ncl$Freq<10],pw.prop$Var1[pw.prop$Freq<0.01] ))
  sels = names(gs_p689)[!names(gs_p689)%in%p2xkd]
  gs_p689 =subset(gs_p689,select = sels)
  
  return(gs_p689)
}

cnmf_top_genes = cbind(extract_program(sid = 'SP04'),
                       extract_program(sid = 'SP05'),
                       extract_program(sid = 'SP11'),
                       extract_program(sid = 'SP13'),
                       extract_program(sid = 'SP15'),
                       extract_program(sid = 'SP20T5'),
                       extract_program(sid = 'SP21T5'),
                       extract_program(sid = 'SP22'),
                       extract_program(sid = 'SP24'),
                       extract_program(sid = 'SP41')
) #scRNAseq

cnmf_top_genes = cbind(extract_program(sid = 'SP04'),
                       extract_program(sid = 'SP05'),
                       extract_program(sid = 'SP11'),
                       extract_program(sid = 'SP13'),
                       extract_program(sid = 'SP15'),
                       extract_program(sid = 'SP20T5'),
                       extract_program(sid = 'SP21T5'),
                       extract_program(sid = 'SP22'),
                       extract_program(sid = 'SP24'),
                       extract_program(sid = 'SP41')
)



cnmf_intersect <- sapply(cnmf_top_genes, function(x) sapply(cnmf_top_genes, function(y) length(intersect(x, y)))) 
cnmf_intersect = as.data.frame(cnmf_intersect)
# cnmf_intersect.ann = data.frame(Program = names(cnmf_intersect),Sample = '', Protocol = '')
# rownames(cnmf_intersect.ann) = names(cnmf_intersect)
# cnmf_intersect.ann$Sample = substr(cnmf_intersect.ann$Program, 1,4)
# cnmf_intersect.ann$Protocol = ifelse(cnmf_intersect.ann$Sample %in% c('SP20','SP21'),'5\' scRNA', ifelse(cnmf_intersect.ann$Sample %in% c('X172','X177','X747'),'snRNA','3\'scRNA'))
library(pheatmap)
cnmf_intersect.ph<-pheatmap(cnmf_intersect,color = custom_magma, show_colnames = F,border_color = NA,fontsize_row = 5)
cnmf_intersect.ph

library(dendextend)
get_subdendrograms2 <- function(dend, k, ...) {
  clusters <- cutree(dend, k, ...)
  dend_list <- lapply(unique(clusters), function(cluster.id) {
    # bugfix: Added `names(clusters)[]` here
    find_dendrogram(dend, names(clusters)[which(clusters == cluster.id)])
  })
  class(dend_list) <- "dendlist"
  dend_list
}

#' Get the average inter-program similarity (since intra-program similarity = 100%)
avg_similarity <- function(dend) {
  
  # subset the similarity matrix to programs in the provided dendrogram
  x <- cnmf_intersect[labels(dend), labels(dend)]
  
  # set the diagonal to NA to not count intra-program similarity
  diag(x) <- NA
  
  # calculate the mean similarity in the rest of the matrix
  mean( as.matrix(x), na.rm = TRUE)
  
}

#' Given the tree produced by pheatmap::pheatmap(), a function to extract
#' all the metaprograms from the tree
#'
#' We have arbitrarily initialized the thresholds for defining metaprograms.
#' NOTE: some programs will *not* be successfully identified within metaprograms,
#' if they don't meet the criteria. Thus, the total number of programs in the output
#' will be fewer than in number of programs in the input.
#'
#' @param tree Dendrogram
#' @param K Numeric, number of subtrees to cut \code{tree} into at the first cut
#' @param min_similarity Numeric, minimum average similarity of programs within
#' a subtree to define it as a metaprogram
#' @param min_programs Numeric, minimum number of programs within a subtree to
#' define it as a metaprogram
#'
#' @return A list, with one element per metaprogram identified. Each element is a
#' character vector containing the names of the programs in the metaprogram.
define_metaprograms <- function(tree, K = 5, min_similarity = 10, min_programs = 3) {
  
  define_metaprograms_recursive <- function(subtree,
                                            min_programs,
                                            min_similarity,
                                            debug = FALSE) {
    
    # 1. if there are fewer leaves in the tree than the minimum number of
    # programs required to define a metaprogram, drop this subtree
    if (attr(subtree, "members") < min_programs) {
      
      return(NULL)
      
      # 2. if the average similarity within this subtree is greater than the
      # minimum similarity, define this subtree as a metaprogram,
      # and add it to the list
    } else if (avg_similarity(subtree) >= min_similarity) {
      
      metaprograms[[i]] <<- labels(subtree)
      # increment the counter outside the sub-function (scoping assignment)
      i <<- i+1
      
      return(labels(subtree))
      
      # 3. if the subtree is large enough but not similar enough, cut the sutree
      # in 2, and recurse down each child/subsubtree
    } else {
      
      subsubtrees <- get_subdendrograms2(subtree, 2)
      
      lapply(subsubtrees, define_metaprograms_recursive,
             min_programs = min_programs,
             min_similarity = min_similarity)
      
    }
  }
  
  # initialize counter & list
  i <- 1
  metaprograms <- list()
  
  # get the first set of subtrees to initialize S, by cutting it into K subtrees
  S <- get_subdendrograms2(tree, K)
  
  # recurse!
  x <- lapply(S, define_metaprograms_recursive,
              min_programs   = min_programs,
              min_similarity = min_similarity)
  
  # enumerate metaprograms
  names(metaprograms) <- seq_along(metaprograms)
  return(metaprograms)
  
}

cnmf_metaprograms_filt <- define_metaprograms(as.dendrogram(cnmf_intersect.ph$tree_col),K = 3,min_similarity = 10,min_programs = 5)
cnmf_metaprograms_filt

cnmf_intersect.ann = data.frame(Program = names(cnmf_intersect),Sample = '')
rownames(cnmf_intersect.ann) = names(cnmf_intersect)
cnmf_intersect.ann$Sample = substr(cnmf_intersect.ann$Program, 1,nchar(cnmf_intersect.ann$Program)-3)
cnmf_intersect.ann$metaprogram= ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`1`, 'M1',
                                       ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`2`, 'M2',
                                              ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`3`, 'M3',
                                                     ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`4`, 'M4',
                                                            ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`5`, 'M5',
                                                                   ifelse(cnmf_intersect.ann$Program %in%cnmf_metaprograms_filt$`6`, 'M6','others'))))))
#otl = cnmf_intersect.ann$Program[cnmf_intersect.ann$metaprogram=='outlier']
#cnmf_intersect = cnmf_intersect[!rownames(cnmf_intersect)%in%otl,!names(cnmf_intersect)%in%otl]
cnmf_hm_filt <- pheatmap(cnmf_intersect,color = custom_magma, show_colnames = F,border_color = NA,fontsize_row = 5,annotation_col = cnmf_intersect.ann[,-1])

hm_programs_filt <- cnmf_hm_filt$tree_col$labels[cnmf_hm_filt$tree_col$order]

# get column indices
cnmf_metaprograms_filt_idx <- map(cnmf_metaprograms_filt, ~ which(hm_programs_filt %in% .x))

# sort so they're from left to right
cnmf_metaprograms_filt_order <- names(sort(map_dbl(cnmf_metaprograms_filt_idx, 1)))

# put in the right order & rename
cnmf_metaprograms_filt_idx <- cnmf_metaprograms_filt_idx[cnmf_metaprograms_filt_order]

# rename metaprograms
names(cnmf_metaprograms_filt) <- plyr::mapvalues(names(cnmf_metaprograms_filt),
                                                 from = names(cnmf_metaprograms_filt_idx),
                                                 to = seq_along(cnmf_metaprograms_filt_idx))

# rename idx
names(cnmf_metaprograms_filt_idx) <- seq_along(cnmf_metaprograms_filt_idx)

# convert to long data frame
cnmf_metaprograms_filt_df <- imap_dfr(cnmf_metaprograms_filt,
                                      ~ data.frame(Metaprogram = as.numeric(.y), Program = .x,
                                                   stringsAsFactors = FALSE)) %>% 
  arrange(Metaprogram)

hm_metaprograms_filt <- hm_programs_filt[unname(unlist(cnmf_metaprograms_filt_idx))]
length(hm_metaprograms_filt)

#cnmf_intersect.ann$metaprogram = paste0('M',match(cnmf_intersect.ann$Program,cnmf_metaprograms_filt_df$Program))
pheatmap(cnmf_intersect[hm_metaprograms_filt,hm_metaprograms_filt],
         cluster_rows = F, cluster_cols = F,
         gaps_row = cumsum(map(cnmf_metaprograms_filt_idx, length)),
         gaps_col = cumsum(map(cnmf_metaprograms_filt_idx, length)),
         annotation_col = cnmf_intersect.ann[hm_metaprograms_filt,-1],
         color = custom_magma, show_colnames = F,border_color = NA,
         fontsize_row = 5)



#
# correlate with QC metrics
cnmf_intersect.ann = merge(cnmf_intersect.ann, cor_qcpw[,-1], by = 'Program')
rownames(cnmf_intersect.ann) = cnmf_intersect.ann$Program
pheatmap(cnmf_intersect[hm_metaprograms_filt,hm_metaprograms_filt],
         cluster_cols = F, cluster_rows = F, 
         gaps_row = cumsum(map(cnmf_metaprograms_filt_idx, length)),
         gaps_col = cumsum(map(cnmf_metaprograms_filt_idx, length)),
         annotation_col = cnmf_intersect.ann[hm_metaprograms_filt,c(2:7,9,8)],
         color = custom_magma, show_colnames = F,border_color = NA,
         fontsize_row = 5,show_rownames = F)
##
#the markers of each metaprogram  
##
ngenes  = 50
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`1`])) ->tb1;  tb1=tb1[!grepl('^RP[LS]', names(tb1))]; tb1=tb1[!grepl('^A[CL][0-9]{6}', names(tb1))];names(sort(tb1, decreasing = T))[1:ngenes]->mpg1#
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`2`])) ->tb2;  tb2=tb2[!grepl('^RP[LS]', names(tb2))]; tb2=tb2[!grepl('^A[CL][0-9]{6}', names(tb2))];names(sort(tb2, decreasing = T))[1:ngenes]->mpg2
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`3`])) ->tb3;  tb3=tb3[!grepl('^RP[LS]', names(tb3))]; tb3=tb3[!grepl('^A[CL][0-9]{6}', names(tb3))];names(sort(tb3, decreasing = T))[1:ngenes]->mpg3
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`4`])) ->tb4;  tb4=tb4[!grepl('^RP[LS]', names(tb4))]; tb4=tb4[!grepl('^A[CL][0-9]{6}', names(tb4))];names(sort(tb4, decreasing = T))[1:ngenes]->mpg4
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`5`])) ->tb5;  tb5=tb5[!grepl('^RP[LS]', names(tb5))]; tb5=tb5[!grepl('^A[CL][0-9]{6}', names(tb5))];names(sort(tb5, decreasing = T))[1:ngenes]->mpg5
table(unlist(cnmf_top_genes[,cnmf_metaprograms_filt$`6`])) ->tb6;  tb6=tb6[!grepl('^RP[LS]', names(tb6))]; tb6=tb6[!grepl('^A[CL][0-9]{6}', names(tb6))];names(sort(tb6, decreasing = T))[1:ngenes]->mpg6

mpgs <-cbind(mpg1, mpg2, mpg3,mpg4, mpg5, mpg6)

##what are these states?
filbin = read.delim('../../../Filbin_metaprograms.txt')


apply(filbin, 2, function(x) length(intersect(x, mpg1))) #G2M
apply(filbin, 2, function(x) length(intersect(x, mpg2))) #OC-like
apply(filbin, 2, function(x) length(intersect(x, mpg3))) #S
apply(filbin, 2, function(x) length(intersect(x, mpg4))) #AC-like
apply(filbin, 2, function(x) length(intersect(x, mpg5))) #MES-like
apply(filbin, 2, function(x) length(intersect(x, mpg6))) #OPC-like2
#apply(filbin, 2, function(x) length(intersect(x, mpg7))) #OPC-like1

#
cor_with_filbin <- rbind(apply(filbin, 2,function(x) length(intersect(x, mpg1))),
      apply(filbin, 2,function(x) length(intersect(x, mpg2))),
      apply(filbin, 2,function(x) length(intersect(x, mpg3))),
      apply(filbin, 2,function(x) length(intersect(x, mpg4))),
      apply(filbin, 2,function(x) length(intersect(x, mpg5))),
      apply(filbin, 2,function(x) length(intersect(x, mpg6))))

cor_with_filbin = as.data.frame(cor_with_filbin)
rownames(cor_with_filbin) = paste0('M', 1:6)
library(pheatmap)
pheatmap(t(cor_with_filbin[c(3,1,4,5,2,6),c("S","G2M","AC.like","MES.like","OC.like","OPC.like.1","OPC.like.2","OPC.like.3")]),
         cluster_rows = F,cluster_cols = F, border_color = 'white',display_numbers = T,number_format = '%.0f')

write.table(mpgs, file = 'K27scDMG.cellular.states.6sates.202502.txt', quote = F, sep = "\t", row.names = F)

#merge OPC-like 1-3 and overlap, plot heatmap
mat = as.data.frame(matrix(0,8,6))
mg = as.data.frame(cbind(mpg1,mpg2,mpg3,mpg4,mpg5,mpg6))
names(mat) = names(mg)
rownames(mat) = names(filbin)
for (i in 1:6){
  for (j in 1:8){
    mat[j,i] = length(intersect(mg[,i], filbin[,j]))
  }
}

mat = as.data.frame(t(mat))
mat$OPC.like = mat$OPC.like.1+ mat$OPC.like.2 + mat$OPC.like.3

#mat = mat[,names(mg)]


pheatmap::pheatmap(mat[c(3,1,4,5,6,2),c('S','G2M','AC.like','MES.like','OPC.like','OC.like')],clustering_method = 'average',
                   display_numbers=F,number_format = "%.0f",cluster_rows = F,cluster_cols = F,border_color = 'white',legend = T)

#assign cell states using the identified metaprograms
metaprograms = list()
for (i in 1:6){
  metaprograms[[i]] = mpgs[,i]
}
names(metaprograms) = paste('MP',1:6,sep = "_")

metaprograms_df = as.data.frame(metaprograms)
names(metaprograms_df) = c('G2M','OClike','S','AClike','MESlike','OPClike')
write.table(metaprograms_df, file = '~/Documents/Projects/H3DSG/Revision/metaprograms.cNMF.n6.20250208.txt',row.names = F, quote =F, sep = "\t")

library(Seurat)
load('~/Dropbox/communter/seurat_integrated.Rda')
#s0 = seurat_integrated
seurat_integrated = subset(seurat_integrated, subset = orig.ident!='SP20B5' & orig.ident!='SP21B5')
m0 = seurat_integrated@meta.data
seurat_integrated = subset(seurat_integrated, subset = cell.state %in% c('AC-like','OC-like','MES-like','OPC-like','G2M','S') )


seurat_integrated <- AddModuleScore(seurat_integrated, features = metaprograms,name = 'MP')


m = seurat_integrated@meta.data
m$TumorGroup = ifelse(m$orig.ident %in% c('SP24',"SP41","SP15","SP04"),'Pons-like','Thalamus-like')


mpnames = c("G2M","OC-like","S","AC-like","MES-like","OPC-like")
m$State2025 = apply(m[,28:33],1, function(x) mpnames[which.max(x)])
m$State2025 = factor(m$State2025,levels = c('S','G2M','OPC-like','AC-like','MES-like','OC-like'))
m1 = m[m$cell.state %in% c('AC-like','OC-like','MES-like','OPC-like','G2M','S'),]

pheatmap(table(m1$State2025, m1$cell.state),scale = 'row')


prop.table(table(m1$cell.state,m1$TumorGroup),margin = 2)
prop.table(table(m1$State2025,m1$TumorGroup),margin = 2)

seurat_integrated$State2025 = m1$State2025

DimPlot(seurat_integrated, group.by = 'State2025')

m0$State2025 = ifelse(m0$cell.barcode %in% m1$cell.barcode, as.character(m1$State2025)[match(m0$cell.barcode, m1$cell.barcode)], m0$cell.state)
m0$State2025[m0$State2025 %in% c('mural','endothellial')] = 'vascular'
m0$State2025[m0$State2025 %in% c('S','G2M')] = 'cycling'
m0$State2025[m0$State2025 %in% c('OPC-like1','OPC-like2')] = 'OPC-like'
m0$State2025[m0$State2025 %in% c('unknown')] = NA
m0$TumorGroup = ifelse(m0$orig.ident %in% c('SP24',"SP41","SP15","SP04"),'Pons-like','Thalamus-like')

prop.table(table(m0$State2025,m0$TumorGroup),margin = 2)

df1 = as.data.frame(prop.table(table(m0$State2025,m0$TumorGroup),margin = 2))
df1$Var1=factor(df1$Var1, levels = c('cycling','OPC-like','AC-like','MES-like','OC-like','vascular','T','B','myeloid','neutrophil','oligodendrocyte'))

ggplot(df1, aes(x = Var2, y = 100*Freq, fill = Var1))+
  geom_bar(stat = 'identity',width = 0.7)+
  scale_fill_manual(values = c("#d53e4f","#f46d43","#fdae61",'#fcc5c0',"#fee08b",'#35978f','#3c75af','#519e3e','#66c2a5','#a6dba0','#ffffb3'))+
  scale_y_continuous(expand = c(0,0))+theme_classic()+
  theme(axis.text = element_text(color = 'black'),legend.key.size = unit(0.5, "cm"))+
  labs(x = '',y='Proportion of cells (%)',fill = '')

library(scCustomize)
random_cells_150 <- Random_Cells_Downsample(seurat_object = seurat_integrated, num_cells = 150,
                                            allow_lower = T,group.by = 'State2025')
#
ramdom1kcells = sample(Cells(seurat_integrated),2000)
gns = unique(c(mpg3,mpg1,mpg6,mpg4,mpg5,mpg2))
library(RColorBrewer)
my_palette <- rev(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"))
seurat_integrated$State20252 = ifelse(seurat_integrated$State2025 %in% c('S','G2M'),'cycling',as.character(seurat_integrated$State2025))
seurat_integrated$State20252 = factor(seurat_integrated$State20252,levels=c('cycling','OPC-like','AC-like','MES-like','OC-like'))
DoHeatmap(seurat_integrated,features = gns,group.by = 'State2025',cells = ramdom1kcells,
          angle = 30,size = 3, #slot = 'data',
          group.colors = c("#d53e4f","#d53e4f","#f46d43","#fdae61",'#fcc5c0',"#fee08b"))+
  scale_fill_gradientn(colours = my_palette)+
  theme(axis.text.y = element_blank(),legend.position = "bottom")+guides(color = 'none')

#dot plot on all cells
load('~/Dropbox/communter/seurat_integrated.Rda')
seurat_integrated$State2025 = m0$State2025[match(seurat_integrated$cell.barcode, m0$cell.barcode)]

cells_pons = Cells(seurat_integrated)[!is.na(seurat_integrated$State2025) &seurat_integrated$orig.ident %in% c('SP24',"SP41","SP15","SP04") ]
cells_thalamus = Cells(seurat_integrated)[!is.na(seurat_integrated$State2025) &seurat_integrated$orig.ident %in% c('SP05','SP11','SP13','SP20T5','SP21T5','SP22') ]

seurat_integrated$State20252 = ifelse(seurat_integrated$State2025 %in% c('AC-like','MES-like'),'AC/MES-like',
                                      ifelse(seurat_integrated$State2025 %in% c('OPC-like','OC-like'),'OPC/OC-like',seurat_integrated$State2025))
p1 <-DotPlot(subset(seurat_integrated,cells  = cells_pons ),
        group.by = 'State20252', features = c('CCL8','CCL19','CCL3','CCL4','CCL2','CXCL14'))

p2 <-DotPlot(subset(seurat_integrated,cells  = cells_thalamus ),
        group.by = 'State20252', features = c('CCL8','CCL19','CCL3','CCL4','CCL2','CXCL14'))
data_p1 = p1$data
data_p2 = p2$data

data_p1$tumor = 'Pons-like'
data_p2$tumor = 'Thalamus-like'

data_p12 = rbind(data_p1, data_p2)
data_p12$id = factor(data_p12$id, levels = rev(c('OPC/OC-like','AC/MES-like','cycling','vascular','T','B','myeloid','neutrophil','oligodendrocyte')))
ggplot(data_p12[data_p12$features.plot=='CXCL14',])+
  geom_point(aes(x = tumor, y = id, color = avg.exp, size =pct.exp ))+
  scale_color_gradient2(low = 'white',high = 'blue')+
  scale_size_continuous(range = c(0.1,6),limits = c(0,30),breaks = c(0,15,30))+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'),axis.text.x = element_text(angle = 30, hjust=1))+
  labs(x = '', y = '')+
  guides(color = 'none')



library(scCustomize)
seurat_integrated$tumorGroup = ifelse(seurat_integrated$orig.ident %in% c('SP'))
FeaturePlot_scCustom(seurat_object = seurat_integrated, features = "CXCL14",order = T,split.by = '')

