## Load libraries and datasets
library("plotgardener")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("plotgardenerData")
library("AnnotationHub")

#
setwd('~/Dropbox/Mac (2)/Documents/Research/SpinalcordGlioma/HOX/')
#the reference genome was hg19
#hoxa: chr7:27130000-27239800
#hoxb: chr17:46600000-46780000
#hoxc: chr12:54330000-54451000
#hoxd: chr2:176950000-177060000
hoxs = data.frame(family=c('HOXA','HOXB','HOXC','HOXD'),chr = c('chr7','chr17','chr12','chr2'),strand =c('-','-','+','+'),
                  st = c(27130000,46600000,54330000,176950000),end = c(27239800,46810000,54451000,177060000))
hoxgns = list(HOXA = paste0('HOXA',c(1:7,9:11,13)),
              HOXB = paste0('HOXB',c(1:9,13)),
              HOXC = paste0('HOXC',c(4:6,8:13)),
              HOXD = paste0('HOXD',c(1,3,4,8:13)))

atac_file = list.files('bulkATAC/',pattern = 'bw',full.names = T)
k27me3_file = list.files('bulkH3K27Me3//',pattern = 'bw',full.names = T)
atac_file = atac_file[c(3,5)] # H170008 group A,H177793 group B
k27me3_file = k27me3_file[c(3,5)] # H170008 group A,H177793 group B

pageCreate(width = 9, height = 4.5, default.units = "inches")
for (i in 1:4){ #HOXA, HOXB, HOXC, HOXD
  gni = hoxs$family[i]
  chri = hoxs$chr[hoxs$family==gni];sti = hoxs$st[hoxs$family==gni];endi = hoxs$end[hoxs$family==gni]
  sdi = hoxs$strand[hoxs$family==gni]; if(sdi=="+"){colsi=c('#669fd9','white')}else{colsi=c('white','#669fd9')}
  
  params_hoxc <- pgParams(chrom = chri,chromstart = sti,chromend = endi, 
                          assembly = "hg19", default.units = "inches",x =0.1+2.2*(i-1))
  
  for (j in 1:2){
    ## Plot ATAC signal
    tdf_atac = readBigwig(file = atac_file[j],
                          chrom = chri,chromstart = sti,chromend = endi)
    tc_range <- pgParams(range = c(0, max(tdf_atac$score)),assembly = "hg19")
    atac_imr <- plotSignal(data = tdf_atac, params = c(params_hoxc, tc_range),
                           fill = "#008f00", linecolor = "#008f00",
                           y = 0.1+1.5*(j-1), height = 0.6, width = 2, just = c("left", "top"))
    ## Plot H3K27Me3 signal
    tdf_k27me3 = readBigwig(file = k27me3_file[j],
                            chrom = chri,chromstart = sti,chromend = endi)
    hk_range <- pgParams(range = c(0, max(tdf_k27me3$score)),assembly = "hg19")
    k27me3_imr <- plotSignal(data = tdf_k27me3, params = c(params_hoxc, hk_range),
                             fill = "#253494", linecolor = "#253494",
                             y = 0.8+1.5*(j-1), height = 0.6, width = 2,just = c("left", "top"))
  }
  
  ## Plot gene track
  genes_imr <- plotGenes(params = params_hoxc, stroke = 1, fontsize = 6,
                         geneOrder = hoxgns[gni],
                         strandLabels = FALSE,fontcolor = colsi,fill = colsi,
                         y = 3.1, height = 0.4, width = 2, just = c("left", "top"))
  ## Annotate genome label
  annoGenomeLabel(plot = genes_imr, params = params_hoxc, 
                  scale = "Kb", fontsize = 7, digits = 0,
                  y = 3.6, just = c("left", "top"))
}

## Hide page guides
pageGuideHide()


# #
# # PAX3 and NKX6-1
# #
pax3 =list('chr2',223029037, 223181035) 
nkx61 = list('chr4',85408683,85425782)

aa=c('bulkATAC/H182545_new.sorted.noorg_Norm_new_CPM_cover.bw',
     'bulkATAC/H179186_new.sorted.noorg_Norm_new_CPM_cover.bw',
'bulkATAC/H164921_new.sorted.noorg_Norm_new_CPM_cover.bw',
'bulkATAC/H177793_new.sorted.noorg_Norm_new_CPM_trans_cover.bw')

km =c('bulkH3K27Me3/H182545_bowtie2.sorted.rmdup.mapped.rmChrM.Norm_new_CPM_cover.bw',
      'bulkH3K27Me3//H179186_bowtie2.sorted.rmdup.mapped.rmChrM.Norm_new_CPM_cover.bw',
      'bulkH3K27Me3/H164921_bowtie2.sorted.rmdup.mapped.rmChrM.Norm_new_CPM_cover.bw',
      'bulkH3K27Me3/H177793_bowtie2.sorted.rmdup.mapped.rmChrM.Norm_new_CPM_cover.bw')


pageCreate(width = 4.5, height = 3.2, default.units = "inches")
params_pax3 <- pgParams(chrom = 'chr2',chromstart = 223030000,chromend = 223180000, 
                        assembly = "hg19", default.units = "inches",x =0.1)
params_nkx6.1 <- pgParams(chrom = 'chr4',chromstart = 85410000,chromend = 85425000, 
                        assembly = "hg19", default.units = "inches",x =2.3)
tcmaxes =c(2,3,1.25,2); kmmaxes=c(24,24,6,20); colsi=c('white','#669fd9')
for (ix in 1:4){ #PAX3
  ## plot atac
  tdf_atac = readBigwig(file = aa[ix],chrom = 'chr2',chromstart = 223030000,chromend = 223180000,)
  tc_range <- pgParams(range = c(0, tcmaxes[ix]),assembly = "hg19")
  atac_imr <- plotSignal(data = tdf_atac, params = c(params_pax3, tc_range),
                         fill = "#008f00", linecolor = "#008f00",#scale=T,
                         y = 0.1+0.7*(ix-1), height = 0.3, width = 2, just = c("left", "top"))
  ## Plot H3K27Me3 signal
  tdf_k27me3 = readBigwig(file = km[ix],chrom = 'chr2',chromstart = 223030000,chromend = 223180000,)
  hk_range <- pgParams(range = c(0, kmmaxes[ix]),assembly = "hg19")
  k27me3_imr <- plotSignal(data = tdf_k27me3, params = c(params_pax3, hk_range),
                           fill = "#253494", linecolor = "#253494",#scale=T,
                           y = 0.45+0.7*(ix-1), height = 0.3, width = 2,just = c("left", "top"))
}
## Plot gene track
genes_imr <- plotGenes(params = params_pax3, stroke = 1, fontsize = 6,
                       strandLabels = FALSE,fontcolor = colsi,fill = colsi,
                       y = 2.85, height = 0.2, width = 2, just = c("left", "top"))
## Annotate genome label
annoGenomeLabel(plot = genes_imr, params = params_pax3, 
                scale = "Kb", fontsize = 7, digits = 0,
                y = 3.1, just = c("left", "top"))

for (ix in 1:4){ #NKX6-1
  ## plot atac
  tdf_atac = readBigwig(file = aa[ix],chrom = 'chr4',chromstart = 85410000,chromend = 85425000,)
  tc_range <- pgParams(range = c(0, tcmaxes[ix]),assembly = "hg19")
  atac_imr <- plotSignal(data = tdf_atac, params = c(params_nkx6.1, tc_range),
                         fill = "#008f00", linecolor = "#008f00",#scale=T,
                         y = 0.1+0.7*(ix-1), height = 0.3, width = 2, just = c("left", "top"))
  ## Plot H3K27Me3 signal
  tdf_k27me3 = readBigwig(file = km[ix],chrom = 'chr4',chromstart = 85410000,chromend = 85425000,)
  hk_range <- pgParams(range = c(0, kmmaxes[ix]),assembly = "hg19")
  k27me3_imr <- plotSignal(data = tdf_k27me3, params = c(params_nkx6.1, hk_range),
                           fill = "#253494", linecolor = "#253494",#scale=T,
                           y = 0.45+0.7*(ix-1), height = 0.3, width = 2,just = c("left", "top"))
}
## Plot gene track
genes_imr <- plotGenes(params = params_nkx6.1, stroke = 1, fontsize = 6,
                       strandLabels = FALSE,fontcolor = colsi,fill = colsi,
                       y = 2.85, height = 0.2, width = 2, just = c("left", "top"))
## Annotate genome label
annoGenomeLabel(plot = genes_imr, params = params_nkx6.1, 
                scale = "Kb", fontsize = 7, digits = 0,
                y = 3.1, just = c("left", "top"))
pageGuideHide()


