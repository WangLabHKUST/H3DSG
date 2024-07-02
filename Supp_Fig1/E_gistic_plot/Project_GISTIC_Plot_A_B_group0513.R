rm(list=ls())
library(rstudioapi)
library("ggplot2")
library("gridExtra")
library(grid)
# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())
library(maftools)
all.lesions <-'all_lesions.conf_90.txt' #system.file("extdata", "all_lesions.conf_99.txt", package = "maftools")
amp.genes <- 'amp_genes.conf_90.txt'#system.file("extdata", "amp_genes.conf_99.txt", package = "maftools")
del.genes <- 'del_genes.conf_90.txt'#system.file("extdata", "del_genes.conf_99.txt", package = "maftools")
scores.gis <- 'scores-2.gistic'#system.file("extdata", "scores.gistic", package = "maftools")

laml.gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, 
                         gisticDelGenesFile = del.genes,
                         gisticScoresFile = scores.gis,
                         isTCGA = TRUE)
pdf('GISTIC_maftools_thresh.pdf',width = 12,height = 6)
gisticChromPlot(gistic = laml.gistic, markBands = "all",
                txtSize = 1,fdrCutOff=0.12,
                cytobandTxtSize = 0.8
                )
dev.off()






