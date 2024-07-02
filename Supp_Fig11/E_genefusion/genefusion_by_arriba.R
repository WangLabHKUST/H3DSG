rm(list = ls())
library(rstudioapi)
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())


fusions_GA <- readRDS('fusion_GA_rbind.RDS')
fusions_GB <- readRDS('fusion_GB_rbind.RDS')
source('support_function.R')
library(circlize)
library(GenomicRanges)
cyto_file <- 'cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv'
cytobands <- read.table(cyto_file , header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
minConfidenceForCircosPlot="medium"




fusions_GA_rbind <- fusions_GA
pdf('Group_1_gene_fusion.pdf',width=6,height = 6)
circos.clear()
circos.initializeWithIdeogram(cytoband=cytobands, plotType="labels")
circos.genomicIdeogram(cytoband=cytobands)
circosColors <- c(translocation="grey80",#"#7E6148B2",#"#000000", 
                  duplication="#00A087B2",#"#00bb00", 
                  deletion="#DC0000B2",#"#ff0000", 
                  inversion="slateblue2"#"#0000ff"
)
for(i in 1:dim(fusions_GA_rbind)[1])
{
  if(fusions_GA_rbind$type[i]=="translocation")
  {
    co <- "grey80"
  }else{
    if(fusions_GA_rbind$type[i]=="duplication")
    {
      co <- "#00A087B2"
    }else{
      if(fusions_GA_rbind$type[i]== "deletion")
      {
        co <-"#DC0000B2"
      }else{
        if(fusions_GA_rbind$type[i]== "inversion")
        {
          co <-"slateblue2"
        }
      }
    }
  }
  circos.link(
    (fusions_GA_rbind$contig1[i]),as.numeric(fusions_GA_rbind$breakpoint1[i]),
    (fusions_GA_rbind$contig2[i]),as.numeric(fusions_GA_rbind$breakpoint2[i]),
    lwd=2, col=co
  )
}
dev.off()

fusion_GB_rbind <- fusions_GB# (fusion_rbind(fusions_GB))
pdf('Group_2_gene_fusion.pdf',width=6,height = 6)
circos.clear()
circos.initializeWithIdeogram(cytoband=cytobands, plotType="labels")
circos.genomicIdeogram(cytoband=cytobands)
circosColors <- c(translocation="grey80",#"#7E6148B2",#"#000000", 
                  duplication="#00A087B2",#"#00bb00", 
                  deletion="#DC0000B2",#"#ff0000", 
                  inversion="slateblue2"#"#0000ff"
)
for(i in 1:dim(fusion_GB_rbind)[1])
{
  if(fusion_GB_rbind$type[i]=="translocation")
  {
    co <- "grey80"#"#7E6148B2" #"#000000"##grey
  }else{
    if(fusion_GB_rbind$type[i]=="duplication")
    {
      co <- "#00A087B2"#"#00bb00"##green
    }else{
      if(fusion_GB_rbind$type[i]== "deletion")
      {
        co <-"#DC0000B2"#"#ff0000"##pink 
      }else{
        if(fusion_GB_rbind$type[i]== "inversion")
        {
          co <-"slateblue2"#"#4DBBD5B2"# "#0000ff"## purple
        }
      }
    }
  }
  circos.link(
    (fusion_GB_rbind$contig1[i]),as.numeric(fusion_GB_rbind$breakpoint1[i]),
    (fusion_GB_rbind$contig2[i]),as.numeric(fusion_GB_rbind$breakpoint2[i]),
    lwd=2, col=co
  )
}
dev.off()
