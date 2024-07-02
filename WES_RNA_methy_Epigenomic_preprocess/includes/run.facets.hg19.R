args=commandArgs(T)
sampleid<- args[1]
cval.choose<- args[2]
# sampleid<- "TCGA-12-3648-01A-VS-TCGA-12-3648-10A"

library("facets")
set.seed(2020)

XX.dat<- readSnpMatrix(paste(sampleid,".mpileup",sep = ""),perl.pileup=T)
xx=preProcSample(XX.dat,gbuild="hg19")

# cval.choose<- 350
oo=procSample(xx,cval=as.numeric(cval.choose))

# oo=procSample(xx,cval=150,dipLogR=-0.02)
pdf(file = paste(sampleid,".cval_",cval.choose,".logRlogORspider.pdf",sep=""),width = 10,height = 10)
logRlogORspider(oo$out,oo$dipLogR)
dev.off()

oo$dipLogR
# oo$dipLogR<- -0.18
oo_out<- oo$out
oo_out$dipLogR<- oo$dipLogR

if(!is.null(oo$flags)){
  oo_out$flags<- paste(as.character(oo$flags),collapse="; ")
}else{
  oo_out$flags<- "No flags"
}


write.table(oo_out,file = paste(sampleid,".cval_",cval.choose,".oo.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)


fit=emcncf(oo,trace=T)
fit_cncf<- fit$cncf
fit_cncf$seglen<- fit$seglen
fit_cncf$dipLogR<- fit$dipLogR
fit_cncf$purity<- fit$purity
fit_cncf$ploidy<- fit$ploidy
fit_cncf<- fit_cncf[,c("seg","chrom","start","end","num.mark","nhet","dipLogR","purity","ploidy","cnlr.median","mafR","segclust","cnlr.median.clust","mafR.clust","cf.em","tcn.em","lcn.em")]#,"seglen"
write.table(fit_cncf,file = paste(sampleid,".cval_",cval.choose,".fit_cncf.txt",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)

fit_cncf.seg<- fit_cncf[,c("seg","chrom","start","end","num.mark","cnlr.median")]
fit_cncf.seg$seg<- sampleid
write.table(fit_cncf.seg,file = paste(sampleid,".cval_",cval.choose,".seg",sep = ""),col.names = T,row.names = F,sep = "\t",quote = F)

fit_purity<- fit$purity
fit_ploidy<- fit$ploidy
cat(file = paste(sampleid,".cval_",cval.choose,".purity.and.ploidy.txt",sep = ""),paste("sampleid:",sampleid,"purity:",fit_purity,"ploidy:",fit_ploidy,"\n",sep = "\t"))

pdf(file = paste(sampleid,".cval_",cval.choose,".CNV.profile.pdf",sep=""),width = 15,height = 15)
sname<- sprintf('%s; ploidy= %.2f; purity= %.2f',sampleid,fit_ploidy , fit_purity)
# plotSample(x= facets$proc_out, emfit= facets$emcncf_fit, sname= sname)
plotSample(x=oo,emfit=fit,sname= sname)
dev.off()

png(file = paste(sampleid,".cval_",cval.choose,".CNV.profile.png",sep=""), units="px", width=1600, height=1600, res=300)
sname<- sprintf('%s; ploidy= %.2f; purity= %.2f',sampleid,fit_ploidy , fit_purity)
# plotSample(x= facets$proc_out, emfit= facets$emcncf_fit, sname= sname)
plotSample(x=oo,emfit=fit,sname= sname)
dev.off()

