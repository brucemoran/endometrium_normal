##from: http://www.bioconductor.org/packages/release/bioc/vignettes/Starr/inst/doc/Starr.jpeg
.libPaths("/home/bmoran/R/x86_64-pc-linux-gnu-library/3.3")
library(Starr)

wd<-c("/store2/bmoran/cardio/chris_meDip/Chris_Watson_5MeC_Array_data_posted/")
owd<-paste0(wd,"/ringo-starr_meDIP")
dir.create(owd,showWarnings=FALSE)
setwd(owd)

##read data
bpmap <- readBpmap(paste0(wd,"/0419UCD\ CEL\ files/Hs_PromPR_v02-3_NCBIv36.bpmap"))

cels <- c(paste0(wd,"/0419UCD\ CEL\ files/0419UCD_MetDNA_SampleA_IP_human_prom_090113.CEL"),
  paste0(wd,"/0419UCD\ CEL\ files/0419UCD_MetDNA_SampleB_IP_human_prom_090113.CEL"),
  paste0(wd,"/0419UCD\ CEL\ files/0419UCD_MetDNA_SampleC_IP_human_prom_090113.CEL"),
  paste0(wd,"/0419UCD\ CEL\ files/0419UCD_MetDNA_SampleD_IP_human_prom_090113.CEL"))

names <- c("A","B","C","D")
type <- c("NORMOXIA","NORMOXIA","HYPOXIA","HYPOXIA")

##Annotations: make probeAnno, transcriptAnno
probeAnno <- bpmapToProbeAnno(bpmap,verbose=T,uniqueSeq=T)
transcriptAnno <- read.gffAnno(paste0(wd,"/human_hg18_2010_11_10/hg18_gene-cpg-tss.gff"))
transcriptAnno[,11]<-paste0("Hs:NCBIv36;",transcriptAnno[,1])
transcriptAnno[,1]<-transcriptAnno[,11];
transcriptAnno<-transcriptAnno[,-11]
tssAnno<-transcriptAnno[transcriptAnno$feature=="Transcription_Start_Site",]

##load CEL files
celExpst <- readCelFile(bpmap, cels, names, type, featureData=T, log.it=T)
##LOESS normalise
celExpst_loess <- normalize.Probes(celExpst, method="loess")

##QC plots showing
pdf("plotDensity-Boxes_Starr_ABCD.pdf")
par(mfcol=c(1,2))
plotDensity(celExpst, oneDevice=T, main="")
plotBoxes(celExpst)
dev.off()

pdf("plotGC-Pos_Biases_Starr_ABCD.pdf")
par(mfcol=c(1,2))
plotGCbias(exprs(celExpst)[,1], featureData(celExpst)$seq, main="")
plotPosBias(exprs(celExpst)[,1], featureData(celExpst)$seq)
dev.off()

pdf("scatter-density_Starr_ABCD.pdf")
plotScatter(celExpst, density=T, cex=0.5)
dev.off()

#########################
## HYPOXIA vs NORMOXIA ##
#########################
##set conditions
ipHypoxia <- celExpst$type == "HYPOXIA"
controlNormoxia <- celExpst$type == "NORMOXIA"

##calculate ratios between ip, controls on LOESS normalised data
celExpst_loess_ratio_hypoxia <- getRatio(celExpst_loess, ipHypoxia, controlNormoxia, description="Hypoxia-vs-Normoxia", fkt=median, featureData=F)

##smooth ratios and name appropriately
celExpst_ratio_smooth_hypoxia <- computeRunningMedians(celExpst_loess_ratio_hypoxia, probeAnno=probeAnno, verbose=T, modColum="type")
sampleNames(celExpst_ratio_smooth_hypoxia) <- paste(sampleNames(celExpst_loess_ratio_hypoxia),"smoothed")

##QC plots for normalisation
pdf("hypoxia_plotMA_unLoessStarr_ABCD.pdf")
plotMA(celExpst, ip=ipHypoxia, control=controlNormoxia)
description <- c("Hypoxia vs. Normoxia (not Loess normalised)")
dev.off()

pdf("hypoxia_plotMA_LoessStarr_ABCD.pdf")
plotMA(celExpst_loess, ip=ipHypoxia, control=controlNormoxia)
description <- c("Hypoxia vs. Normoxia (Loess normalised)")
dev.off()

##remove all unknowns in smoothed ratios (e.g. where values are not available for one condition)
y0_hypoxia <- apply(exprs(celExpst_ratio_smooth_hypo), 2, upperBoundNull)

##define how far apart CHER (chromatin enriched regions) must be to collapse into a single
distCutOff_hypoxia <- 100000
findChers_hypoxia <- findChersOnSmoothed(celExpst_ratio_smooth_hypoxia, probeAnno=probeAnno, thresholds=y0_hypoxia, distCutOff=distCutOff_hypoxia, verbose=T, minProbesInRow = 10)
relateChers_hypoxia <- relateChers(findChers_hypoxia, transcriptAnno, upstream=500)
chers_hypoxia <- as.data.frame.cherList(relateChers_hypoxia)
chers_hypoxia <- chers_hypoxia[order(chers_hypoxia$maxLevel, decreasing=TRUE),]

write.table(chers_hypoxia,"_hypoxia_chers.txt",sep="\t",quote=F)

#########################
## NORMOXIA vs HYPOXIA ##
#########################
##set conditions
ipNormoxia <- celExpst$type == "NORMOXIA"
controlHypoxia <- celExpst$type == "HYPOXIA"
##set conditions
ipHypoxia <- celExpst$type == "HYPOXIA"
controlNormoxia <- celExpst$type == "NORMOXIA"

##calculate ratios between ip, controls on LOESS normalised data
celExpst_loess_ratio_normoxia <- getRatio(celExpst_loess, ipNormoxia, controlnormoxia, description="Normoxia-vs-Hypoxia", fkt=median, featureData=F)

##smooth ratios and name appropriately
celExpst_ratio_smooth_normoxia <- computeRunningMedians(celExpst_loess_ratio_normoxia, probeAnno=probeAnno, verbose=T, modColum="type")
sampleNames(celExpst_ratio_smooth_normoxia) <- paste(sampleNames(celExpst_loess_ratio_normoxia),"smoothed")

##QC plots for normalisation
pdf("normoxia_plotMA_unLoessStarr_ABCD.pdf")
plotMA(celExpst, ip=ipNormoxia, control=controlHypoxia)
description <- c("Normoxia vs. Hypoxia (not Loess normalised)")
dev.off()

pdf("normoxia_plotMA_LoessStarr_ABCD.pdf")
plotMA(celExpst_loess, ip=ipNormoxia, control=controlHypoxia)
description <- c("Normoxia vs. Hypoxia (Loess normalised)")
dev.off()

##remove all unknowns in smoothed ratios (e.g. where values are not available for one condition)
y0_normoxia <- apply(exprs(celExpst_ratio_smooth_hypo), 2, upperBoundNull)

##define how far apart CHER (chromatin enriched regions) must be to collapse into a single
distCutOff_normoxia <- 100000
findChers_normoxia <- findChersOnSmoothed(celExpst_ratio_smooth_normoxia, probeAnno=probeAnno, thresholds=y0_normoxia, distCutOff=distCutOff_normoxia, verbose=T, minProbesInRow = 10)
relateChers_normoxia <- relateChers(findChers_normoxia, transcriptAnno, upstream=500)
chers_normoxia <- as.data.frame.cherList(relateChers_normoxia)
chers_normoxia <- chers_normoxia[order(chers_normoxia$maxLevel, decreasing=TRUE),]

write.table(chers_normoxia,"_normoxia_chers.txt",sep="\t",quote=F)

##save above in RData file
save.image("_Starr_analysis_hypoxia-normoxia.RData")

##annotation using sh script of the same name
system("sh Ringo_Starr_meDIP.sh")
