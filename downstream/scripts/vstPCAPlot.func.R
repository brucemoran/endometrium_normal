##produce PCA plot from dds after applying VST

vstPCAPlot <- function(vst, intgroup, PLOTDIR, TAG, ntop=NULL, savePlot=NULL){
  if(class(vst)!="DESeqTransform"){
    print(paste0(vst, " is not of type DESeqTransform, please run rlog() or varianceStabilizingTransformation()"))
  }
  else{
  if(is.null(ntop)){
    ntop <- 500
  }

  if(!is.numeric(ntop)){
    print("Require a numeric for ntop")
    break
  }

  if(is.null(savePlot)){
    savePlot <- TRUE
  }

  print(paste0("PCA of VST on: ", TAG , "..."))
  pcaData <- DESeq2::plotPCA(vst, intgroup = intgroup, returnData=T, ntop=ntop)
  for (x in 1:length(intgroup)){
    pcaData[intgroup[x]] <- as.vector(unlist(pcaData[intgroup[x]]))
    if(is.numeric(pcaData[intgroup[x]][,1])){
      pcaData[intgroup[x]] <- as.character(pcaData[intgroup[x]][,1])
    }
  }
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  colnames(pcaData)[colnames(pcaData) %in% intgroup] <- paste("intgroup",1:length(intgroup), sep="_")
  if(length(intgroup)>=2){
    pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, shape = intgroup_1, colour = intgroup_2)) +
                scale_shape_manual(values=c(1:length(unique(c(pcaData$intgroup_1))))+1) +
                geom_point(size=3) +
                xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                coord_fixed() +
                labs(title=TAG, subtitle=paste0("top ", ntop, " genes by variance"), shape = intgroup[1], colour = intgroup[2])
  }
  if(length(intgroup)==1){
    pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, shape = intgroup_1)) +
                scale_shape_manual(values=c(1:length(unique(c(pcaData$intgroup_1))))+1) +
                geom_point(size =3) +
                xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                coord_fixed() +
                labs(title=TAG, subtitle=paste0("top ", ntop, " genes by variance"), shape=intgroup[1])
  }
  if(savePlot==TRUE){
    ggsave(filename=paste0(PLOTDIR, "/", TAG, ".ntop_", ntop, ".deseq2_pca.pdf"), plot=pcaPlot)
  }
  return(pcaPlot)
}}
