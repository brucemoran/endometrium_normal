##general utilities for DE analysis and other RNAseq things

libs <- c("DESeq2", "apeglm", "ggplot2", "tidyverse", "dplyr", "biomaRt", "gtools")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annoBM <- function(attrs, mart, filts=NULL, vals=NULL,  retone=NULL, arrone=NULL){
  if(is.null(vals)){
    vals <- ""
  }
  if(is.null(filts)){
    filts <- ""
  }

  gBM <- biomaRt::getBM(attributes=attrs, filters=filts, values=vals, mart = mart)
  if(is.null(retone)){
    if(is.null(arrone)){
      return(as_tibble(gBM))
    }
    if(!is.null(arrone)){
      return(as_tibble(gBM) %>%
      dplyr::arrange_(arrone))
    }
  }
  if(!is.null(retone)){
    if(is.null(arrone)){
      return(as_tibble(gBM) %>%
        dplyr::filter_(retone))
    }
    if(!is.null(arrone)){
      return(as_tibble(gBM, col_names=list(vals)) %>%
        dplyr::select_(retone) %>%
        dplyr::arrange_(arrone))
    }
  }
}

#########
## TPM ##
#########
##countData is a data frame or tibble
##geneCol is name of column with geneIDs to look up
convertToTpm <- function(countData, geneCol) {
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  if(!is.tibble(countData)){
    countData <- as_tibble(countData, rownames=geneCol)
  }
  if(is_tibble(countData)){
    vals <- as.vector(unlist(countData[,colnames(countData) %in% geneCol]))
    lengthsAll <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position"),
                              values=countData[[geneCol]], mart = mart) %>%
               dplyr::arrange(external_gene_name) %>%
               dplyr::mutate(gene_length = end_position - start_position) %>%
               dplyr::select(!!geneCol, gene_length)
    lengthData <- left_join(countData, lengthsAll, by = geneCol) %>%
                   dplyr::select(genecol=!!geneCol, gene_length) %>%
                   dplyr::group_by(genecol) %>%
                   dplyr::filter(gene_length == max(gene_length)) %>%
                   unique()

    countDataInt <- countData %>% mutate_if(is.numeric, as.integer) %>%
                                  dplyr::rename(genecol=!!geneCol) %>%
                                  dplyr::filter(genecol %in% c(lengthData$genecol)) %>%
                                  dplyr::select(-genecol)

    rateData <- countDataInt / c(lengthData$gene_length)
    rateData[is.na(rateData)] <- 0
    rateData <- rateData / sum(rateData) * 1e6
    tpmData <- as_tibble(cbind(lengthData[["genecol"]],
                               rateData)) %>%
               dplyr::rename(!!geneCol := "lengthData[[\"genecol\"]]")
    return(tpmData)
  }
}

###############
## Factorise ##
###############
factorise <- function(f){
  for(i in 1:ncol(f)){
    f[,i] <- as.factor(as.character(f[,i]))
  }
  return(f)
}


################
## vstPCAPlot ##
################
#http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#exploratory-analysis-and-visualization
vstReturnPCAPlot <- function(dds, intgroup, PLOTDIR, TAG, ntop=NULL, plotsOnly=NULL){
  if(class(dds)!="DESeqDataSet"){
    print(paste0(dds, " is not a DESeqDataSet"))
  }
  else{
  if(is.null(ntop)){
    ntop <- 500
  }
  if(is.null(plotsOnly)){
    plotsOnly <- FALSE
  }
  if(!is.numeric(ntop)){
    print("Require a numeric for ntop")
    break
  }
  vstdds <- varianceStabilizingTransformation(dds)
  pcaData <- DESeq2::plotPCA(vstdds, intgroup = intgroup, returnData=T, ntop=ntop)
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
  if(plotsOnly==TRUE){
    return(pcaPlot)
  }
  if(plotsOnly==FALSE){
    ggsave(filename=paste0(PLOTDIR, "/", TAG, ".ntop_", ntop, ".deseq2_pca.pdf"), plot=pcaPlot)
    return(vstdds)
  }
}}

##create ExpressionSet from data.frame, clinical data

createEset <- function(dataframe, pheno){

  ##uses first column as rownames if given a tibble as clinical data
  ##else assumes rownames are sampleIDs, if not make it so
  ##reorders on pData rownames

  if(is_tibble(pheno)){
    pheno <- pheno %>% as.data.frame()
    rownames(pheno) <- pheno[,1]
    pheno <- pheno[,-1]
  }
  pData <- new("AnnotatedDataFrame",
                data=data.frame(pheno))
  datamatrix <- dataframe %>% dplyr::select(rownames(pData)) %>% as.matrix()
  combined.eset <- ExpressionSet(assayData=datamatrix,
                                 phenoData=pData)
  return(combined.eset)
}
