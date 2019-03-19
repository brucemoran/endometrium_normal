##survival for input targets
##inputs:
#NZCOUNTS <- counts in tibble, remove all zero-sum
#GENCOL <- named column in counts with IDs to use in survival

library(survival)
library(tidyverse)
library(biomaRt)
library(DESeq2)

#########
## TPM ##
#########
GENOME <- "hsapiens_gene_ensembl"
convertToTpm <- function(countData, geneCol) {
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = GENOME)
  if(!is.tibble(countData)){
    print("Not a tibble, retry with tibble")
    break
  }
  if(is.tibble(countData)){
    attributes <- c("ensembl_gene_id",
                    "external_gene_name",
                    "start_position",
                    "end_position")
    countData <- countData %>% dplyr::arrange_(geneCol)
    lengthsAll <- biomaRt::getBM(attributes = attributes,
                                 mart = mart) %>%
               dplyr::filter(!!as.name(geneCol) %in% unlist(countData)) %>%
               dplyr::arrange_(geneCol) %>%
               dplyr::mutate(gene_length = end_position - start_position) %>%
               dplyr::select(geneCol, gene_length)
    lengthData <- left_join(countData, lengthsAll) %>%
                   dplyr::select(gene_length)
    countDataNum <- countData %>% select_if(is.numeric)
    rateData <- countDataNum / unlist(lengthData)
    rateData1 <- rateData / sum(rateData) * 1e6
    tpmData <- as_tibble(cbind(countData[,geneCol], rateData1))
    return(tpmData)
  }
}

##convert string to numeric for survival
##0 <- alive (no event)
areYouSurvive <- function(OSUSE){
  ##keywords suggestive of good/bad
  kw0Vec <- c("iv","no","No")
  kw1Vec <- c("yes","Yes","ea")
  ##is osuse[,1] numeric
  if(!is.numeric(as.vector(OSUSE))){
    set0 <- 0;
    set1 <- 0;
    for(x in 1:3){
      if(length(grep(kw0Vec[x],as.vector(OSUSE[1])))>0) set0 <- x
    }
    if(set0==0){
      for(x in 1:3){
        if(length(grep(kw1Vec[x],as.vector(OSUSE[1])))>0) set1 <- x
      }
    }
    if(sum(set1,set0)==0){
      return(0)
    }
    if(set0!=0){
      gsuse <- as.numeric(! as.vector(OSUSE) %in% as.vector(OSUSE[1]))
      return(gsuse)
    }
    if(set1!=0){
      gsuse <- as.numeric(as.vector(OSUSE) %in% as.vector(OSUSE[1]))
      return(gsuse)
    }
  }
}

#######################
## Survival Function ##
#######################
survivalAnalysis <- function(GENES, CLINICAL, CLINCOLSURV, CLINCOLDAYS, CLINSAMPLECOL, EXPR, EXPRGENECOL, TAG){
  ##genes of interest:
  clin <- read_tsv(CLINICAL)
  osCol <- CLINCOLSURV
  daysCol <- CLINCOLDAYS

  clinsamples <- clin %>% dplyr::select_(CLINSAMPLECOL) %>% c() %>% unlist()
  exprsamples <- colnames(EXPR)[colnames(EXPR) %in% clinsamples]
  clin <- clin %>% dplyr::filter(clinsamples %in% exprsamples) %>%
                   arrange_(CLINSAMPLECOL)
  os <- clin %>% dplyr::select(!!dssCol) %>% c() %>% unlist()
  os <- areYouSurvive(os)
  days <- clin %>% dplyr::select(!!daysCol) %>% unlist() %>% as.integer()
  days[is.na(days)] <- 0
  ##normalise to library by DESeq2
  EXPRCLINS <- EXPR %>%
          dplyr::select(EXPRGENECOL, exprsamples)
  ddsExpr <- as.data.frame(EXPRCLINS[,-grep(EXPRGENECOL,colnames(EXPRCLINS))])
  rownames.ddsExpr <- unlist(EXPRCLINS[,EXPRGENECOL])
  ddsClin <- as.data.frame(clin)
  rownames(ddsClin) <- unlist(clin[,CLINSAMPLECOL])
  dds <- DESeqDataSetFromMatrix(countData = ddsExpr,
                                colData = ddsClin,
                                design =~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)

  EXPRNORM <- tibble::as_tibble(counts(dds, normalized = TRUE)) %>%
              dplyr::mutate(rownames.ddsExpr) %>%
              dplyr::select(rownames.ddsExpr, exprsamples)
  colnames(EXPRNORM) <- c(EXPRGENECOL, exprsamples)
  EXPRGENESNORM <- EXPRNORM %>% dplyr::filter(!!as.name(EXPRGENECOL) %in% GENES)

  EXPRGENESTPM <- convertToTpm(EXPRGENESNORM, EXPRGENECOL)
  tpmexpr <- as.data.frame(EXPRGENESTPM)
  rownames(tpmexpr) <- unlist(EXPRGENESTPM[,EXPRGENECOL])
  tpmexpr <- tpmexpr[,-grep(EXPRGENECOL, colnames(tpmexpr))]

  for (GENEX in 1:length(GENES)){

    ##set dir, load data
    wd <- dir.create(paste0(getwd(),"/survival_",GENES[GENEX]), showWarnings=F)
    setwd(wd)

    ###############################
    ## Survival Analysis Running ##
    ###############################
    ##survival function iteratively over the above sets to make plots:
    gene <- GENES[GENEX]
    osuse <- data.frame(row.names=colnames(vstexpr1), days, os)

    geneExpr <- unlist(vstexpr1[gene,])
    medianExpr[as.numeric(as.vector(geneExpr)) >= median(as.numeric(as.vector(geneExpr)))] <- 1

    table(medianExpr)

    ##include age, stage
    survObj <- Surv(time=osuse[,1], event=osuse[,2])
    survFit <- survfit(survObj ~ medianExpr)
    svd <- survdiff(survObj ~ medianExpr)
    coxPH <- coxph(survObj ~ medianExpr)
    capt <- capture.output(survdiff(formula = survObj ~ medianExpr))
    pval <- strsplit(capt[length(capt)],",")[[1]][2]

    pdf(paste0("TCGA_UCEC","-",gene,"-median_survival.pdf"))
    plot(survFit, lty = 1:4, xlab = "Days To Event", ylab = "Proportion Surviving",
    lwd = 2,col=c("forestgreen","red"),main=paste0("TCGA_UCEC - ", gene," median survival"))
    legend(x = 6, y = 1, legend = c(paste0(gene," below median (n=",svd$n[1],")"),paste0(gene," above median (n=",svd$n[2],")")), lty = c(1:4), lwd = 2, cex = 0.8,col=c("forestgreen","red"))
    text(8, 0.7, pval, pos = 4)
    dev.off()
  }
}
