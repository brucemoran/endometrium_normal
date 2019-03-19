##R script to run GOseq, CAMERA, ROAST based on:
##https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-gene-set-testing.nb.html
options(scipen=999)
options(stringAsFactors=FALSE)

#libraries
libs <- c("edgeR", "fgsea", "CAMERA", "tidyverse", "dplyr", "biomaRt", "DEFormats")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

##inputs:
##res = tibble of topTags, n=inf format, colname FDR
##geneCol = colname for gene IDs, NB must be Ensembl ENS... format
##sigp = level for significance testing
##DGE, DESIGN, CONTRAST all from limma
##MSIGDIR = location of MSigDB pathways
runCAMERA <- function(RES, GENECOL, DGE, DESIGN, CONTRAST, MSIGDIR, GENOME=NULL, SIGP=NULL, MSIGMATCH=NULL){

  ##test for Ensembl gene ID
  if(!length(grep("^ENS",c(unlist(RES[1, GENECOL]))))==1){
    print("Gene IDs must be from Ensembl")
    break
  }
  ##suppose human if not stated otherwise, else use to convert to human
  if(is.null(GENOME)){
    GENOME <- "hsapiens_gene_ensembl"
  }
  ##significance level arbitrary default
  if(is.null(SIGP)){
    SIGP <- 0.05
  }
  ##MSigDB pathways matched to this vector are used in analysis
  if(is.null(MSIGMATCH)){
    MSIGMATCH <- c(".")
  }

  #############
   ## CAMERA #
    #########
  pathways <- unique(unlist(lapply(MSIGMATCH,function(f){grep(f,dir(pattern="gmt", MSIGDIR), value=T)})))
  pathwaysCameraResList <- lapply(pathways, function(pway){
    pathwayList <- gmtPathways(paste0(MSIGDIR, "/", pway))
    print(paste0("Working on: ", pway, "; length: ", length(pathwayList)))

    ##if genome is not human
    if(GENOME != "hsapiens_gene_ensembl"){
      mart <- biomaRt::useMart(biomart = "ensembl", dataset = GENOME)
      gene2hgene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                          "hsapiens_homolog_ensembl_gene"),
                                             mart = mart)) %>%
                    dplyr::rename(human_ensembl_gene_id="hsapiens_homolog_ensembl_gene")
      humart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      gene2hgene2hname <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                                "external_gene_name"),
                                                   mart = humart)) %>%
                          dplyr::rename(human_ensembl_gene_id="ensembl_gene_id") %>%
                          left_join(gene2hgene, .) %>%
                          na.omit()
    }
    ##if human already, find external_gene_name and rename
    if(GENOME == "hsapiens_gene_ensembl"){
      humart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      gene2hgene2hname <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                                "external_gene_name"),
                                                   mart = humart)) %>%
                          na.omit()
    }

    ##now have two tibbles, each with ensembl_id (any) mapped to external_gene_id (human)
    ##this is used to prune the pathwayList
    prunedEnsPathwayList <- lapply(pathwayList, function(f){
      gene2hgene2hname %>% dplyr::filter(external_gene_name %in% f) %>%
                                  dplyr::pull(ensembl_gene_id) %>% unique() %>% sort()
    })
    prunedExtPathwayList <- lapply(pathwayList, function(f){
      gene2hgene2hname %>% dplyr::filter(external_gene_name %in% f) %>%
                                dplyr::pull(external_gene_name) %>% unique() %>% sort()
    })
    ##index
    indexd <- ids2indices(gene.sets=prunedEnsPathwayList, identifiers=rownames(DGE$counts))
    if(length(indexd)==0){
      indexd <- ids2indices(gene.sets=prunedExtPathwayList, identifiers=rownames(DGE$counts))
    }
    ##estimate dispersion, even if estimated
    DGED <- estimateDisp(DGE, design=DESIGN)
    ##attach genes from the pathway for further investigation
    camera.res <- as_tibble(camera.DGEList(DGED,
                                           index=indexd,
                                           design=DESIGN,
                                           contrast=CONTRAST,
                                           inter.gene.cor=0.05),
                                           rownames="MSIGDB_ID") %>%
                  dplyr::filter(FDR < SIGP)

    if(length(pull(camera.res, MSIGDB_ID))>0){
        camera.res$MSIGDB_GENESET <- unlist(lapply(pull(camera.res, MSIGDB_ID), function(f){
          paste(unlist(prunedExtPathwayList[[f]]),collapse=",")
        }))
        camera.res$MSIGDB_LINK <- paste0("http://software.broadinstitute.org/gsea/msigdb/cards/", pull(camera.res, MSIGDB_ID))
    }

    ##attach
    return(camera.res)
  })
  names(pathwaysCameraResList) <- pathways
  return(pathwaysCameraResList)
}

####################
## DESeq2 CAMERA ##
#################

##inputs:
##res = tibble of DESeq2 results
##geneCol = colname for gene IDs in RES, NB must be Ensembl ENSx format
##sigp = level for significance testing
##DDS = DESeq2 Data Set; DESIGN, CONTRAST all converted to limma format with DEFormats
##MSIGDIR = location of MSigDB pathways
runDESeq2CAMERA <- function(RES, GENECOL, DDS, DESIGN, CONTRAST, MSIGDIR, GENOME=NULL, SIGP=NULL, MSIGMATCH=NULL){

  ##test for Ensembl gene ID
  if(!length(grep("^ENS",c(unlist(RES[1, GENECOL]))))==1){
    print("Gene IDs must be from Ensembl")
    break
  }
  ##suppose human if not stated otherwise, else use to convert to human
  if(is.null(GENOME)){
    GENOME <- "hsapiens_gene_ensembl"
  }
  ##significance level arbitrary default
  if(is.null(SIGP)){
    SIGP <- 0.05
  }
  ##MSigDB pathways matched to this vector are used in analysis
  if(is.null(MSIGMATCH)){
    MSIGMATCH <- c(".")
  }

  ##convert DDS to DGE
  DGE <- as.DGEList(DDS)

  #############
   ## CAMERA #
    #########
  pathways <- unique(unlist(lapply(MSIGMATCH,function(f){grep(f,dir(pattern="gmt", MSIGDIR), value=T)})))
  pathwaysCameraResList <- lapply(pathways, function(pway){
    pathwayList <- gmtPathways(paste0(MSIGDIR, "/", pway))
    print(paste0("Working on: ", pway, "; length: ", length(pathwayList)))

    ##if genome is not human
    if(GENOME != "hsapiens_gene_ensembl"){
      mart <- biomaRt::useMart(biomart = "ensembl", dataset = GENOME)
      gene2hgene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                          "hsapiens_homolog_ensembl_gene"),
                                             mart = mart)) %>%
                    dplyr::rename(human_ensembl_gene_id="hsapiens_homolog_ensembl_gene")
      humart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      gene2hgene2hname <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                                "external_gene_name"),
                                                   mart = humart)) %>%
                          dplyr::rename(human_ensembl_gene_id="ensembl_gene_id") %>%
                          left_join(gene2hgene, .) %>%
                          na.omit()
    }
    ##if human already, find external_gene_name and rename
    if(GENOME == "hsapiens_gene_ensembl"){
      humart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
      gene2hgene2hname <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",
                                                                "external_gene_name"),
                                                   mart = humart)) %>%
                          na.omit()
    }

    ##now have two tibbles, each with ensembl_id (any) mapped to external_gene_id (human)
    ##this is used to prune the pathwayList
    prunedEnsPathwayList <- lapply(pathwayList, function(f){
      gene2hgene2hname %>% dplyr::filter(external_gene_name %in% f) %>%
                                  dplyr::pull(ensembl_gene_id) %>% unique() %>% sort()
    })
    prunedExtPathwayList <- lapply(pathwayList, function(f){
      gene2hgene2hname %>% dplyr::filter(external_gene_name %in% f) %>%
                                dplyr::pull(external_gene_name) %>% unique() %>% sort()
    })
    ##index
    indexd <- ids2indices(prunedEnsPathwayList, rownames(DGE$counts))
    ##estimate dispersion, even if estimated
    DGED <- estimateDisp(DGE, design=DESIGN)
    ##attach genes from the pathway for further investigation
    camera.res <- as_tibble(camera.DGEList(DGED,
                                           index=indexd,
                                           design=DESIGN,
                                           contrast=CONTRAST,
                                           inter.gene.cor=0.05),
                                           rownames="MSIGDB_ID") %>%
                  dplyr::filter(FDR < SIGP)

    if(length(pull(camera.res, MSIGDB_ID))>0){
        camera.res$MSIGDB_GENESET <- unlist(lapply(pull(camera.res, MSIGDB_ID), function(f){
          paste(unlist(prunedExtPathwayList[[f]]),collapse=",")
        }))
        camera.res$MSIGDB_LINK <- paste0("http://software.broadinstitute.org/gsea/msigdb/cards/", pull(camera.res, MSIGDB_ID))
    }

    ##attach
    return(camera.res)
  })
  names(pathwaysCameraResList) <- pathways
  return(pathwaysCameraResList)
}
