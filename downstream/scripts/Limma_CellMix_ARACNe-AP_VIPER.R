#! /usr/bin/R

##create DESeqDataSet objects, from which TPM are made per gene per sample
##use TPMs for deconRNAseq with early, mid scretory groups to classify DC_normal, SRR_normal, tumour
##based on Filipovic 2018, Molecular definition of group 1 innate lymphoid cells in the mouse uterus
##also run for TPMs in Kallisto output at transcript level
options(scipen=999)
argsIn <- commandArgs(trailingOnly=TRUE)

##libraries
print("Loading libraries...")
libs <- c("httr", "limma", "CAMERA", "DEFormats", "apeglm", "ggplot2", "tidyverse", "plyr", "dplyr", "biomaRt", "sleuth", "CellMix", "reshape2", "roots", "viper", "mixtools", "bcellViper", "aracne.networks", "Biobase", "rgexf", "fgsea", "gtools", "tximport")
library("BiocManager")
libsLoaded <- lapply(libs,function(l){
  if(!l %in% installed.packages()){
    BiocManager::install(l)
  }
  suppressWarnings(suppressMessages(library(l, character.only = T)))
})

##input parameters, load clinical and functions
BASEDIR <- argsIn[1]
CLINDIR <- argsIn[2]
CGCTABLE <- argsIn[3] ##full path of CGC table tsv from website (requires login)
DATADIR <- paste0(BASEDIR, "/analysis/RNAseq")
SCRIPTDIR <-  paste0(BASEDIR, "/scripts")
# TX2GENE <- paste0(CLINDIR, "/Homo_sapiens.GRCh37.cdna.all.fa.tx2gene")
OUTDIR <- paste0(BASEDIR, "/analysis/RNAseq/Limma_CellMix_ARACNe")
RDATDIR <- paste0(OUTDIR,"/RData")
dir.create(RDATDIR, showWarnings=FALSE, recursive=TRUE)

##source functions
print("Loading clinical data and sourcing functions...")
funcs <- dir(pattern="func.R$", SCRIPTDIR, full.names=T)
for (x in 1:length(funcs)){
  print(paste0("Sourcing: ", funcs[x])); source(funcs[x])
}

##load clinical data
ifelse(length(grep("MMMT_TCGA-UCEC-UCS_SRP104165_SRP105769_DCRNASEQ.clin.RData", dir(CLINDIR)))==1,
  load(paste0(CLINDIR, "/MMMT_TCGA-UCEC-UCS_SRP104165_SRP105769_DCRNASEQ.clin.RData")),
  system(paste0("Rscript --vanilla ", SCRIPTDIR, "/clinical_data.R ", CLINDIR))
  load(paste0(CLINDIR, "/MMMT_TCGA-UCEC-UCS_SRP104165_SRP105769_DCRNASEQ.clin.RData"))
)

##annotation
##tx2gene
GENOME <- "hsapiens_gene_ensembl"
mart <- biomaRt::useMart(biomart = "ensembl", dataset = GENOME)
tx2gene <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
gene2name <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = mart))
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

##sample, count inputs
if(length(grep("ddsPCAList.RData",dir(RDATDIR)))==0){
  combined_clinIn <- combined_clin %>% dplyr::filter(sampleID %in% dir(DATADIR))
  sampleIDs <- unlist(combined_clin$sampleID)
  sampleIDsIn <- unlist(combined_clinIn$sampleID)

  ##tximport
  files <- file.path(DATADIR, sampleIDs, "kallisto", "hg19", "abundance.h5")
  names(files) <- sampleIDs
  files <- files[file.exists(files)==TRUE]
  sampleIDuse <- names(files)
  txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = TRUE)
  counts <- txi.kallisto.tsv$counts
  tpms <- txi.kallisto.tsv$abundance
  tpmstsv <- as_tibble(tpms, rownames="ensembl_gene_id")
  write_tsv(tpmstsv, path=paste0(OUTDIR, "/kallisto.gene.tpms.tsv"))

  ##Reduce to tibble, removing all rows summing to zero (tidyverse <3)
  print("Annotating and filtering count data...")
  nzCountsAnnoTr <- left_join(annoBM(attrs=c("ensembl_gene_id", "external_gene_name"),
                                 filts="ensembl_gene_id",
                                 arrone="ensembl_gene_id",
                                 mart=mart),
                              as_tibble(counts, rownames="ensembl_gene_id")) %>%
                     dplyr::mutate(Mean = apply(dplyr::select(.,unlist(sampleIDuse)),1,mean)) %>%
                     dplyr::filter(Mean > 5) %>%
                     dplyr::select(-Mean) %>%
                     arrange(ensembl_gene_id) %>%
                     dplyr::mutate_if(is.numeric, round, 0) %>%
                     na.omit()

  ##Normal and TCGA data into dds objects, SF-normalise, PCA plots, VST outputs
  ddsPCAList <- ddsPCA(nzCountsAnnoTr, OUTDIR, topn=500)
  ddsList <- list(ddsPCAList[[5]], ddsPCAList[[6]], ddsPCAList[[7]])
  ddsPCAList <- list(ddsPCAList[[1]], ddsPCAList[[2]], ddsPCAList[[3]], ddsPCAList[[4]])
  save(ddsPCAList, nzCountsAnnoTr, combined_clin, file=paste0(RDATDIR,"/ddsPCAList.RData"))
  save(ddsList, nzCountsAnnoTr, combined_clin, file=paste0(RDATDIR,"/ddsList.RData"))
}
if(length(grep("ddsPCAList.RData",dir(RDATDIR)))==1){
  load(paste0(RDATDIR,"/ddsPCAList.RData"))
}

nz.co.all.vst <- ddsPCAList[[1]]
nz.co.norm40.vst <- ddsPCAList[[2]]
cond.norm40 <- ddsPCAList[[3]]
nz.co.norm40 <- ddsPCAList[[4]]

##download and preprocess GSE97929 data
if(length(grep("GSE97929_preprocess.RData",dir(RDATDIR)))==0){
  ss.allreads.vst.med <- GSE97929preprocess(GSE97929_clin,
                                            nzCountsAnnoTr,
                                            OUTDIR=OUTDIR)
  save(ss.allreads.vst.med, file=paste0(RDATDIR, "/GSE97929_preprocess.RData"))
}
if(length(grep("GSE97929_preprocess.RData",dir(RDATDIR)))==1){
  load(paste0(RDATDIR,"/GSE97929_preprocess.RData"))
}

##CellMix to determine fractions and 'absolute' of cellphase
if(length(grep("cellMixList.RData",dir(RDATDIR)))==0){
  cellMixList <- suppressWarnings(cellMixplot(nz.co.all.vst,
                                              nz.co.norm40.vst,
                                              ss.allreads.vst.med,
                                              combined_clin,
                                              BASEDIR=BASEDIR,
                                              OUTDIR=OUTDIR))
  save(cellMixList, file=paste0(RDATDIR, "/cellMixList.RData"))
}
if(length(grep("cellMixList.RData",dir(RDATDIR)))==1){
  load(paste0(RDATDIR,"/cellMixList.RData"))
}

datasets <- cellMixList[[1]]
signatures <- cellMixList[[2]]
combined_clin_res <- cellMixList[[3]]

##DESeq2 on cellphases to determine DE genes, run GSEA
##NB creates plots, files, nothing used ongoing herein
if(length(grep("cellphaseLimmaCAMERA.RData", dir(RDATDIR)))==0){
  cellphaseLimmaCAMERA(combined_clin_res,
                       nz.co.norm40,
                       gene2name,
                       OUTDIR=OUTDIR,
                       PATHWAYDIR=paste0(BASEDIR,"/data/pathway"))
}

##run ARACNE-AP in PBS, and msVIPER with those results back in R
##NB once ARACNE-AP runs, those bootstraps are reused unless dir deleted
datasetsUCS <- datasets[,colnames(datasets) %in% c(levels(cond.norm40$sampleID),
                                                 c(NIHMS891671_clin$bcr_patient_barcode))]
datasetsUCEC <- datasets[,colnames(datasets) %in% c(levels(cond.norm40$sampleID),
                                                   c(NATURE12113_clin$bcr_patient_barcode))]
combined_clin_res.epi <- combined_clin_res %>%
                           dplyr::filter(sampleID %in% colnames(datasetsUCEC)) %>%
                           dplyr::filter(whichmax %in% grep("epi",whichmax,value=TRUE))
combined_clin_res.str <- combined_clin_res %>%
                           dplyr::filter(sampleID %in% colnames(datasetsUCEC)) %>%
                           dplyr::filter(whichmax %in% grep("str",whichmax,value=TRUE))
datasetsUCEC.epi <- datasetsUCEC[,colnames(datasetsUCEC) %in% combined_clin_res.epi$sampleID,]
datasetsUCEC.str <- datasetsUCEC[,colnames(datasetsUCEC) %in% combined_clin_res.str$sampleID,]

##inputs to msViper iteration
msviperTags <- c("norm40_vs_TCGA-UCS",
                 "norm40_vs_TCGA-UCEC",
                 "norm40_vs_TCGA-UCEC_epi",
                 "norm40_vs_TCGA-UCEC_str")

msviperDatasets <- list(datasetsUCS, datasetsUCEC, datasetsUCEC.epi, datasetsUCEC.str)
msviperClinres <- list(combined_clin_res, combined_clin_res, combined_clin_res.epi, combined_clin_res.str)
msviperList <- as.list(1:4)
names(msviperList) <- names(msviperDatasets) <- names(msviperClinres) <- msviperTags

##iterate
for(x in 1:length(msviperTags)){

  if(length(grep(paste0("msViper.", msviperTags[x], ".RData"),dir(RDATDIR)))==0){
    print(paste0("Working on: ", msviperTags[x] , "..."))
    dataset <- msviperDatasets[[msviperTags[x]]]
    runAracneAP(dataset, signatures, msviperTags[x], BASEDIR, OUTDIR)
    msviperList[[msviperTags[x]]] <- runMsViper(combined_clin_res=msviperClinres[[msviperTags[x]]],
                            TAG=msviperTags[x],
                            BASEDIR=BASEDIR,
                            OUTDIR=OUTDIR,
                            DERES=paste0(OUTDIR, "/DESeq2_FGSEA/res.str_early.vs.mid.nz.txt"))
  }
  else{
    print(paste0("Found completed: ", msviperTags[x] , "..."))
  }
}
save(msviperList, file=paste0(RDATDIR, "/msViper.all.RData"))

##final output
print(paste0("Saving full output RData..."))
save.image(paste0(RDATDIR, "/DESeq2_CellMix_ARACNe-AP_VIPER.RData"))
