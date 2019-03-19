#! /apps/software/R/3.4.0/bin/R

##collate all clinical data into single file
#files are:
## TCGA UCS data: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5599133/bin/NIHMS891671-supplement-2.xlsx from DOI:10.1016/j.ccell.2017.02.010
## GSE98386_series_matrix.txt.gz from ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE98nnn/GSE98386/matrix/
argsIn <- commandArgs(trailingOnly=TRUE)

libs <- c("httr", "GEOquery", "tidyverse", "readxl" ,"downloader")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = T)))})

##inputs
CLINDIR <- argsIn[1]
setwd(CLINDIR)

##TCGA UCEC clinical data
if(length(grep("NATURE12113.csv", dir(CLINDIR)))==0){
  ucec_url <- "https://media.nature.com/original/nature-assets/nature/journal/v497/n7447/extref/nature12113-s2.zip"
  download(ucec_url, dest="/tmp/tmp.zip", mode="wb")
  unzip("/tmp/tmp.zip", exdir=CLINDIR)
  ucec_clin <- read_excel(paste0(CLINDIR,"/datafile.S1.1.KeyClinicalData.xls"))
  write_csv(ucec_clin, "NATURE12113.csv")
}

##TCGA UCS clinical data
if(length(grep("NIHMS891671.csv", dir(CLINDIR)))==0){
  ucs_url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5599133/bin/NIHMS891671-supplement-2.xlsx"
  GET(ucs_url, write_disk(tf <- tempfile("")))
  ucs_clin <- read_excel(tf, skip=1)
  unlink(tf)
  write_csv(ucs_clin, "NIHMS891671.csv")
}

##map SRRs to GSEs (projects), then return map of all GSMs (samples) for each SRR (also samples)
##all SRRs
if(length(grep("series.matrix.csv", dir(CLINDIR)))==0){

  SRRs <- grep("SRR",dir("../fastq/"), value=TRUE)

  ##get current map file (~4GB)
  sra_gsm_map_ftp <- "ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab"
  GET(sra_gsm_map_ftp, write_disk(tf <- tempfile("")))
  sra_gsm_map <- read_tsv(tf)

  ##GSEs of interest
  names(SRP_GSEs) <- SRP_GSEs <- c("SRP104165", "SRP105769")
  names(GSEs) <- GSEs <- c("GSE97929", "GSE98386")

  ##SRRs of GSEs, and their GSMs
  SRR_GSMs <- sra_gsm_map %>%
                dplyr::filter(Study %in% SRP_GSEs) %>%
                dplyr::select(Accession, Alias, Study) %>%
                separate(., Alias, "geo_accession", "_", extra="drop")
  SRR_GSMs <- SRR_GSMs[-starts_with("SRX",vars=SRR_GSMs$Accession),]

  unlink(tf)

  ##GEO clincal data for all
  GSE_clin <- lapply(GSEs, function(g){
    gse <- getGEO(g)
    clin <- as_tibble(gse[[1]])
    ##remap names
    clin <- inner_join(SRR_GSMs, clin)
    return(clin)
  })

  for (x in 1:length(names(GSE_clin))){
    write_csv(GSE_clin[[x]], paste0(names(GSE_clin)[x], ".series.matrix.csv"))
  }
}

###########################
## MANUAL REVIEW OF DATA ##
###########################
##source("clin_data_cols.R")
GSE97929_use <- c("Accession",
                  "geo_accession",
                  "Study",
                  "title",
                  "menstrual.cycle.phase.ch1",
                  "library.ch1",
                  "tissue.ch1",
                  "country.of.origin.ch1")
GSE98386_use <- c("Accession",
                  "geo_accession",
                  "Study",
                  "title",
                  "receptivity.time.point.ch1",
                  "sequencing.batch.ch1",
                  "tissue.ch1",
                  "tissue.state.ch1",
                  "country.of.origin.ch1")
NIHMS891671_use <- c("bcr_patient_barcode",
                     "clinical_stage",
                     "primary_therapy_outcome_success",
                     "new_tumor_event_after_initial_treatment",
                     "days_to_new_tumor_event_after_initial_treatment",
                     "person_neoplasm_cancer_status",
                     "vital_status",
                     "DiseaseStatus",
                     "InferredTumorStatus",
                     "Prog/Recurrence",
                     "PFS",
                     "OS",
                     "days_to_death",
                     "days_to_last_followup",
                     "race",
                     "ethnicity",
                     "histological_type",
                     "Histologic classification")
NATURE12113_use <- c("bcr_patient_barcode",
                     "patient_type",
                     "age",
                     "2009stagegroup",
                     "vital_status",
                     "disease_status_at_lfu",
                     "recurred_progressed",
                     "pfs_days",
                     "os_days",
                     "race",
                     "ethnicity",
                     "BMI",
                     "histology",
                     "histology_grade",
                     "IntegrativeCluster",
                     "msi_status_7_marker_call",
                     "cna_cluster_k4",
                     "mrna_expression_cluster")

dcrosby_use <- c("Sample", "Group", "Batch", "Embryo")
dcrosby_clin_use <- dcrosby_clin <- read_csv("dcrosby_pregnancy_pheno.csv")
GSE97929_clin <- read_csv("GSE97929.series.matrix.csv")
GSE98386_clin <- read_csv("GSE98386.series.matrix.csv")
NIHMS891671_clin <- read_csv("NIHMS891671.csv")
NATURE12113_clin <- read_csv("NATURE12113.csv")
GSE97929_clin_use <- GSE97929_clin %>% dplyr::select(one_of(GSE97929_use))
GSE98386_clin_use <- GSE98386_clin %>% dplyr::select(one_of(GSE98386_use))
NIHMS891671_clin_use <- NIHMS891671_clin %>% dplyr::select(one_of(NIHMS891671_use))
NATURE12113_clin_use <- NATURE12113_clin %>% dplyr::select(one_of(NATURE12113_use))

##combine
combined_clin <- dcrosby_clin %>%
                 dplyr::mutate(sampleID = Sample,
                               Tissue = rep("Endometrium", times=dim(dcrosby_clin)[1]),
                               Strand=rep(2, times=dim(dcrosby_clin)[1]),
                               Source="UCD_DC") %>%
                 dplyr::select(sampleID, Group, Tissue, Batch, Strand, Source)

combined_clin <- rbind(combined_clin, GSE97929_clin %>%
                 dplyr::mutate(sampleID = Accession,
                        Group = menstrual.cycle.phase.ch1,
                        Tissue = tissue.ch1,
                        Batch = library.ch1 + max(combined_clin$Batch),
                        Strand = rep(3, times=dim(GSE97929_clin)[1]),
                        Source="GSE97929") %>%
                 dplyr::select(sampleID, Group, Tissue, Batch, Strand, Source))

combined_clin <- rbind(combined_clin, GSE98386_clin %>%
                 dplyr::mutate(sampleID = Accession,
                        Group = unlist(lapply(receptivity.time.point.ch1,function(f){
                          if(f == "LH+2"){"early secretory"}
                          else{"mid secretory"}})),
                        Tissue = tissue.ch1,
                        Batch = sequencing.batch.ch1 + max(combined_clin$Batch),
                        Strand=rep(2, times=dim(GSE98386_clin)[1]),
                        Source="GSE98386") %>%
                 dplyr::select(sampleID, Group, Tissue, Batch, Strand, Source))

##NB on Batch, bioinformatics.mdanderson.org/tcgambatch/ gives DSC 0, so batch effect within UCS is minimal
combined_clin <- rbind(combined_clin, NIHMS891671_clin %>%
                 dplyr::mutate(sampleID = bcr_patient_barcode,
                        Group = clinical_stage,
                        Tissue = unlist(lapply(NIHMS891671_clin$histological_type,function(f){
                          if(length(grep("Het",f))==1){return("HET")}
                          if(length(grep("Hom", f))==1){return("HOM")}
                          if(length(grep("NOS", f))==1){return("NOS")}})),
                        Batch = 1 + max(combined_clin$Batch),
                        Strand=rep(0, times=dim(NIHMS891671_clin)[1]),
                        Source="TCGA_UCS") %>%
                 dplyr::filter(race == "WHITE", ethnicity != "HISPANIC OR LATINO") %>%
                 dplyr::select(sampleID, Group, Tissue, Batch, Strand, Source))

##NB on Batch, bioinformatics.mdanderson.org/tcgambatch/ gives DSC 0, so batch effect within UCS is minimal
combined_clin <- rbind(combined_clin, NATURE12113_clin %>%
                 dplyr::mutate(sampleID = bcr_patient_barcode,
                        Group = `2009stagegroup`,
                        Tissue = NATURE12113_clin$histology,
                        Batch = 1 + max(combined_clin$Batch),
                        Strand=rep(0, times=dim(NATURE12113_clin)[1]),
                        Source="TCGA_UCEC") %>%
                 dplyr::filter(race == "WHITE", ethnicity != "HISPANIC OR LATINO") %>%
                 dplyr::select(sampleID, Group, Tissue, Batch, Strand, Source))

##MMMT
NIHMS891671_clin <- NIHMS891671_clin %>%
  dplyr::mutate(MMMTtype = unlist(lapply(histological_type, function(f){
    strsp <- strsplit(f,": ")[[1]][2]
    strsp <- gsub("Heterologous Type", "HET", strsp)
    strsp <- gsub("Homologous Type", "HOM", strsp)
    return(strsp)
  })))

##UCEC retain all
save(dcrosby_clin, GSE97929_clin, GSE98386_clin, NIHMS891671_clin, NATURE12113_clin, dcrosby_clin_use, GSE97929_clin_use, GSE98386_clin_use, NIHMS891671_clin_use, combined_clin, file="MMMT_TCGA-UCEC-UCS_SRP104165_SRP105769_DCRNASEQ.clin.RData")
