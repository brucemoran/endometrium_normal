#! R

##preprocesses cell-type data from GSE97929
##clinical data already preprocessed
##geneVec holds genes to be used in the analysis, e.g. from

GSE97929preprocess <- function(GSE97929_clin){

  #############
  ## CellMix ##
  #############
  print("Reading in cell-type counts...")
  allreads_https  <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97929&format=file&file=GSE97929%5Freads%5Fall%2Etxt%2Egz"

  httr::GET(allreads_https, httr::write_disk(tf <- tempfile("")))
  allreads <- as_tibble(read.table(tf, header=T)) %>% dplyr::filter(as.vector(Gene) %in% grep("RNA_SPIKE", Gene, value=TRUE, invert=TRUE))
  unlink(tf)
  ss.all <- allreads
  colnames(ss.all) <- c("external_gene_name",
                           as.character(unlist(lapply(colnames(ss.all)[2:dim(ss.all)[2]], function(f){
                             strsp <- strsplit(f,"\\.")[[1]][3]
                             strsp <- gsub("_1", "_2", strsp)
                             strsp <- paste0(strsp,"_1")
                             strsp <- gsub("_2_1", "_2", strsp)
                             return(strsp)
                           }))))
  #remove zero-sum
  ss.all <- ss.all %>%
            dplyr::mutate(Median = rowMedians(as.matrix(.[grep("exter",names(.), invert=TRUE)]))) %>%
            dplyr::filter(Median > 0) %>%
            dplyr::select(-Median) %>%
            na.omit() %>%
            dplyr::select(-"NA_1") %>%
            column_to_rownames(var="external_gene_name") %>%
            as.data.frame()

  #condition
  ss.cond <- GSE97929_clin %>%
             dplyr::mutate(sampleID = unlist(lapply(title, function(f){
                      paste0(strsplit(f,"_")[[1]][1], "_", strsplit(f,"_")[[1]][2])
                    }))) %>%
             dplyr::mutate(Group_1 = unlist(lapply(characteristics_ch1, function(f){
                      paste0(strsplit(f," ")[[1]][4])
                    }))) %>%
             dplyr::mutate(Group_2 = unlist(lapply(characteristics_ch1.3, function(f){
                      paste0(strsplit(f," ")[[1]][3])
                    }))) %>%
             dplyr::mutate(Group = paste0(Group_1,"_",Group_2)) %>%
             dplyr::select(sampleID, library.ch1, Group) %>%
             dplyr::rename(Batch = library.ch1)
  ss.all <- ss.all %>% dplyr::select(ss.cond$sampleID)
  ss.cond <- ss.cond %>% dplyr::filter(sampleID %in% colnames(ss.all))
  dds.ss.all <- DESeqDataSetFromMatrix(countData=ss.all ,
                                       colData=ss.cond,
                                       design=~1)
  dds.ss.all <- estimateSizeFactors(dds.ss.all)
  ss.all.vst <- varianceStabilizingTransformation(dds.ss.all)

  ##deconvolute mid-secretory samples
  filterTissuePhase <- function(gse, tissue, phase){
    gse %>% dplyr::filter(characteristics_ch1.3 == !!tissue) %>%
            dplyr::filter(characteristics_ch1 == !!phase) %>%
            dplyr::transmute(sampleID=unlist(lapply(title, function(f){
              strsp <- strsplit(f, "_")[[1]]
              stro <- paste0(strsp[1], "_", strsp[2])
              return(stro)
  }))) %>% unique() %>% c()}

  ss.all.tissuephase <- c(early_epithelium =
                        filterTissuePhase(GSE97929_clin,
                                          "tissue: endometrium epithelium",
                                          "menstrual cycle phase: early secretory"),
                         mid_epithelium = filterTissuePhase(GSE97929_clin,
                                          "tissue: endometrium epithelium",
                                          "menstrual cycle phase: mid secretory"),
                         early_stroma = filterTissuePhase(GSE97929_clin,
                                          "tissue: endometrium stroma",
                                          "menstrual cycle phase: early secretory"),
                         mid_stroma = filterTissuePhase(GSE97929_clin,
                                          "tissue: endometrium stroma",
                                          "menstrual cycle phase: mid secretory"))

  #which of those to use downstream?
  ss.all.chosen <- ss.all.tissuephase

  # ##make TPM
  # ss.tpmreads <- convertToTpm(countData=ddss.allreads.esf, geneCol="external_gene_name")
  ##previously made TPM, now using VST as per:
  ##http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization

  ##previoulsy took input for GSE98386, to remove genes not therein
  ##this filtered at a later stage now
  # ##set for each tissue/tissuephase also in GSE98386
  # gene.set.use.nzc <- geneVec
  # ss.all.vst.topk <- as_tibble(assay(ss.all.vst), rownames="external_gene_name") %>%
  #                        dplyr::filter(external_gene_name %in% gene.set.use.nzc)

  #take median of above splits
  ss.all.vst.list.med <- lapply(seq_along(ss.all.chosen), function(f){
              renam <- strsplit(names(ss.all.chosen)[f],"\\.")[[1]][1]

              ##to use  ss.all.vst.topk, change below ss.all.vst to ss.all.vst.topk
              as_tibble(assay(ss.all.vst), rownames="external_gene_name") %>% dplyr::select(external_gene_name, !!ss.all.chosen[[f]]) %>%
                             dplyr::mutate(Median = rowMedians(as.matrix(.[grep("exter",names(.),  invert=TRUE)]))) %>%
                             dplyr::select(external_gene_name, Median) %>%
                             dplyr::rename(!!renam := Median)
            })

  ##remove rows with sum 0
  ss.all.vst.med <- as_tibble(plyr::join_all(ss.all.vst.list.med,
                                            by="external_gene_name",
                                            type="left")) %>%
                                            dplyr::mutate_if(is.numeric, funs(round(., 3)))
  return(list(ss.all.vst, ss.all.vst.med, dds.ss.all))
}
