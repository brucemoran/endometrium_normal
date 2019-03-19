#
cellMixplot <- function(all.vst, norm40.vst, ss.all.vst.med, combined_clin, BASEDIR, OUTDIR) {

  dir.create(paste0(OUTDIR,"/CellMix"), showWarnings=FALSE, recursive=TRUE)

  ##CellMix
  ##VST pre-removed batch 6
  datasets <- all.vst %>%
              group_by(external_gene_name) %>%
              dplyr::mutate(count=n()) %>%
              dplyr::filter(count==1) %>%
              ungroup() %>%
              dplyr::select(-count) %>%
              na.omit %>%
              as.data.frame()

  datasets <- datasets[table(datasets[,1])==1,]
  rownames(datasets) <- datasets[,1]
  datasets <- datasets[,-1]

  ##signatures to matrix (medians)
  signatures <- base::as.data.frame(ss.all.vst.med[,2:dim(ss.all.vst.med)[2]])
  rownames(signatures) <- unlist(ss.all.vst.med[,1])

  ##CellMix as per authors instructions
  #https://support.bioconductor.org/p/59384/
  x <- as.matrix(datasets)
  y <- as.matrix(signatures)
  res <- gedProportions(x, y, rescale = TRUE, log = FALSE, normalize = TRUE)
  scres <- as_tibble(t(scoef(res)), rownames="sampleID")
  combined_clin_res <- left_join(combined_clin, scres, by = "sampleID") %>%
                       dplyr::select(sampleID, Group, early_epithelium, mid_epithelium, early_stroma, mid_stroma) %>%
                       dplyr::filter(!sampleID %in% as.vector(cond.norm40$sampleID[cond.norm40$Batch==6])) %>%
                       na.omit() %>%
                       dplyr::mutate_if(is.numeric, funs(round(., 3)))

  ##create whichmax, i.e. which tissue/stage is max
  combined_clin_res$whichmax <- suppressWarnings(colnames(combined_clin_res)[apply(combined_clin_res,1,which.max)])
  combined_clin_res <- combined_clin_res %>% dplyr::select(1,2,whichmax,3,4,5,6)

  ##compare to histopath
  stroma.ratio <- read_csv(paste0(BASEDIR,"/data/clinical/Scratch_Endometrium_Samples_Ratio.csv")) %>%
                  dplyr::rename(sampleID = "Sample ID",
                                epithelium_histopath = "Glands_Percentage",
                                stroma_histopath = "Stroma_Percentage")
  combined_clin_res <- left_join(combined_clin_res, stroma.ratio, by = "sampleID")
  combined_clin_res.sum <- combined_clin_res %>%
                           dplyr::mutate(epithelium_cellmix = early_epithelium + mid_epithelium ,
                                         stroma_cellmix = early_stroma + mid_stroma) %>%
                           dplyr::select(1,2, epithelium_cellmix, stroma_cellmix, epithelium_histopath, stroma_histopath)
  combined_clin_res.mlt <- combined_clin_res.sum  %>%
                              dplyr::filter(startsWith(sampleID, "s")) %>%
                              melt() %>%
                              dplyr::mutate(class = 1)

  ##create 'class' i.e. CellMix or histopath fraction
  combined_clin_res.mlt$class <- unlist(lapply(as.vector(combined_clin_res.mlt$variable),function(f){
                                                    strsplit(f, "_")[[1]][2]}))

  plot.order <- combined_clin_res.mlt %>%
                dplyr::filter(variable == "epithelium_histopath") %>%
                dplyr::mutate(rank = dense_rank(value)) %>%
                dplyr::select(rank) %>%
                unlist(c())

  combined_clin_res.mlt$rank <- c(plot.order+50, plot.order+40, plot.order, plot.order+20)
  combined_clin_res.mlt <- combined_clin_res.mlt[order(combined_clin_res.mlt$rank),]

  ##factor to allow intelligible plotting
  combined_clin_res.mlt$sampleID <- factor(combined_clin_res.mlt$sampleID, levels=rev((combined_clin_res.mlt$sampleID)[1:(dim(combined_clin_res.mlt)[1]/4)]), ordered=TRUE)

  ##plot
  ggp <- ggplot(combined_clin_res.mlt, aes(x=class, y=value, fill=variable)) +
         geom_bar(stat="identity") +
         coord_flip() +
         facet_grid(sampleID ~. , scales = "free", space = "free") +
         theme(legend.position = "top")
  ggsave(paste0(OUTDIR,"/CellMix/histopath_vs_CellMix.DC.proportions.pdf"), ggp)

  #epithelial fractions (1-this = stromal)
  epi_histop <- combined_clin_res.mlt %>%
                dplyr::filter(variable == "epithelium_histopath") %>%
                dplyr::select(value) %>% unlist(c())
  epi_cellmx <- combined_clin_res.mlt %>%
                dplyr::filter(variable == "epithelium_cellmix") %>%
                dplyr::select(value) %>% unlist(c())

  ##correlation, T-test
  ctres <- cor.test(epi_histop, epi_cellmx)
  ttres <- t.test(epi_histop, epi_cellmx)

  ##data frame to plot
  dfp <- data.frame(histop=epi_histop, cellmx=epi_cellmx)
  dfpo <- dfp[order(dfp$histop),]
  dfpo$y <- seq(from=1, to=length(epi_histop), by=1)
  dfpl <- data.frame(epithelium_fraction=c(dfpo$histop,dfpo$cellmx),
                     x=c(dfpo$y, dfpo$y),
                     class=c(rep("histopath", times=length(epi_histop)),
                             rep("cellmix", times=length(epi_histop))))
  coefp <- coef(lm(dfpo$histop ~ dfpo$cellmx))
  summlm <- summary(lm(dfpo$histop ~ dfpo$cellmx))

  ##plot
  ggpd <- ggplot(dfpl, aes(x=x, y=epithelium_fraction, colour=class)) +
          geom_point() +
          geom_abline(intercept=coefp[1], slope=summlm$r.squared)
  ggsave(paste0(OUTDIR,"/CellMix/histopath_vs_CellMix.DC.epithelium_fraction_rsquared.pdf"), ggpd)

  ##heatmap of Stage vs. classification
  heatmap_clin <- combined_clin_res %>% dplyr::filter(! 1:length(Group) %in% starts_with("Stage", vars=Group)) %>%
                                            dplyr::select(1,4,5,6,7) %>%
                                            as.data.frame()
  rownames(heatmap_clin) <- heatmap_clin[,1]
  heatmap_clin <- heatmap_clin[,-1]
  aheatmap.clin <- aheatmap(t(heatmap_clin), fontsize=8, cexRow=1)

  pdf(paste0(OUTDIR,"/CellMix/CellMix.Norm40.heatmap.pdf"), onefile=F)
    aheatmap.clin
  dev.off()

  ##########
  ## TCGA ##
  ##########
  ##define clinical sets
  ucs_tcga_clin <- combined_clin_res %>%
               dplyr::filter(sampleID %in% NIHMS891671_clin$bcr_patient_barcode) %>%
               left_join(., NIHMS891671_clin, by=c("sampleID" = "bcr_patient_barcode")) %>%
               dplyr::mutate(Stage = paste0("Stage_",unlist(lapply(gsub("[ABC12]","",gsub("Stage ","", clinical_stage)),function(f){as.numeric(as.roman(f))})))) %>%
               dplyr::mutate(MMMT2type = unlist(lapply(MMMTtype, function(f){if(f=="NOS"){return("NOS")}else{return("HET/HOM")}}))) %>%
               dplyr::select(sampleID, MMMT2type, Stage, whichmax, early_epithelium, mid_epithelium, early_stroma, mid_stroma) %>%
               dplyr::arrange(MMMT2type)

  ucec_tcga_clin <- combined_clin_res %>%
               dplyr::filter(sampleID %in% NATURE12113_clin$bcr_patient_barcode) %>%
               left_join(., NATURE12113_clin, by=c("sampleID" = "bcr_patient_barcode")) %>%
               dplyr::select(sampleID, histology, `2009stagegroup`, IntegrativeCluster, whichmax, early_epithelium, mid_epithelium, early_stroma, mid_stroma) %>%
               dplyr::rename(Stage = `2009stagegroup`, Histology = histology) %>%
               dplyr::arrange(Histology)

  ##heatmap to show early, mid stroma and MMMTtype in UCS
  annCols <- list(ucs_tcga_clin$Stage,ucs_tcga_clin$MMMT2type)
  names(annCols) <- c("Stage", "MMMT type")
  aheatmap.ucs <- aheatmap(t(ucs_tcga_clin[,c("early_epithelium", "mid_epithelium", "early_stroma", "mid_stroma")]),
           annCol=annCols,
           annColors=list(c("white", "grey88", "grey55", "black"), c("forestgreen", "purple")))

  pdf(paste0(OUTDIR,"/CellMix/Stage_MMMT.TCGA-UCS.4cellphase.heatmap.pdf"), onefile=F)
    aheatmap.ucs
  dev.off()

  ##heatmap to show early, mid stroma and histology in UCEC
  annCols <- list(ucec_tcga_clin$Stage, ucec_tcga_clin$Histology, ucec_tcga_clin$IntegrativeCluster)
  names(annCols) <- c("Stage", "Histology", "IntegrativeCluster")
  aheatmap.ucec <- aheatmap(t(ucec_tcga_clin[,c("early_epithelium", "mid_epithelium", "early_stroma", "mid_stroma")]),
           annCol=annCols,
           annColors=list(c("white", "grey88", "grey55", "black"), c("forestgreen", "purple", "orange")))

  pdf(paste0(OUTDIR,"/CellMix/Stage_Histology.TCGA-UCEC.4cellphase.heatmap.pdf"), onefile=F)
    aheatmap.ucec
  dev.off()

  ##FET over whichmax vs. MMMT2type, Stage
  ##none significant
  fetUCS.list <- fetUCEC.list <- as.list(1:3)
  fet.ucs.mmmt.wm <- table(ucs_tcga_clin$MMMT2type, ucs_tcga_clin$whichmax)
  fet.ucs.stage.wm <- table(ucs_tcga_clin$Stage, ucs_tcga_clin$whichmax)
  fet.ucs.stage.mmmt <- table(ucs_tcga_clin$Stage, ucs_tcga_clin$MMMT2type)
  fetUCS.list[[1]] <- fisher.test(fet.ucs.mmmt.wm)
  fetUCS.list[[2]] <- fisher.test(fet.ucs.stage.mmmt)
  fetUCS.list[[3]] <- fisher.test(fet.ucs.stage.wm)

  ##Stage vs. histology is significant
  fet.ucec.stage.wm <- table(ucec_tcga_clin$Stage, ucec_tcga_clin$whichmax)
  fet.ucec.stage.hist <- table(ucec_tcga_clin$Stage, ucec_tcga_clin$Histology)
  fet.ucec.hist.wm <- table(ucec_tcga_clin$Histology, ucec_tcga_clin$whichmax)
  fetUCEC.list[[1]] <- chisq.test(fet.ucec.stage.wm)
  fetUCEC.list[[2]] <- fisher.test(fet.ucec.stage.hist)
  fetUCEC.list[[3]] <- fisher.test(fet.ucec.hist.wm)

  ##return
  cellMixList <- list(datasets, signatures, combined_clin_res,
                      list(ctres, ttres),
                      list(ggp, ggpd, aheatmap.clin, aheatmap.ucs, aheatmap.ucec),
                      fetUCS.list,
                      fetUCEC.list)
  return(cellMixList)

}
