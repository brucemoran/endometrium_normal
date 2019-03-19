ddsPCAOne <- function(countsIn, clinIn, OUTDIR, topn=NULL, savePlots=NULL) {
  print("Making PCA plots...")

  if(is.null(topn)){
    topn <- 500
  }
  if(is.null(savePlots)){
    savePlots <- TRUE
  }
  dir.create(paste0(OUTDIR, "/PCA"), showWarnings=FALSE, recursive=TRUE)

  ##matrix input for DESeq2
  counts.all <- base::as.matrix(countsIn[,unlist(lapply(1:dim(countsIn)[2],
                            function(f){is.numeric(unlist(countsIn[,f]))}))],
                            dimnames=list(unlist(countsIn$ensembl_gene_id),
                            colnames(countsIn)[3:dim(countsIn)[2]]))
  rownames(counts.all) <- unlist(countsIn$ensembl_gene_id)

  ##clinical
  cond.all <- clinIn %>% dplyr::filter(sampleID %in% colnames(counts.all)) %>%
                         dplyr::mutate(Group = gsub(" ", "_", Group)) %>%
                         dplyr::mutate(Type = unlist(lapply(Batch,function(f){if(f>7){return("Tumour")}else{return("Normal")}})))

  cond.all$Type <- factor(cond.all$Type)
  cond.all$Batch <- factor(cond.all$Batch)

  samp.norm <- cond.all %>% dplyr::filter(! Batch %in% c(8,9)) %>% dplyr::select(1) %>% unlist()
  samp.tcga <- cond.all %>% dplyr::filter(Batch %in% c(8,9)) %>% dplyr::select(1) %>% unlist()
  samp.norm67 <- cond.all %>% dplyr::filter(Batch %in% c(6,7)) %>% dplyr::select(1) %>% unlist()
  samp.norm.rm6 <- cond.all %>% dplyr::filter(! Batch %in% c(6,8,9)) %>% dplyr::select(1) %>% unlist()

  ##DDS object
  print("Creating dds's...")
  dds.all <- DESeqDataSetFromMatrix(countData=counts.all,
                                    colData=cond.all,
                                    design=~Type)
  print("estimateSizeFactors...")
  dds.all <- estimateSizeFactors(dds.all)
  print("estimateDispersions...")
  dds.loc <- estimateDispersions(dds.all, fitType="local")
  dds.par <- estimateDispersions(dds.all, fitType="parametric")

  ##test parametric vs. local fit of dispersion, choosing lower
  med.loc <- median(abs(log(rowData(loc.dds.all)$dispGeneEst) - log(rowData(loc.dds.all)$dispFit)))
  med.par <- median(abs(log(rowData(par.dds.all)$dispGeneEst) - log(rowData(par.dds.all)$dispFit)))
  dds.chose <- dds.par
  if(med.par > med.loc){dds.chose <- dds.loc}

  print("varianceStabilizingTransformation...")
  vst.dds.all <- varianceStabilizingTransformation(dds.chose, blind=FALSE)

  ##placeholders for output
  all.list <- tcga.list <- norm.list <- norm67.list <- norm.rm6.list <- norm40.list <- as.list(1:4)

  ##run on all samples together (above first to tst working)
  print("All Samples...")
  all.list[[1]] <- dds.chose
  all.list[[2]] <- cond.all
  all.list[[3]] <- vst.dds.all
  all.list[[4]] <- vstPCAPlot(vst=all.list[[3]],
                              intgroup=c("Batch", "Group"),
                              PLOTDIR=paste0(OUTDIR, "/PCA"),
                              TAG="Group-Batch.All",
                              ntop=topn)

  ##just tumour data
  print("TCGA UCS and UCEC Samples...")
  tcga.list[[1]] <- dds.all[,samp.tcga]
  tcga.list[[2]] <- cond.all %>% dplyr::filter(sampleID %in% samp.tcga)
  tcga.list[[3]] <- vst.dds.all[,samp.tcga]
  tcga.list[[4]] <- vstPCAPlot(vst=tcga.list[[3]],
                               intgroup=c("Batch", "Group"),
                               PLOTDIR=paste0(OUTDIR, "/PCA"),
                               TAG="Group-Batch.tcga",
                               ntop=topn)

  ##remove tumour to get a look at normals
  ##we find 6 is separated from 1, 2, 3, 7; remove
  print("Normal Samples...")
  norm.list[[1]] <- dds.all[,samp.norm]
  norm.list[[2]] <- cond.all %>% dplyr::filter(sampleID %in% samp.norm)
  norm.list[[3]] <- vst.dds.all[,samp.norm]
  norm.list[[4]] <- vstPCAPlot(vst=norm.list[[3]],
                         intgroup=c("Batch", "Group"),
                         PLOTDIR=paste0(OUTDIR, "/PCA"),
                         TAG="Group-Batch.Norm",
                         ntop=topn)

  ##split into only 6,7 highlights batch of public data
  print("GSE98386 Samples...")
  norm67.list[[1]] <- dds.all[,samp.norm67]
  norm67.list[[2]] <- cond.all %>% dplyr::filter(sampleID %in% samp.norm67)
  norm67.list[[3]] <- vst.dds.all[,samp.norm67]
  norm67.list[[4]] <- vstPCAPlot(vst=norm.list[[3]],
                               intgroup=c("Batch", "Group"),
                               PLOTDIR=paste0(OUTDIR, "/PCA"),
                               TAG="Group-Batch.norm67",
                               ntop=topn)

  ##we find 4 is separated from 1, 2, 3, 7; remove 6
  print("Norm without Batch6 Samples...")
  norm.rm6.list[[1]] <- dds.all[,samp.norm.rm6]
  norm.rm6.list[[2]] <- cond.all %>% dplyr::filter(sampleID %in% samp.norm.rm6)
  norm.rm6.list[[3]] <- vst.dds.all[,samp.norm.rm6]
  norm.rm6.list[[4]] <- vstPCAPlot(vst=norm.rm6.list[[3]],
                                   intgroup=c("Batch", "Group"),
                                   PLOTDIR=paste0(OUTDIR, "/PCA"),
                                   TAG="Group-Batch.Norm.rm6",
                                   ntop=topn)

  ##diagnostic plot of VST
  ##counts and conds
  print("Filtering normals...")
  norm.rm6.vst.ass <- assay(norm.rm6.list[[3]])
  norm.rm6.vst.mlt <- melt(norm.rm6.vst.ass) %>%
                      dplyr::rename(vst_expr = value, sampleID = Var2)
  rm6.norm.vst <- mean(colSums(norm.rm6.vst.ass))-(sd(colSums(norm.rm6.vst.ass))*2)
  norm40.vst.ass <- norm.rm6.vst.ass[, grep("TRUE", colSums(norm.rm6.vst.ass) > rm6.norm.vst)]
  samp.norm40 <- colnames(norm40.rm6.vst.ass)
  norm40.vst.mlt <- norm.rm6.vst.mlt %>% dplyr::filter(sampleID %in% samp.norm40)


  ##final set used for ongoing analysis
  print("Norm40 Samples...")
  norm40.list[[1]] <- dds.all[,samp.norm40]
  norm40.list[[2]] <- cond.all %>% dplyr::filter(sampleID %in% samp.norm40)
  norm40.list[[3]] <- vst.dds.all[,samp.norm40]
  norm40.list[[4]] <- vstPCAPlot(vst=norm40.list[[3]],
                                 intgroup=c("Batch", "Group"),
                                 PLOTDIR=paste0(OUTDIR, "/PCA"),
                                 TAG="Group-Batch.norm40.rm6",
                                 ntop=topn)

  print("VST distributions...")
  densplot.norm <- ggplot(norm.rm6.vst.mlt, aes(vst_expr, colour=sampleID)) +
  geom_density() +
  labs(title="VST expression per sample", subtitle="Irregular curve indicates filtered samples SRR5488689, SRR548868990")
  ggsave(paste0(OUTDIR, "/supplementary.VST-density.norm.rm6.pdf"), densplot.norm)

  ##remove samples falling twice SD below mean of total VST per sample
  ##2 removed, look like errant
  densplot.norm40 <- ggplot(norm40.vst.mlt, aes(vst_expr, colour=sampleID)) +
  geom_density() +
  labs(title="VST expression per sample", subtitle="Post filtering using VST mean - 2*SD VST")
  ggsave(paste0(OUTDIR,"/supplementary.VST-density.norm40.pdf"), densplot.norm40)

  ##outputs
  print("Return output and saving...")
  ddsPCAList <- list(all.list,
                     tcga.list,
                     norm.list,
                     norm67.list,
                     norm.rm6.list,
                     norm40.list,
                     densplot.norm,
                     densplot.norm40)

  return(ddsPCAList)
}
