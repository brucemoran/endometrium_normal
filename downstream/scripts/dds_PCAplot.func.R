ddsPCA <- function(countsIn, clinIn, OUTDIR, topn=NULL, savePlots=NULL) {
  print("Making PCA plots...")

  if(is.null(topn)){
    topn <- 500
  }
  if(is.null(savePlots)){
    savePlots <- TRUE
  }
  dir.create(paste0(OUTDIR, "/PCA"), showWarnings=FALSE, recursive=TRUE)

  ##matrix input for DESeq2
  nz.co.all <- base::as.matrix(countsIn[,unlist(lapply(1:dim(countsIn)[2],
                            function(f){is.numeric(unlist(countsIn[,f]))}))],
                            dimnames=list(unlist(countsIn$ensembl_gene_id),
                            colnames(countsIn)[3:dim(countsIn)[2]]))
  rownames(nz.co.all) <- unlist(countsIn$ensembl_gene_id)

  ##clinical
  cond.all <- clinIn %>% dplyr::filter(sampleID %in% colnames(nz.co.all)) %>%
                         dplyr::mutate(Group = gsub(" ", "_", Group)) %>%
                         dplyr::mutate(Type = unlist(lapply(Batch,function(f){if(f>7){return("Tumour")}else{return("Normal")}}))) %>%
                         as.data.frame()
  cond.norm <- factorise(cond.all[! cond.all$Batch %in% c(8,9),c(1,2,4,7)])
  cond.norm67 <- factorise(cond.norm[cond.norm$Batch %in% c(6,7),])
  cond.tcga <-  factorise(cond.all[cond.all$Batch %in% c(8,9),c(1,2,4,7)])
  cond.norm.rm6 <- factorise(cond.norm[cond.norm$Batch %in% c(1,2,3,7),])

  ##count subsets
  print("Subsetting count data...")
  nz.co.norm <- nz.co.all[,colnames(nz.co.all) %in% cond.norm$sampleID]
  nz.co.norm67 <- nz.co.all[,colnames(nz.co.all) %in% cond.norm67$sampleID]
  nz.co.tcga <- nz.co.all[,colnames(nz.co.all) %in% cond.tcga$sampleID]
  nz.co.norm.rm6 <- nz.co.all[,colnames(nz.co.all) %in% cond.norm.rm6$sampleID]

  ##DDS objects
  print("Creating dds's...")
  dds.all <- DESeqDataSetFromMatrix(countData=nz.co.all,
                                    colData=cond.all,
                                    design=~Type)
  dds.norm <- DESeqDataSetFromMatrix(countData=nz.co.norm,
                                  colData=cond.norm,
                                  design=~Type)
  dds.norm67 <- DESeqDataSetFromMatrix(countData=nz.co.norm67,
                                      colData=cond.norm67,
                                      design=~Type)
  dds.tcga <- DESeqDataSetFromMatrix(countData=nz.co.tcga,
                                     colData=cond.tcga,
                                     design=~Type)
  dds.norm.rm6 <- DESeqDataSetFromMatrix(countData=nz.co.norm.rm6,
                                         colData=cond.norm.rm6,
                                         design=~Type)
  print("estimateSizeFactors...")
  dds.all <- estimateSizeFactors(dds.all)
  dds.norm <- estimateSizeFactors(dds.norm)
  dds.norm67 <- estimateSizeFactors(dds.norm67)
  dds.tcga <- estimateSizeFactors(dds.tcga)
  dds.norm.rm6 <- estimateSizeFactors(dds.norm.rm6)

  ##we find 4 is separated from 1, 2, 3, 7; remove 6
  print("Norm without Batch6 Samples...")
  nz.co.norm.rm6.vst <- vstPCAPlot(dds=dds.norm.rm6,
                                   intgroup=c("Batch", "Group"),
                                   PLOTDIR=paste0(OUTDIR, "/PCA"),
                                   TAG="Group-Batch.Norm.rm6",
                                   ntop=topn)
  nz.co.norm.rm6.vst[[3]] <- dds.norm.rm6
  nz.co.norm.rm6.vst[[4]] <- cond.norm.rm6

  ##run on all samples together (above first to tst working)
  print("All Samples...")
  nz.co.all.vst <- vstPCAPlot(dds=dds.all,
                              intgroup=c("Batch", "Group"),
                              PLOTDIR=paste0(OUTDIR, "/PCA"),
                              TAG="Group-Batch.All",
                              ntop=topn)
  nz.co.all.vst[[3]] <- dds.all
  nz.co.all.vst[[4]] <- cond.all

  ##remove tumour to get a look at normals
  ##we find 6 is separated from 1, 2, 3, 7; remove
  print("Normal Samples...")
  nz.co.norm.vst <- vstPCAPlot(dds=dds.norm,
                               intgroup=c("Batch", "Group"),
                               PLOTDIR=paste0(OUTDIR, "/PCA"),
                               TAG="Group-Batch.Norm",
                               ntop=topn)
  nz.co.norm.vst[[3]] <- dds.norm
  nz.co.norm.vst[[4]] <- cond.norm

  ##split into only 6,7 highlights batch of public data
  print("GSE98386 Samples...")
  nz.co.norm67.vst <- vstPCAPlot(dds=dds.norm67,
                                 intgroup=c("Batch", "Group"),
                                 PLOTDIR=paste0(OUTDIR, "/PCA"),
                                 TAG="Group-Batch.norm67",
                                 ntop=topn)
  nz.co.norm67.vst[[3]] <- dds.norm67
  nz.co.norm67.vst[[4]] <- cond.norm67

  ##just tumour data
  print("TCGA UCS and UCEC Samples...")
  nz.co.tcga.vst <- vstPCAPlot(dds=dds.tcga,
                               intgroup=c("Batch", "Group"),
                               PLOTDIR=paste0(OUTDIR, "/PCA"),
                               TAG="Group-Batch.tcga",
                               ntop=topn)
  nz.co.tcga.vst[[3]] <- dds.tcga
  nz.co.tcga.vst[[4]] <- cond.tcga

  ##diagnostic plot of VST
  ##counts and conds
  print("Filtering normals...")
  nz.co.norm.rm6.vst.ass <- assay(nz.co.norm.rm6.vst[[2]])
  nz.co.norm.rm6.vst.mlt <- melt(nz.co.norm.rm6.vst.ass) %>%
                            dplyr::rename(vst_expr = value, sampleID = Var2)
  rm.norm.vst <- mean(colSums(nz.co.norm.rm6.vst.ass))-(sd(colSums(nz.co.norm.rm6.vst.ass))*2)
  nz.co.norm40.rm6.vst.ass <- nz.co.norm.rm6.vst.ass[, grep("TRUE", colSums(nz.co.norm.rm6.vst.ass) > rm.norm.vst)]

  ##norm40.rm6
  rm.norm40.vst.ass <- mean(colSums(nz.co.norm40.rm6.vst.ass))-(sd(colSums(nz.co.norm40.rm6.vst.ass))*2)
  nz.co.norm40.rm6.vst.ass <- nz.co.norm40.rm6.vst.ass[, grep("TRUE", colSums(nz.co.norm40.rm6.vst.ass) > rm.norm40.vst.ass)]
  nz.co.norm40.rm6.vst.mlt <- melt(nz.co.norm40.rm6.vst.ass) %>%
                            dplyr::rename(vst_expr = value, sampleID = Var2 )

  ##we find 6 is separated from 1, 2, 3, 7; remove 6
  cond.norm40.rm6 <- factorise(cond.norm[cond.norm$Batch %in% c(1,2,3,7) &
                                         cond.norm$sampleID %in% levels(nz.co.norm40.rm6.vst.mlt$sampleID),])
  nz.co.norm40.rm6 <- nz.co.all[,colnames(nz.co.all) %in% cond.norm40.rm6$sampleID]
  dds.norm40.rm6 <- DESeqDataSetFromMatrix(countData=nz.co.norm40.rm6,
                                  colData=cond.norm40.rm6,
                                  design=~1)
  dds.norm40.rm6 <- estimateSizeFactors(dds.norm40.rm6)

  ##final set used for ongoing analysis
  print("Norm_40 Samples...")
  nz.co.norm40.rm6.vst <- vstPCAPlot(dds=dds.norm40.rm6,
                                     intgroup=c("Batch", "Group"),
                                     PLOTDIR=paste0(OUTDIR, "/PCA"),
                                     TAG="Group-Batch.norm40.rm7",
                                     ntop=topn)
  nz.co.norm40.rm6.vst[[3]] <- dds.norm40.rm6
  nz.co.norm40.rm6.vst[[4]] <- cond.norm40.rm6

  nz.co.norm40.rm6.vst.ass <- assay(nz.co.norm40.rm6.vst[[2]])
  nz.co.norm40.rm6.vst.mlt <- melt(nz.co.norm40.rm6.vst.ass) %>%
                              dplyr::rename(vst_expr = value, sampleID = Var2 )

  print("VST distribution...")
  densplot <- ggplot(nz.co.norm.rm6.vst.mlt, aes(vst_expr, colour=sampleID)) +
  geom_density() +
  labs(title="VST expression per sample", subtitle="Irregular curve indicates filtered samples SRR5488689, SRR548868990")
  ggsave(paste0(OUTDIR, "/supplementary.VST-density.norm.rm6.pdf"), densplot)

  ##remove samples falling twice SD below mean of total VST per sample
  ##2 removed, look like errant
  print("VST distribution Norm_40...")
  densplot_norm40 <- ggplot(nz.co.norm40.rm6.vst.mlt, aes(vst_expr, colour=sampleID)) +
  geom_density() +
  labs(title="VST expression per sample", subtitle="Post filtering using VST mean - 2*SD VST")
  ggsave(paste0(OUTDIR,"/supplementary.VST-density.norm40.rm6.pdf"), densplot_norm40)

  ##outputs
  print("Return output and saving...")
  ddsPCAList <- list(nz.co.all.vst, nz.co.tcga.vst, nz.co.norm.vst, nz.co.norm67.vst, nz.co.norm.rm6.vst, nz.co.norm40.rm6.vst, densplot, densplot_norm40)

  return(ddsPCAList)
}
