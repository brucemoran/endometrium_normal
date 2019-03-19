
cellphaseLimmaCAMERA <- function(combined_clin_res, norm40.counts, gene2name, OUTDIR, PATHWAYDIR){

  DEOUTDIR <- paste0(OUTDIR, "/Limma_CAMERA")
  dir.create(DEOUTDIR, showWarnings=FALSE, recursive=TRUE)

    #################################
   ## Limma epithelium vs. stroma ##
  #################################

  #run DESeq2 on epithelium whichmax vs. stroma whichmax removing tumours
  norm.epi_str <- norm40.counts
  combined_clin_res.norm40.epi_str <- combined_clin_res %>%
                                      dplyr::filter(sampleID %in% colnames(norm.epi_str)) %>%
                                      dplyr::mutate(ConsCellType = unlist(lapply(whichmax, function(f){
                                        strsplit(f,"_")[[1]][2]
                                      }))) %>%
                                      dplyr::select(sampleID, ConsCellType, Group, whichmax) %>%
                                      base::as.data.frame()

  ##interaction of diet and tissue, so can look within Diet groups at Tissue effect on DGE
  design <- model.matrix(~0+whichmax, combined_clin_res.norm40.epi_str)
  dge <- DGEList(counts=norm.epi_str)
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)

#   ##MDS plots
  col.whichmax <- factor(combined_clin_res.norm40.epi_str$whichmax)
  levels(col.whichmax) <-  sample(rainbow(nlevels(col.whichmax)))
  col.whichmax <- as.character(col.whichmax)
  lcpm <- cpm(dge, log=TRUE, prior.count=0.5, group="whichmax")
  whichmaxN <- as.vector(combined_clin_res.norm40.epi_str$whichmax)
#   pdf(paste0(DEOUTDIR,"/MDS.Tissue.limma-voom.Individual.pdf"))
#     plotMDS(lcpm, labels=whichmaxN, col=col.whichmax)
#     title(main="A. whichmax groups")
#   dev.off()

#   ##voom
  dcv <- voom(dge, design, plot=FALSE)
  dcfit <- lmFit(dcv, design)
  dcmc1 <- makeContrasts(
            EarlyS_vs_MidS = whichmaxearly_stroma - whichmaxmid_stroma,
            levels=colnames(design))
  dcmc2 <- makeContrasts(
            EarlyS_vs_MidE = whichmaxearly_stroma - whichmaxmid_epithelium,
            levels=colnames(design))
  dcfitmc1 <- contrasts.fit(dcfit, dcmc1)
  dcfitmc2 <- contrasts.fit(dcfit, dcmc2)
  dcfitmc1 <- eBayes(dcfitmc1, robust=TRUE)
  dcfitmc2 <- eBayes(dcfitmc2, robust=TRUE)

#   summary(decideTests(dcfitmc1, adjust.method="fdr", p.value=0.01, lfc=2))
#   summary(decideTests(dcfitmc2, adjust.method="fdr", p.value=0.01, lfc=2))

  EarlyS_vs_MidS <- as_tibble(topTreat(dcfitmc1, coef=1, n=Inf), rownames="ensembl_gene_id") %>%
                    dplyr::mutate(absFC = 2^logFC) %>%
                    dplyr::arrange(t) %>%
                    dplyr::mutate(fgseaRank = seq(1:length(t))) %>%
                    left_join(., gene2name, by = "ensembl_gene_id") %>%
                    na.omit() %>%
                    dplyr::select(ensembl_gene_id, external_gene_name, logFC, absFC, adj.P.Val, fgseaRank, t)

  EarlyS_vs_MidE <- as_tibble(topTreat(dcfitmc2, coef=1, n=Inf), rownames="ensembl_gene_id") %>%
                    dplyr::mutate(absFC = 2^logFC) %>%
                    dplyr::arrange(t) %>%
                    dplyr::mutate(fgseaRank = seq(1:length(t))) %>%
                    left_join(., gene2name, by = "ensembl_gene_id") %>%
                    na.omit() %>%
                    dplyr::select(ensembl_gene_id, external_gene_name, logFC, absFC, adj.P.Val, fgseaRank, t)

  ##CAMERA
  cameraListo <- runCAMERA(RES=EarlyS_vs_MidS,
                          GENECOL="ensembl_gene_id",
                          DGE=dge,
                          DESIGN=design,
                          CONTRAST=dcmc1,
                          MSIGDIR="/data/genome/reference/pathways/msigdb",
                          GENOME=GENOME,
                          SIGP=0.05,
                          MSIGMATCH=".")

  CAMERAOUTDIR <- paste0(OUTDIR, "/CAMERA-pathways")
  dir.create(CAMERAOUTDIR, showWarnings=FALSE)
  for(x in 1:length(cameraListo)){
    if(dim(cameraListo[[x]])[1]!=0){
      write_tsv(as.data.frame(cameraListo[[x]]), path=paste0(gsub(".gmt",".tsv",paste0(CAMERAOUTDIR,"/nz.co.norm40.", names(cameraListo)[x]))))
    }
  }
  return(list(EarlyS_vs_MidS, EarlyS_vs_MidE, cameraListo))
}
