#! R

##FGSEA analysis on set of genes
##run fgsea on DESeqResults object using particular pathway (typical MSigDB)
#https://stephenturner.github.io/deseq-to-fgsea/
#https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
#also add our own section to find what the sig DE genes are into
fgseamsViperPlot <- function(mrs, mrsstat, geneset=NULL, pathway, qval, TAG, OUTDIR, genecol=NULL) {

  ##mrs is output from msViper, with an es object with specified mrstat as name
  ##mrsstat is statistic used to rank data, comes from msViper
  ##geneset is genes of interest for fgsea (subset to these); use all if not specified
  ##pathway is pathway (from MSigDB in GMT format if *gmt, else named list)
  ##qval is the adjusted p value for all results herein
  ##genecol is the name of the column in tibble with gene names found in pathways

  if(is.null(genecol)){
    genecol <- "external_gene_name"
  }

  if(is.null(geneset)){
    geneset <- names(mrs$es[[mrsstat]])
  }

  if(length(grep(".gmt$", pathway, perl=TRUE))==1){
    pathways.msig <- gmtPathways(pathway)
  }

  ##create ranks
  ranks <- as_tibble(data.frame(stat=mrs$es[[mrsstat]]), rownames="external_gene_name") %>%
           dplyr::filter(external_gene_name %in% geneset) %>%
           dplyr::select(genecol, stat) %>%
           na.omit() %>%
           tibble::deframe()

  ##run, order fgsea
  fgseaRes <- fgsea(pathways=pathways.msig, stats=ranks, nperm=1000000)
  fgseaResLE <- fgseaRes[lapply(fgseaRes$leadingEdge,length)>25,]
  fgseaResTidy <- as_tibble(fgseaResLE) %>%
                  dplyr::mutate(FDR = p.adjust(pval,method="BH")) %>%
                  dplyr::filter(padj < 0.1) %>%
                  dplyr::arrange(desc(NES))

  ##
  pathways.msig.DEsig.absFC2 <- pathways.msig %>% enframe("pathway", genecol) %>%
                    unnest() %>%
                    inner_join(DESeqResults.t, by=genecol) %>%
                    dplyr::mutate(absFC = logratio2foldchange(log2FoldChange)) %>%
                    dplyr::filter(padj < qval) %>%
                    dplyr::filter(abs(absFC) > 2) %>%
                    left_join(., fgseaResTidy, by="pathway") %>%
                    dplyr::select(pathway, genecol, padj.y, NES, padj.x, absFC) %>%
                    dplyr::rename(padj_pathway = "padj.y", padj_gene = "padj.x") %>%
                    dplyr::arrange(NES)

  ##plot
  pathways.msig.DEsig.absFC2.plt <- pathways.msig.DEsig.absFC2 %>%
                                    dplyr::select(pathway, NES, padj_pathway) %>%
                                    dplyr::rename(padj = "padj_pathway") %>%
                                    dplyr::filter(padj < qval) %>%
                                    distinct()

  ggplot(pathways.msig.DEsig.absFC2.plt, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<qval)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(TAG, " pathways NES from GSEA")) +
    theme_minimal() +
    theme(axis.text.y=element_text(size=5))
  ggsave(paste0(OUTDIR, "/fgseaDESeq.", TAG, ".sig.absFC2.plt.pdf"))

  return(pathways.msig.DEsig.absFC2)
}

##run fgsea on each MR set individually, returning list same length as MRs, table of FGSEA
fgseamsViperIndiv <- function(mrs, mrsstat, pathway, qval, TAG, OUTDIR, genecol=NULL) {

  ##mrs is output from msViper, with an es object with specified mrstat as name
  ##mrsstat is statistic used to rank data, comes from msViper
  ##pathway is pathway (from MSigDB in GMT format if *gmt, else named list)
  ##qval is the adjusted p value for all results herein
  ##genecol is the name of the column in tibble with gene names found in pathways

  if(is.null(genecol)){
    genecol <- "external_gene_name"
  }

  if(is.null(geneset)){
    geneset <- names(mrs$es[[mrsstat]])
  }

  if(length(grep(".gmt$", pathway, perl=TRUE))==1){
    pathways.msig <- gmtPathways(pathway)
  }

  ##create ranks
  ranks <- as_tibble(data.frame(stat=mrs$es[[mrsstat]]), rownames="external_gene_name") %>%
           dplyr::filter(external_gene_name %in% geneset) %>%
           dplyr::select(genecol, stat) %>%
           na.omit() %>%
           tibble::deframe()

  ##run, order fgsea
  fgseaRes <- fgsea(pathways=pathways.msig, stats=ranks, nperm=1000000)
  fgseaResLE <- fgseaRes[lapply(fgseaRes$leadingEdge,length)>25,]
  fgseaResTidy <- as_tibble(fgseaResLE) %>%
                  dplyr::mutate(FDR = p.adjust(pval,method="BH")) %>%
                  dplyr::filter(padj < 0.1) %>%
                  dplyr::arrange(desc(NES))

  ##
  pathways.msig.DEsig.absFC2 <- pathways.msig %>% enframe("pathway", genecol) %>%
                    unnest() %>%
                    inner_join(DESeqResults.t, by=genecol) %>%
                    dplyr::mutate(absFC = logratio2foldchange(log2FoldChange)) %>%
                    dplyr::filter(padj < qval) %>%
                    dplyr::filter(abs(absFC) > 2) %>%
                    left_join(., fgseaResTidy, by="pathway") %>%
                    dplyr::select(pathway, genecol, padj.y, NES, padj.x, absFC) %>%
                    dplyr::rename(padj_pathway = "padj.y", padj_gene = "padj.x") %>%
                    dplyr::arrange(NES)

  ##plot
  pathways.msig.DEsig.absFC2.plt <- pathways.msig.DEsig.absFC2 %>%
                                    dplyr::select(pathway, NES, padj_pathway) %>%
                                    dplyr::rename(padj = "padj_pathway") %>%
                                    dplyr::filter(padj < qval) %>%
                                    distinct()

  ggplot(pathways.msig.DEsig.absFC2.plt, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=padj<qval)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(TAG, " pathways NES from GSEA")) +
    theme_minimal() +
    theme(axis.text.y=element_text(size=5))
  ggsave(paste0(OUTDIR, "/fgseaDESeq.", TAG, ".sig.absFC2.plt.pdf"))

  return(pathways.msig.DEsig.absFC2)
}
