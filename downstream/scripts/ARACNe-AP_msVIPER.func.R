#
runAracneAP <- function(datasets, signatures, TAG, BASEDIR, OUTDIR){

  ##TAG creates dir name in ARACNE-AP/ in OUTDIR
  dir.create(paste0(OUTDIR, "/ARACNe-AP/", TAG), showWarnings=FALSE, recursive=TRUE)

  ###############
  ## ARACNe-AP ##
  ###############
  ##runs shell script in BASEDIR/scripts/
  #output VST, signature genelist for use by ARACNe-AP
  datasetscol <- data.frame(gene=rownames(datasets))
  datasetscol <- cbind(datasetscol, datasets)
  write.table(datasetscol,
              file=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/dataset.vst.tsv"),
              quote=F, sep="\t", row=F)
  write.table(rownames(signatures),
              file=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/signature.genes.txt"),
              quote=F, sep="\t", row=F, col=F)


  ##test for presence of a network
  ##ARACNe-AP still requires some heavy lifting, so this stymies accidental rerun
  if(length(grep("network_ready.txt", dir(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/bootstraps"))))==0){
     print("Network file for ARACNe-AP not found, creating...")
     system(paste0("sh ", BASEDIR, "/scripts/ARACNe_AP_run.sh ", OUTDIR, " ", TAG))
  }

  ##the following replaces last '}' and was used interactively to allow running of msViper\
  ##following generation of the nextwork.txt file
  ##with multiple runs of ARACNe required, we allow the jobs to run, then use an \
  ##if/else for presence of the network.txt file to run msViper

  #   while (!file.exists(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/bootstraps/network.txt"))) {
  #     Sys.sleep(300)
  #   }
  # }
}

runMsViper <- function(combined_clin_res, TAG, BASEDIR, OUTDIR, DERES, CGCTABLE, minconx=NULL){

  if(is.null(minconx)){
    minconx <- 0
  }
  ###########
  ## VIPER ##
  ###########
  ## create regulon, waiting for it to be made
  while(!file.exists(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/bootstraps/network_ready.txt"))){
    Sys.sleep(300)
  }
  ##NB that 'network.txt' file is moved to network_ready.txt
  ##possible that ARACNe-AP 'consolidate' produces intermediate network.txt,
  ##so mv to new file to specify consolidate fully complete

  regulonaracne <- aracne2regulon(afile=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/bootstraps/network_ready.txt"),
                                  eset=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/dataset.vst.tsv"))

  ## create eset
  vsts <- read.table(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/dataset.vst.tsv"), row=1, header=T)
  pheno <- combined_clin_res %>%
           dplyr::mutate(sampleID = gsub("-",".",sampleID)) %>%
           dplyr::filter(sampleID %in% colnames(vsts)) %>%
           as.data.frame()
  pheno$description <- c(rep("normal", time=length(grep("Stage", pheno$Group, invert=T))),
                         rep("tumour", time=length(grep("Stage", pheno$Group))))
  pData <- new("AnnotatedDataFrame",
                data=data.frame(row.names = colnames(vsts),
                      description = as.vector(pheno$description)))
  combined.eset <- ExpressionSet(assayData=as.matrix(vsts),
                             phenoData=pData)

  ## generate signature for normal and tumour
  signature <- viper::rowTtest(combined.eset, "description", "normal", "tumour")
  rownamesp <- rownames(signature$p.value)
  signature$p.value <- as.matrix(p.adjust(signature$p.value, method="BH"))
  rownames(signature$p.value) <- rownamesp
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  nullmodel <- viper::ttestNull(combined.eset, "description", "normal", "tumour", per=5000, repos=TRUE, verbose=TRUE)

  ## msViper, leading edge analyses + bootstrapped
  mrs <- msviper(signature, regulonaracne, nullmodel, minsize=10, verbose=TRUE)
  mrs <- viper::ledge(mrs)
  mrs <- viper::shadow(mrs, regulators=0.01, shadow=0.01, verbose=FALSE, per=5000)
  mrsc <- msviperCombinatorial(mrs, regulators=0.01, verbose=FALSE)

  ##if statement to deny error from synergy analysis which otherwise exits from the script
  if(length(mrsc$regulon)>length(mrs$regulon)){
    mrs <- msviperSynergy(mrsc, verbose = FALSE)
    ##BH adjust p-values
    mrs$es$q.value <- p.adjust(mrs$es$p.value,method="BH")
  }
  if(length(mrsc$regulon)==length(mrs$regulon)){
    ##BH adjust p-values
    mrs$es$q.value <- p.adjust(mrs$es$p.value, method="BH")
  }

  ## outputs
  pdf(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/", TAG, ".msViper.syn.results.pdf"))
    plot(mrs, pval=mrs$es$q.value, cex = .7)
  dev.off()

  mrstab <- data.frame(NES=mrs$es$nes,
                       Size=mrs$es$size,
                       p.value=mrs$es$p.value,
                       q.value=mrs$es$q.value)

  write.table(mrstab, paste0(OUTDIR, "/ARACNe-AP/", TAG, "/", TAG, ".msViper.results.tsv"),
              quote=F, row=T, col=T, sep="\t")

  ##write output for significant 'mrs' in terms of leading edge (targets)
  ##levels of significance
  qvals <- c(0.01)
  for (pp in 1:length(qvals)){

    ##subset of signifcant MRs
    mrsnms <- sort(mrs$es$q.value) < qvals[pp]

    ##leading edge (ledge) output
    ##NB ledge is used to determine the set of targets in the signature of MRs
    ledgenms <- sort(grep("--", names(mrsnms[mrsnms=="TRUE"]), invert=T, value=T))
    mrsldgList <- lapply(ledgenms, function(f){
      if(! is.na(names(mrs$ledge[f]))){
        # if(names(mrs$ledge[f]) %in% as.vector(unlist(CGC["Gene Symbol"])) | length(unique(c(unlist(mrs$ledge[f])) %in% as.vector(unlist(CGC["Gene Symbol"])))) == 2 | length( unique(c(unlist(mrs$ledge[f])) %in% as.vector(unlist(CGC["Gene Symbol"]))) == "TRUE") == 1) {
          return(as.vector(unlist(mrs$ledge[f])))
        # }
      }
    })
    names(mrsldgList) <- ledgenms
    ledgeall <- sort(c(unique(unlist(mrsldgList)), ledgenms))

    ##FGSEA all and tabulate
    
    ##remove elements with no entries
    mrsldgList <- lapply(mrsldgList,function(f){if(identical(f, character(0))=="FALSE"){return(f)}})
    mrsldgList[sapply(mrsldgList, is.null)] <- NULL

    ##write those ledges out
    # for(x in 1:length(mrsldgList)){
    #   dir.create(paste0(OUTDIR, "/ARACNe-AP/", TAG, "/analysis/msViper/leading_edge/", TAG, "/", qvals[pp]), showWarnings=F, recursive=TRUE)
    #   write.table(mrsldgList[[x]],
    #               paste0(OUTDIR, "/ARACNe-AP/", TAG, "/analysis/msViper/leading_edge/", TAG, "/", qvals[pp], "/", names(mrsldgList)[x], ".ledge.msViper.genes.p_",qvals[pp],".txt"),
    #               quote=F, row=F, col=F, sep="\t")
    # }

    ## read DE results, parse ledge data for nodes
    de.results <- read_tsv(DERES)
    colnames(de.results)[1] <- "external_gene_name"
    de.padj.results <- de.results %>% dplyr::filter(padj < 0.01)
    nodenms <-  as.vector(unlist(de.padj.results$external_gene_name[de.padj.results$external_gene_name %in% names(cgcmrsldgList)]))
    ##or not requiring DE
    nodenms <- names(mrsldgList)
    edgesList <- lapply(seq_along(mrsldgList), function(f){
                               if(names(mrsldgList)[f] %in% nodenms){
                                 edg <- as.vector(unlist(mrsldgList[[f]]))
                                 edg <- edg[edg %in% de.results$external_gene_name]
                                 return(data.frame(rep(names(mrsldgList)[f],length=length(edg)),edg))
                               }
                             })
    edges <- do.call(rbind,edgesList)
    edgeGenes <- sort(unique(as.vector(unlist(edges))))
    resEdges <- de.results %>% dplyr::filter(external_gene_name %in% edgeGenes) %>%
                            dplyr::arrange(external_gene_name) %>%
                            dplyr::mutate(absFC = 2^log2FoldChange)
    absFC <- resEdges$absFC
    sumabsFC <- summary(absFC)
    for(x in 1:length(absFC)){
      if(absFC[x] < c(sumabsFC[2][[1]] - abs(sumabsFC[3][[1]]))) {
        absFC[x] <- sumabsFC[2][[1]] - abs(sumabsFC[3][[1]])
      }
      if(absFC[x] > c(sumabsFC[5][[1]] + abs(sumabsFC[3][[1]]))){
        absFC[x] <- sumabsFC[5][[1]] + abs(sumabsFC[3][[1]])
      }
    }
    absFCCol <- colorRampPalette(colors=c("red", "white", "forestgreen"), space="rgb")(length(absFC))
    absFCrgb <- cbind(t(col2rgb(absFCCol)), rep("1",length(absFC)))

    numConnects <- data.frame(gene=names(table(as.vector(unlist(edges)))),
                              size=as.vector(log2(table(as.vector(unlist(edges))))+.1)*3)
    resnodesVizAtt <- list(color=absFCrgb,
                           size=numConnects$size,
                           shape=rep("disk",length(absFC)))
    colnames(edges) <- c("X1","X2")

    nodes <- data.frame(matrix(c(edgeGenes, edgeGenes),ncol=2))
    resNodesAtt <- data.frame(X1=edgeGenes, Gene=edgeGenes, Size=numConnects$size)

    ## gexf
    if(minconx != 0){

      ##pick nodes to include by minimum connections
      whichResnodes <- resnodesVizAtt[[2]] > minconx
      resnodesVizAttMC <- list(color=absFCrgb[whichResnodes,],
                               size=numConnects$size[whichResnodes],
                               shape=rep("disk",length(absFC[whichResnodes])))
      colnames(edges) <- c("X1","X2")

      nodesMC <- data.frame(matrix(c(edgeGenes[whichResnodes], edgeGenes[whichResnodes]),ncol=2))
      resNodesAttMC <- data.frame(X1=edgeGenes[whichResnodes], Gene=edgeGenes[whichResnodes], Size=numConnects$size[whichResnodes])
      edgesMC <- edges[edges$X2 %in% levels(unlist(nodesMC)),]

      write.gexf(nodesMC,
                 edgesMC,
                 nodesAtt=resNodesAttMC,
                 nodesVizAtt=resnodesVizAttMC,
                 output=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/", TAG, ".all.results.p", qvals[pp], ".minconx_", minconx, ".gexf"))
    }
    if(minconx == 0){
      write.gexf(nodes,
                 edges,
                 nodesAtt=resNodesAtt,
                 nodesVizAtt=resnodesVizAtt,
                 output=paste0(OUTDIR, "/ARACNe-AP/", TAG, "/", TAG, ".results.p", qvals[pp], ".gexf"))
    }
  }

  retList <- list(mrs, mrstab)
  print(paste0("Saving: ", TAG))
  rm(c(nullmodel, rownamesp, mrstab))

  save(list=ls(envir=environment(), all.names = TRUE),
            file=paste0(OUTDIR, "/RData/msViper.", TAG, ".RData"),
            envir=environment())
  return(retList)
}
