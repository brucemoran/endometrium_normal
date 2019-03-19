##making classifier
library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)

allreads_https  <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE97929&format=file&file=GSE97929%5Freads%5Fall%2Etxt%2Egz"
GET(allreads_https, write_disk(tf <- tempfile("")))
allreads <- as.tibble(read.table(tf, header=T))
unlink(tf)
ss.allreads <- allreads
colnames(ss.allreads) <- c("external_gene_name",
                         as.character(unlist(lapply(colnames(ss.allreads)[2:dim(ss.allreads)[2]], function(f){
                           strsp <- strsplit(f,"\\.")[[1]][3]
                           strsp <- gsub("_1", "_2", strsp)
                           strsp <- paste0(strsp,"_1")
                           strsp <- gsub("_2_1", "_2", strsp)
                           return(strsp)
                         }))))
#remove zero-sum
ss.allreads <- ss.allreads %>%
               dplyr::mutate(Median = rowMedians(as.matrix(.[grep("exter",names(.), invert=TRUE)]))) %>%
               dplyr::filter(Median > 0) %>%
               dplyr::select(-Median) %>%
               na.omit() %>%
               dplyr::select(-"NA_1") %>%
               as.data.frame()
rownames(ss.allreads) <- ss.allreads[,1]
ss.allreads <- ss.allreads[,-1]

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
ss.allreads <- ss.allreads %>% dplyr::select(ss.cond$sampleID)
ss.cond <- ss.cond %>% dplyr::filter(sampleID %in% colnames(ss.allreads))
dds.allreads <- DESeqDataSetFromMatrix(countData=ss.allreads ,
                                colData=ss.cond,
                                design=~1)
ss.allreads.vst <- assay(vstReturnPCAPlot(dds.allreads,
                                intgroup=c("Batch", "Group"),
                                PLOTDIR=PLOTDIR,
                                TAG="Group-Batch.GSE97929"))

# from https://www.biostars.org/p/86981/
ffun <- filterfun(pOverA(p = 0.2, A = 5), cv(a = 0.7, b = 10))
ss.allreads.vst.filt <- genefilter(2^ss.allreads.vst,ffun)
ss.allreads.filt <- ss.allreads[rownames(ss.allreads) %in% grep("RNA_SPIKE_ERCC", names(grep("TRUE",ss.allreads.vst.filt, value=T)), invert=TRUE, value=TRUE),]

##now remove those filtered, back into vstPCA
dds.allreads.filt <- DESeqDataSetFromMatrix(countData=ss.allreads.filt,
                                colData=ss.cond,
                                design=~1)
ss.allreads.filt.vst <- assay(vstReturnPCAPlot(dds.allreads.filt,
                                intgroup=c("Batch", "Group"),
                                PLOTDIR=PLOTDIR,
                                TAG="Group-Batch.GSE97929"))
##some negative values, move to 0
ss.allreads.filt.vst[ss.allreads.filt.vst < 0] <- 0
ss.allreads.filt.vst.eset <- createEset(as.data.frame(ss.allreads.filt.vst), ss.cond)
##https://cran.r-project.org/web/packages/NNLM/NNLM.pdf
res <- nmf(ss.allreads.filt.vst.eset , rank=4)
summary(res, class=ss.allreads.filt.vst.eset$Group)
w <- basis(res)
h <- coef(res)
s <- featureScore(res)

pdf("Gene.nnmf.heatmap.pdf")
heatmap(decomp$W, Colv = NA, xlab ='Meta-gene', ylab ='Gene', margins = c(2,2),labRow ='', labCol ='', scale ='column', col = cm.colors(100));
dev.off()
pdf("Patient.nnmf.heatmap.pdf")
heatmap(decomp$H, Rowv = NA, ylab ='Meta-gene', xlab ='Patient', margins = c(2,2),labRow ='', labCol ='', scale ='row', col = cm.colors(100));
dev.off()
target <- factor(ss.cond$Group)
names(target) <- ss.cond$sampleID

tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)
rf_output <- randomForest(x=predictor_data, y=target, importance=TRUE, ntree=10001, proximity=TRUE, sampsize=sampsizes)
rf_importances <- importance(rf_output, scale=FALSE)
confusion <- rf_output$confusion
sensitivity <- (confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity <- (confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error <- rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy <- 1-overall_error
class1_error <- paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error <- paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy <- 100-overall_error
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
pdf("MDS_pdffile.pdf")
target_labels <- as.vector(target)
target_labels[target_labels=="early_epithelium"]="E"
target_labels[target_labels=="mid_epithelium"]="M"
target_labels[target_labels=="early_stroma"]="S"
target_labels[target_labels=="mid_stroma"]="M"
MDSplot(rf_output, target, k=2, pch=target_labels, palette=c("lightblue", "dodgerblue", "lightgreen", "forestgreen"), main="MDS plot")
dev.off()
