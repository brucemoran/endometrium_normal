##use kallisto output as input for pipe
##this will massibely reduce footprint of memory and disk usage
options(scipen=999)

##libraries
libs <- c("tximport", "tidyverse")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = T)))})

##input parameters, load clinical and functions
BASEDIR <- "/data/genome/projects/david_crosby/MMMT_TCGA-UCS_SRP104165_SRP105769_DCRNASEQ_ENS"
CLINDIR <- paste0(BASEDIR, "/data/clinical")
DATADIR <- paste0(BASEDIR, "/analysis/RNAseq")
TX2GENE <- paste0(CLINDIR, "/Homo_sapiens.GRCh37.cdna.all.fa.tx2gene")

##based on: https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#use-with-downstream-bioconductor-dge-packages
tx2gene <- read.table(TX2GENE, header=TRUE)
files <- grep("abundance.tsv", dir(DATADIR, recursive=TRUE, full.names=TRUE), value=TRUE)
names(files) <- unlist(lapply(files, function(f){rev(strsplit(f,"/")[[1]])[4]}
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
nzCounts <- as_tibble(txi.kallisto.tsv$counts, rownames="ensembl_gene_id") %>%
            dplyr::mutate_if(is.numeric, round, 0) %>%
            dplyr::mutate_if(is.numeric, as.integer) %>%
            dplyr::mutate(Mean = rowMeans(as.matrix(.[grep("ensembl_gene_id",names(.), invert=TRUE)]))) %>%
            dplyr::filter(Mean > 5) %>%
            dplyr::select(-Mean) %>%
            dplyr::mutate(ensembl_gene_id = sub("\\.[^\\.]+$","", ensembl_gene_id)) %>%
            arrange(ensembl_gene_id)


##add annotation info and into a tibble
print("Reading count data and joining...")
nzCountsAnnoTr <- left_join(annoBM(attrs=c("ensembl_gene_id", "external_gene_name"),
                               filts="ensembl_gene_id",
                               arrone="ensembl_gene_id",
                               mart=mart),
                               nzCounts,
                               by = c(ensembl_gene_id = "ensembl_gene_id")) %>%
                   dplyr::group_by(ensembl_gene_id) %>%
                   ungroup() %>%
                   na.omit()

##Kallisto counts
# tx2gene <- read.table(TX2GENE, header=TRUE)
# files <- grep("abundance.tsv", dir(DATADIR, recursive=TRUE, full.names=TRUE), value=TRUE)
#
# txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
# nzCounts <- as_tibble(txi.kallisto.tsv$counts, rownames="ensembl_gene_id") %>%
#             dplyr::mutate_if(is.numeric, round, 0) %>%
#             dplyr::mutate_if(is.numeric, as.integer) %>%
#             dplyr::mutate(Mean = rowMeans(as.matrix(.[grep("ensembl_gene_id",names(.), invert=TRUE)]))) %>%
#             dplyr::filter(Mean > 5) %>%
#             dplyr::select(-Mean) %>%
#             dplyr::mutate(ensembl_gene_id = sub("\\.[^\\.]+$","", ensembl_gene_id)) %>%
#             arrange(ensembl_gene_id) %>%
#             dplyr::select("ensembl_gene_id",sampleIDsIn)

##Kallisto annotated
# print("Reading count data and joining...")
# nzCountsAnnoTr <- left_join(annoBM(attrs=c("ensembl_gene_id", "external_gene_name"),
#                                filts="ensembl_gene_id",
#                                arrone="ensembl_gene_id",
#                                mart=mart),
#                                nzCounts,
#                                by = c(ensembl_gene_id = "ensembl_gene_id")) %>%
#                    dplyr::group_by(ensembl_gene_id) %>%
#                    ungroup() %>%
#                    na.omit()
