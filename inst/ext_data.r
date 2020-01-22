#' 
#' load RefSeq transcript data fro hg19 and hg38
#' @title refseq data type
#' @details A refseq object is created containing transcript annotations
#' including exons coordinates
#' 

require(data.table)  # contains fread
require(taRifx)  # contains remove.factors



setwd("svcnvplus/")

refseq_hg19_df <- as.data.frame(remove.factors(fread("inst/extdata//NCBI_RefSeq_Feb2009_GRCh37-hg19_UCSC_TableBrowser.csv.gz")))
refseq_hg38_df <- as.data.frame(remove.factors(fread("inst/extdata//NCBI_RefSeq_Feb2013_GRCh38-hg38_UCSC_TableBrowser.csv.gz")))

exonStarts_hg19 = lapply(unname(sapply(refseq_hg19_df$exonStarts, function(x)  strsplit(x,","))),as.numeric)
exonEnds_hg19 = lapply(unname(sapply(refseq_hg19_df$exonEnds, function(x)  strsplit(x,","))),as.numeric)
exonFrames_hg19 = lapply(unname(sapply(refseq_hg19_df$exonFrames, function(x)  strsplit(x,","))),as.numeric)
df_hg19 <- refseq_hg19_df[,setdiff(colnames(refseq_hg19_df),c("exonStarts","exonEnds","exonFrames"))]

names(exonStarts_hg19) <-names(exonEnds_hg19) <-names(exonFrames_hg19) <- refseq_hg19_df$name


refseq_hg19<-list(df=df_hg19,
                  exonStarts = exonStarts_hg19,
                  exonEnds = exonEnds_hg19,
                  exonFrames = exonFrames_hg19)


exonStarts_hg38 = lapply(unname(sapply(refseq_hg38_df$exonStarts, function(x)  strsplit(x,","))),as.numeric)
exonEnds_hg38 = lapply(unname(sapply(refseq_hg38_df$exonEnds, function(x)  strsplit(x,","))),as.numeric)
exonFrames_hg38 = lapply(unname(sapply(refseq_hg38_df$exonFrames, function(x)  strsplit(x,","))),as.numeric)
df_hg38 <- refseq_hg38_df[,setdiff(colnames(refseq_hg38_df),c("exonStarts","exonEnds","exonFrames"))]

names(exonStarts_hg38) <-names(exonEnds_hg38) <-names(exonFrames_hg38) <-refseq_hg38_df$name

refseq_hg38 <- list(df=df_hg38,
                  exonStarts = exonStarts_hg38,
                  exonEnds = exonEnds_hg38,
                  exonFrames = exonFrames_hg38)

save(refseq_hg38, refseq_hg19, file="data/refseq.rda")

## end
#refseq_hg38 <-list(
#    symbol = as.character(refseq_hg38_df$name2),
#    name = as.character(refseq_hg38_df$name),
#    chrom = as.character(refseq_hg38_df$chrom),
#    strand = as.character(refseq_hg38_df$strand),
#    txStart = as.numeric(refseq_hg38_df$txStart),
#    txEnd = as.numeric(refseq_hg38_df$txEnd),
#    cdsStart = as.numeric(refseq_hg38_df$cdsStart),
#    cdsEnd = as.numeric(refseq_hg38_df$cdsEnd),
#    exonCount = as.numeric(refseq_hg38_df$exonCount),
#    exonStarts = exonStarts_hg38,
#    exonEnds = exonEnds_hg38,
#    score = as.numeric(refseq_hg38_df$score),
#    exonFrames = exonFrames_hg38,
#    numtx=nrow(refseq_hg38_df)
#)
#
