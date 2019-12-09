#' 
#'
#' This function calculates the percent genome changed by CNV
#' @param seg (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation
#' @export
#' @examples
#' validate.seg()


validate.seg <- function(seg){
  
  require(taRifx)  # contains remove.factors
  
  stopifnot(ncol(seg) >= 6)
  seg <- remove.factors(data.frame(seg[,1:6]))
  colnames(seg) <- c("sample","chrom","start","end","probes","segmean")
  if(length(grep("chr",seg[1,2])) == 0) seg[,"chrom"] <- paste("chr",seg$chrom,sep="")
  stopifnot(is.numeric(seg$start))
  stopifnot(is.numeric(seg$end))
  stopifnot(is.numeric(seg$segmean))
  stopifnot(is.character(seg$sample))
  stopifnot(is.character(seg$chrom))
  
  seg<- seg[order(seg$start),]
  options(warn=-1)
  seg<- seg[order(as.numeric(gsub("chr","",seg$chrom))),]
  options(warn=0)
  seg<- seg[order(seg$sample),]
  
  return(seg)
  
}
