#' 
#'
#' Identify read-depth breakpoints provided a segmentation file
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param min.seg.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @keywords CNV, segmentation
#' @export
#' @examples
#' read.depth.breaks()


seg.breaks <- function(seg,
                      fc.pct = 0.2,
                      min.seg.size = NULL,
                      min.num.probes=NULL,
                      low.cov=NULL,
                      verbose=TRUE){

  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  
  segdat <- validate.seg(seg)
  if(!is.null(min.seg.size)) segdat <- segdat[which(segdat$end - segdat$start >= min.seg.size),]
  if(!is.null(min.num.probes)) segdat <- segdat[which(segdat$probes  >= min.num.probes),]
  if(verbose){
    message("Calculating breakpoints")
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(unique(segdat[,1]))
  }
  
  cn_breaks <- list()
  for(i in unique(segdat$sample)){
    sample_id <- i
    segdata_i <- segdat[which(segdat$sample == i),]
    for(chr in paste("chr",c(1:22,"X"),sep="") ){
      segdata_i_chr <- segdata_i[which(segdata_i[,2] == chr),]
      if(nrow(segdata_i_chr) > 1){
        delta_all <-  (2^segdata_i_chr[1:nrow(segdata_i_chr)-1,6]) / (2^segdata_i_chr[2:nrow(segdata_i_chr),6])
        break_idx <- c(which( log2(delta_all) >= log2(1+fc.pct)),which( log2(delta_all) < log2(1 - fc.pct)))
        brkpos1 <- segdata_i_chr[break_idx , 4]
        brkpos2 <- segdata_i_chr[break_idx + 1, 3]
        delta <- delta_all[break_idx]
        if(length(brkpos1) > 1) {
          cn_breaks[[paste(i,chr)]]<- cbind(rep(sample_id,length(brkpos1)),rep(chr,length(brkpos1)),brkpos1,brkpos2,delta)
        }else if(length(brkpos1) == 1){
          cn_breaks[[paste(i,chr)]]<- cbind(sample_id,chr,brkpos1,brkpos2,delta)
        }
      }
    }
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)
  }
  if(verbose) close(pb)
  breakpoints <- data.frame(do.call(rbind,cn_breaks))
  colnames(breakpoints) <- c("sample","chrom","start","end","FC")
  segdataB<-breakpoints[which(!breakpoints[,"chrom"] %in% c("chrY","chrM")),]
  breakpoints[,"chrom"] <- as.character(breakpoints[,"chrom"] )
  breakpoints[,"start"] <- as.numeric(as.character(breakpoints[,"start"] ))
  breakpoints[,"end"] <- as.numeric(as.character(breakpoints[,"end"] ))
  breakpoints[,"sample"] <- as.character(breakpoints[,"sample"] )
  breakpoints[,"FC"] <- as.numeric(as.character(breakpoints[,"FC"] ))
  
  if(!is.null(low.cov)){
    message("Filtering breakpoints in low coverage regiomns")
    colnames(low.cov) <- c("chrom","start","end")
    low_cov_GR = with(low.cov, GRanges(chrom, IRanges(start=start, end=end)))
    breakpoints_GR = with(breakpoints, GRanges(chrom, IRanges(start=start, end=end)))
    overlapgr <- GenomicAlignments::findOverlaps(breakpoints_GR,low_cov_GR,ignore.strand=TRUE)
    breakpoints <- breakpoints[setdiff(1:nrow(breakpoints),queryHits(overlapgr)),]
  }
  
  return(breakpoints)
}
