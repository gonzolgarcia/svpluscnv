#' 
#'
#' Identify SV call breakpoints and filter low cov regions if provided
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @keywords CNV, segmentation
#' @export
#' @examples
#' read.depth.breaks()



sv.breaks <- function(sv,low.cov=NULL){
  
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  
  svdat <- validate.sv(sv)
  svdat.breaks <- data.frame(c(svdat$sample,svdat$sample),c(svdat$chrom1,svdat$chrom2),c(svdat$pos1,svdat$pos2),c(1:nrow(svdat),1:nrow(svdat)))
  colnames(svdat.breaks) <- c("sample","chrom","pos","id")
  if(!is.null(low.cov)){
    low.cov.df <- data.frame(low.cov[,1:3])
    colnames(low.cov.df) <- c("chrom","start","end")
  
    sv_ranges <- with(svdat.breaks, GRanges(chrom, IRanges(start=pos, end=pos)))
    low.cov_ranges <- with(low.cov.df, GRanges(chrom, IRanges(start=start, end=end)))
  
    low.cov_ranges = GenomicAlignments::findOverlaps(sv_ranges,low.cov_ranges)
    
    return(svdat.breaks[which(!svdat.breaks$id %in% queryHits(low.cov_ranges)),])
  }else{
    return(svdat.breaks)
  }
}
