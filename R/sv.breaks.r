#' 
#'
#' Identify SV call breakpoints and filter low cov regions if provided
#' @param sv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from SV calls data
#' sv <- validate.sv(svdat_lung_ccle)
#' 
#' sv.breaks(sv)



sv.breaks <- function(sv,low.cov=NULL){
  
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  
  svdat <- validate.sv(sv)
  
  brk.burden <- rep(0,length(unique(svdat$sample)))
  names(brk.burden) <- unique(svdat$sample)

  svdat.breaks <- data.frame(c(svdat$sample,svdat$sample),
                             c(svdat$chrom1,svdat$chrom2),
                             c(svdat$pos1,svdat$pos2),
                             c(svdat$pos1,svdat$pos2),
                             c(svdat$strand1,svdat$strand2),
                             c(svdat$svclass,svdat$svclass),
                             c(rownames(svdat),rownames(svdat)))
  colnames(svdat.breaks) <- c("sample","chrom","start","end","strand","svclass","id")
  if(!is.null(low.cov)){
    low.cov.df <- data.frame(low.cov[,1:3])
    colnames(low.cov.df) <- c("chrom","start","end")
  
    sv_ranges <- with(svdat.breaks, GRanges(chrom, IRanges(start=start, end=end)))
    low.cov_ranges <- with(low.cov.df, GRanges(chrom, IRanges(start=start, end=end)))
  
    low.cov_ranges = GenomicAlignments::findOverlaps(sv_ranges,low.cov_ranges)
    
    svdat.breaks <- remove.factors(svdat.breaks[which(!svdat.breaks$id %in% queryHits(low.cov_ranges)),])
  }else{
    svdat.breaks <- remove.factors(svdat.breaks)
  }
  brk.burden.sub <- table(svdat.breaks$sample)
  brk.burden[names(brk.burden.sub)] <- brk.burden.sub
  
  return(list(breaks=svdat.breaks,
              brk.burden=brk.burden))
  
}


