#' 
#'
#' Identify CNV breakpoints provided a segmentation file
#' @param seg (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param segment.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param verbose (logical) 
#' @keywords CNV, segmentation
#' @export
#' @examples
#' break.density()



match.variant.breaks <- function(seg, sv, 
                              maxgap=50000,
                              fc.pct = 0.2,
                              low.cov=NULL,
                              clean.brk=NULL,
                              verbose=FALSE){

  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)

  segdat <- validate.seg(seg)
  segdat.breaks <- seg.breaks(segdat, fc.pct = fc.pct, low.cov= low.cov, clean.brk=clean.brk, verbose=verbose) 
  
  svdat <- validate.sv(sv)
  svdat.breaks <- sv.breaks(svdat,low.cov=low.cov)
  
  
  common_samples <- intersect(svdat$sample,segdat$sample)
  
  stopifnot(length(common_samples) > 0, local = TRUE) 
  if(verbose){
    message(paste("Finding common breaks for",length(common_samples),"common samples",sep=" "))
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(unique(segdat[,1]))
  }
  
  sv_results <- seg_results <- restab <- list()
  for(id in common_samples){
    
    svdat.breaks_i <- data.frame(svdat.breaks[which(svdat.breaks$sample == id),])
    sv_ranges <- with(svdat.breaks_i, GRanges(chrom, IRanges(start=start, end=end)))
    
    segdat.breaks_i <- data.frame(segdat.breaks[which(segdat.breaks$sample == id),])
    seg_ranges <- with(segdat.breaks_i, GRanges(chrom, IRanges(start=start, end=end)))

    sv_seg = GenomicAlignments::findOverlaps(sv_ranges,seg_ranges,maxgap=maxgap)

    sv_match <- sort(unique(queryHits(sv_seg)))
    seg_match <- sort(unique(subjectHits(sv_seg)))
   
    restab[[id]] <- c(length(sv_match), nrow(svdat.breaks_i), length(seg_match), nrow(segdat.breaks_i))
    names(restab[[id]]) <- c("matched.sv", "total.sv", "matched.seg", "total.seg")
    
    sv_results[[id]] <- sv[svdat.breaks_i[sv_match,"id"],]
    seg_results[[id]] <- segdat.breaks_i[seg_match,]
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)
  }
  if(verbose) close(pb)
  
  restab <- data.frame(do.call(rbind,restab))

  return(list(
    sv_validated = do.call(rbind,sv_results),
    seg_validated = do.call(rbind,seg_results),
    restab= restab))
}

