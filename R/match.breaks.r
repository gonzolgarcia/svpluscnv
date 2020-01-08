#' 
#'
#' Identify CNV breakpoints provided a segmentation file
#' @param brk1 (data.frame) brakpoint coordinates 1 (sample, chrom, pos), as returned by `sv.breaks` and `seg.breaks`
#' @param brk2 (data.frame) brakpoint coordinates 1 (sample, chrom, pos), as returned by `sv.breaks` and `seg.breaks`
#' @param maxgap (numeric) distance (base pairs) limit for nreakpoints to be consider colocalized 
#' @param verbose (logical) 
#' @keywords CNV, SV, genomic breakpoints
#' @export
#' @examples
#' match.breaks()



match.breaks <- function(brk1, brk2, 
                          maxgap=50000,
                          verbose=FALSE){
  
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)
 
  common_samples <- intersect(names(brk1$brk.burden),names(brk2$brk.burden))
  stopifnot(length(common_samples) > 0, local = TRUE) 
  
  if(verbose){
    message(paste("Finding common breaks for",length(common_samples),"common samples",sep=" "))
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(common_samples)
  }
  
  brk1_match <- brk2_match <- restab <- list()
  for(id in common_samples){
    
    brk1_i <- brk1$breaks[which(brk1$breaks$sample == id),]
    brk_ranges1 <- with(brk1_i, GRanges(chrom, IRanges(start=start, end=end)))

    brk2_i <- brk2$breaks[which(brk2$breaks$sample == id),]
    brk_ranges2 <- with(brk2_i, GRanges(chrom, IRanges(start=start, end=end)))
    

    options(warn=-1)
    seg_seg = GenomicAlignments::findOverlaps(brk_ranges1,brk_ranges2,maxgap=maxgap)
    options(warn=0)
    
    brk_match1 <- sort(unique(queryHits(seg_seg)))
    brk_match2 <- sort(unique(subjectHits(seg_seg)))
    
    restab[[id]] <- c(length(brk_match1), nrow(brk1_i), length(brk_match2), nrow(brk2_i))
    names(restab[[id]]) <- c("matched.brk1", "total.brk1", "matched.brk2", "total.brk2")
    
    brk1_match[[id]] <- brk1_i[brk_match1,]
    brk2_match[[id]] <- brk2_i[brk_match2,]
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)
  }
  if(verbose) close(pb)
  
  restab <- data.frame(do.call(rbind,restab))
  return(list(
    brk1_match = do.call(rbind,brk1_match),
    brk2_match = do.call(rbind,brk2_match),
    restab= restab))
}

