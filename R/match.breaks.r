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
#' 
#' ## Obtain breakpoints from segmentation data
#' seg <- validate.seg(segdat_lung_ccle)
#' brk1 <- seg.breaks(seg)
#' 
#' ## Obtain breakpoints from SV calls data
#' sv <- validate.sv(svdat_lung_ccle)
#' brk2 <- sv.breaks(sv)
#' 
#' common.brk <- match.breaks(brk1, brk2)
#' 
#' ## average percentage of colocalizing breaks
#' restab <- data.frame(common.breaks$restab)[order(common.breaks$restab$total.brk2),]
#' m2 <- sprintf("%.1f",100*mean(restab$matched.brk2/restab$total.brk2)) 
#' 
#' ## Plot the proportion of SV breakpoints that have colocalizing CNV breakpoints
#' barplot(rbind(restab$matched.brk2, restab$total.brk2 - restab$matched.brk2),
#'         border=NA,las=2,xlab="",horiz=FALSE,cex.main=.7,cex.names=.4, names=rownames(restab))
#' legend("top",paste("SV breaks matched by CNV breaks\n","Average = ",m2,"%",sep=""),bty='n')
#' grid(ny=NULL,nx=NA)




match.breaks <- function(brk1, 
                         brk2, 
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

