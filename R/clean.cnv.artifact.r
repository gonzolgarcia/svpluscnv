#' 
#'
#' Obtain a matrix with the weighted average CN per chromosome arm 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param n.reps (numeric) number of samples with identical segment to consider artifact
#' @param seg.size (numeric) only smaller segments will be modified in the seg data.frame
#' @param pc.overlap (numeric) minimun percentage overlap for a pair of segments to be consider identical 
#' @param fill (logical) whether to fill gaps from the segmentaed file after filtering artifacts
#' @param verbose  (logical)
#' @keywords CNV, segmentation, filter
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' clean.cnv.artifact(seg)

clean.cnv.artifact<- function(seg,
                              n.reps=4,
                              seg.size=5000000,
                              pc.overlap=0.99,
                              fill=TRUE,
                              verbose=TRUE){

  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)
  require(tidyr,quietly = TRUE,warn.conflicts = FALSE)
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)

      
segdat <- validate.seg(seg)

all_artifacts_l <-list()

segdat_short <- segdat[which(segdat$end - segdat$start < seg.size),]

for(chr in unique(segdat$chrom)){

  if(verbose) cat("\r",chr)

  segchr <- segdat_short[which(segdat_short$chrom == chr),]
  segchr.gr <- with(segchr, GRanges(chrom, IRanges(start=start, end=end)))
  hits = GenomicAlignments::findOverlaps(segchr.gr,segchr.gr)
  overlaps <- pintersect(segchr.gr[queryHits(hits)], segchr.gr[subjectHits(hits)])
  
  percentOverlapA <- width(overlaps) / width(segchr.gr[queryHits(hits)])
  percentOverlapB <- width(overlaps) / width(segchr.gr[subjectHits(hits)])
  hits_p <- as.data.frame(hits[intersect(which(percentOverlapA >= pc.overlap),which(percentOverlapB >= pc.overlap)),])
  reps <- aggregate(subjectHits~queryHits,hits_p,paste,simplify=FALSE)
  reps_list <- reps$subjectHits
  names(reps_list) <- reps$queryHits
  reps_list_collapse <- lapply(lapply(reps_list,sort),paste,collapse=" ")
  groups_a <- table(unlist(reps_list_collapse))
  all_artifacts <- as.numeric(unlist(strsplit(names(which(groups_a > n.reps))," ")))
  all_artifacts_l[[chr]] <- segchr[all_artifacts,]
}

all_artifacts <- do.call(rbind,unname(all_artifacts_l))
toremove <- unite(all_artifacts, newcol, c(sample,chrom,start,end), remove=FALSE,sep=":")$newcol
allsegids <- unite(segdat, newcol, c(sample,chrom,start,end), remove=FALSE,sep=":")$newcol
segdat_clean <- segdat[which(!allsegids %in% toremove),]

if(fill){ 
    segclean_fill <-  segment.gap(segdat_clean, verbose=verbose)
    return(segclean_fill)
  }else{
    return(segdat_clean)
  }

}
