#' 
#'
#' Obtain a matrix with the weighted average CN per chromosome arm 
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param n.reps (numeric) number of samples with identical segment to consider artifact
#' @param cnv.size (numeric) only smaller segments will be modified in the cnv data.frame
#' @param pc.overlap (numeric) minimun percentage overlap for a pair of segments to be consider identical 
#' @param fill (logical) whether to fill gaps from the segmentaed file after filtering artifacts
#' @param verbose  (logical)
#' @keywords CNV, segmentation, filter
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' clean.cnv.artifact(cnv)

clean.cnv.artifact<- function(cnv,
                              n.reps=4,
                              cnv.size=5000000,
                              pc.overlap=0.99,
                              fill.gaps=TRUE,
                              verbose=TRUE){

  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)
  require(tidyr,quietly = TRUE,warn.conflicts = FALSE)
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)

      
cnvdat <- validate.cnv(cnv)

all_artifacts_l <-list()

cnvdat_short <- cnvdat[which(cnvdat$end - cnvdat$start < cnv.size),]

for(chr in unique(cnvdat$chrom)){

  if(verbose) cat("\r",chr)

  segchr <- cnvdat_short[which(cnvdat_short$chrom == chr),]
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
allsegids <- unite(cnvdat, newcol, c(sample,chrom,start,end), remove=FALSE,sep=":")$newcol
cnvdat_clean <- cnvdat[which(!allsegids %in% toremove),]

if(fill.gaps){ 
    segclean_fill <-  segment.gap(cnvdat_clean, verbose=verbose)
    return(segclean_fill)
  }else{
    return(cnvdat_clean)
  }

}


#' 
#' 
#' Fills the gaps in a segmentation data.frame. Chromosome limits are defined for the complete segmentation dataset then segments fill the missing terminal regions. 
#' The CN log-ratio of the added segments is set to the average of the closest neighbours in each sample.
#' 
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @param verbose (logical)
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' segment.gap()

segment.gap <- function(cnv,
                        chrlist=NULL,
                        verbose=TRUE){
    
    cnvdat <- validate.cnv(cnv)
    if(is.null(chrlist)) chrlist <- unique(cnvdat$chrom)
    
    chrlims <- chromosome.limit.coords(cnvdat)
    
    cnvdat <- cnvdat[order(cnvdat$start),]
    cnvdat <- cnvdat[order(match(cnvdat$chrom, chrlist)),]
    cnvdat <- cnvdat[order(cnvdat$sample),]
    
    
    if(verbose){
        message("Filling gaps is the segmentation data.frame")
        pb <- txtProgressBar(style=3)
        cc <-0
        tot <- nrow(cnvdat)
    }
    newsegments<-list()
    if(cnvdat[1,"start"] > chrlims[cnvdat[1,"chrom"],"begin"]){ 
        newsegments[["1"]] <- data.frame(cnvdat[1,c("sample","chrom")],chrlims[cnvdat[1,"chrom"],"begin"],cnvdat[1,"start"]-1,0,cnvdat[1,"segmean"])
    }
    
    for(i in 2:nrow(cnvdat)){
        if(cnvdat[i,"chrom"] == cnvdat[i-1,"chrom"] ){
            if( cnvdat[i,"start"] - cnvdat[i-1,"end"] > 100){
                newsegments[[as.character(i)]] <- data.frame(cnvdat[i,c("sample","chrom")],cnvdat[i-1,"end"]+1,cnvdat[i,"start"]-1,0,mean(cnvdat[c(i,i-1),"segmean"]) )
            }
        }else{
            if(cnvdat[i,"start"] > chrlims[cnvdat[i,"chrom"],"begin"]){ 
                newsegments[[as.character(i)]] <- data.frame(cnvdat[i,c("sample","chrom")],chrlims[cnvdat[i,"chrom"],"begin"],cnvdat[i,"start"]-1,0,cnvdat[i,"segmean"])
            }
            if(cnvdat[i-1,"end"] < chrlims[cnvdat[i-1,"chrom"],"end"]){
                newsegments[[as.character(i)]] <- data.frame(cnvdat[i-1,c("sample","chrom")],cnvdat[i-1,"end"]+1,chrlims[cnvdat[i-1,"chrom"],"end"],0,cnvdat[i,"segmean"])
            }
        }
        if(verbose) cc <- cc+1
        if(verbose) setTxtProgressBar(pb, cc/tot)
    }
    if(cnvdat[i,"end"] < chrlims[cnvdat[i,"chrom"],"end"]) newsegments[[as.character(i)]] <- data.frame(cnvdat[i,c("sample","chrom")],cnvdat[i,"end"]+1,chrlims[cnvdat[i,"chrom"],"end"],0,cnvdat[i-1,"segmean"])
    if(verbose) close(pb)
    
    newsegments<- lapply(newsegments, setNames, colnames(cnvdat))
    
    segout <- rbind(cnvdat, do.call(rbind,newsegments))
    segout <- segout[order(segout$start),]
    segout <- segout[order(segout$sample),]
    
    return(segout)
}

