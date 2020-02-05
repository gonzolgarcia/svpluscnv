#' Class to store breakpoint annotations in association with genomic features (e.g. gene loci)
#' @param breaks (data.frame): the breakpoint info containing data.frame, this will be occupied by the CNV segmentation data in the case of cnv.break.annot or SV for sv.break.annot. Unique random string rownames are added to the returned breaks data.frame.
#' @param burden (numeric): a vector containing the total number of breakpoints in each sample 
#' @param param (list): a list of parametres provided 
#' @return an instance of the class 'breaks' containing breakpoint and breakpoint burden information
#' @export
breaks <- setClass("breaks", representation(
                            breaks  = "data.frame",
                            burden = "numeric",
                            param = "list"
                        ))


setMethod("show","breaks",function(object){
    writeLines(paste("An object of class breaks from svcnvplus containing the following stats:
                \nNumber of samples=",length(object@burden),
                "\nTotal number of breakpoints =",nrow(object@breaks)))
})



#' Identify read-depth breakpoints provided a segmentation file
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param break.width (numeric) the maximum distance between a segment end and the subsequent segment start positions beyond which breakpoints are discarded
#' @param min.cnv.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param clean.brk (numeric) identical breakpoints across multiple samples tend to be artifacts; remove breaks > N 
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from segmentation data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' cnv.breaks(cnv)


cnv.breaks <- function(cnv,
                       fc.pct = 0.2,
                       break.width = NULL,
                       min.cnv.size = NULL,
                       min.num.probes=NULL,
                       chrlist=NULL,
                       low.cov=NULL, 
                       clean.brk=NULL,
                       verbose=TRUE){
    

    cnvdat <- validate.cnv(cnv)
    
    if(is.null(chrlist)) chrlist <- unique(cnvdat$chrom)
    chrlist <- chr.sort(chrlist)
    
    brk.burden <- rep(0,length(unique(cnvdat$sample)))
    names(brk.burden) <- unique(cnvdat$sample)
    
    if(!is.null(min.cnv.size)) cnvdat <- cnvdat[which(cnvdat$end - cnvdat$start >= min.cnv.size),]
    if(!is.null(min.num.probes)) cnvdat <- cnvdat[which(cnvdat$probes  >= min.num.probes),]
    
    lastrow <- nrow(cnvdat)
    pos <- round(apply(cbind(cnvdat[2:(lastrow),"start"], cnvdat[1:(lastrow-1),"end"]),1,mean))
    chrom <- cnvdat[2:(lastrow),"chrom"]
    sample <- cnvdat[2:(lastrow),"sample"]
    width <- cnvdat[2:(lastrow),"start"] - cnvdat[1:(lastrow-1),"end"]
    FC <-  (2^cnvdat[1:(lastrow-1),"segmean"]) / (2^cnvdat[2:lastrow,"segmean"])
    
    breakpoints <- remove.factors(data.frame(sample,chrom,pos,width,FC))
    
    break_idx <- c(which( log2(FC) >= log2(1+fc.pct)),which( log2(FC) < log2(1 - fc.pct)))
    
    samechr <- which(apply(cbind(cnvdat[1:(lastrow-1),"chrom"],cnvdat[2:(lastrow),"chrom"]),1,anyDuplicated) == 2)

    samesample <-  which(apply(cbind(cnvdat[1:(lastrow-1),"sample"],cnvdat[2:(lastrow),"sample"]),1,anyDuplicated) == 2)

    if(is.null(break.width)) break.width <- Inf
    brwidthin <- which(width < break.width)
    
    breakpoints <- breakpoints[Reduce(intersect, list(break_idx,samechr,samesample,brwidthin)),]
    
    rownames(breakpoints) <- createRandomString(nrow(breakpoints),10)
    
    if(!is.null(low.cov)){
        message("Filtering breakpoints in low coverage regiomns")
        colnames(low.cov) <- c("chrom","start","end")
        low_cov_GR = with(low.cov, GRanges(chrom, IRanges(start=start, end=end)))
        breakpoints_GR = with(breakpoints, GRanges(chrom, IRanges(start=start, end=end)))
        overlapgr <- GenomicAlignments::findOverlaps(breakpoints_GR,low_cov_GR,ignore.strand=TRUE)
        breakpoints <- breakpoints[setdiff(1:nrow(breakpoints),queryHits(overlapgr)),]
    }
    
    if(!is.null(clean.brk)){
        breakids <- unite(breakpoints[,c(2:4)],newcol)$newcol
        breakids.freq <- sort(table(breakids),decreasing=T)
        breakpoints <- breakpoints[which(breakids %in% names(which(breakids.freq < clean.brk))),]
    }
    
    brk.burden.sub <- table(breakpoints$sample)
    brk.burden[names(brk.burden.sub)] <- brk.burden.sub
    
    return(breaks(breaks=remove.factors(breakpoints),
                burden=brk.burden,
                param=list(
                    datatype="cnv",
                    fc.pct = fc.pct,
                    min.cnv.size = min.cnv.size,
                    min.num.probes=min.num.probes,
                    low.cov=low.cov, 
                    clean.brk=clean.brk
                )
                )
           )
}



#' Identify SV call breakpoints and filter low cov regions if provided
#' 
#' @param svc (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @keywords Structural variants
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from SV calls data
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' svc.breaks(svc)



svc.breaks <- function(svc,low.cov=NULL){
    
    svcdat <- validate.svc(svc)
    
    brk.burden <- rep(0,length(unique(svcdat$sample)))
    names(brk.burden) <- unique(svcdat$sample)
    
    svcdat.breaks <- data.frame(c(svcdat$sample,svcdat$sample),
                               c(svcdat$chrom1,svcdat$chrom2),
                               c(svcdat$pos1,svcdat$pos2),
                               c(svcdat$strand1,svcdat$strand2),
                               c(svcdat$svclass,svcdat$svclass),
                               c(rownames(svcdat),rownames(svcdat)))
    colnames(svcdat.breaks) <- c("sample","chrom","pos","strand","svclass","id")
    if(!is.null(low.cov)){
        low.cov.df <- data.frame(low.cov[,1:3])
        colnames(low.cov.df) <- c("chrom","start","end")
        
        svc_ranges <- with(svcdat.breaks, GRanges(chrom, IRanges(start=pos, end=pos)))
        low.cov_ranges <- with(low.cov.df, GRanges(chrom, IRanges(start=start, end=end)))
        
        low.cov_ranges = GenomicAlignments::findOverlaps(svc_ranges,low.cov_ranges)
        
        svcdat.breaks <- remove.factors(svcdat.breaks[which(!svcdat.breaks$id %in% queryHits(low.cov_ranges)),])
    }else{
        svcdat.breaks <- remove.factors(svcdat.breaks)
    }
    brk.burden.sub <- table(svcdat.breaks$sample)
    brk.burden[names(brk.burden.sub)] <- brk.burden.sub
    
    rownames(svcdat.breaks) <- createRandomString(nrow(svcdat.breaks),10)
    
    return(breaks(breaks=svcdat.breaks,
                burden=brk.burden,
                param=list(
                    datatype="cnv",
                    low.cov=low.cov
                )))
    
}




#' Generating a genomic map based on a defined bin size and sliding window and counts the number of breakpoints mapped onto each bin. This function is used internally by svcnvplus::shattered.regions and svcnvplus::shattered.regions.cnv
#' 
#' @param brk (list) This object of the class 'breaks' obtained from CNV segmentation data (svcnvplus::cnv.breaks) and Structural Variant calls (svcnvplus::svc.breaks). 
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) indicating the chromosome most distal coordinates with coverage. Also returned by the function svcnvplus::chromosome.limit.coords.
#' @param window.size (numeric) size in megabases of the genmome bin onto which breakpoints will be mapped 
#' @param slide.size (numeric) size in megabases of the sliding genomic window; if slide.size < window.size the genomic bins will overlap
#' @param chrlist (character) vector containing chromosomes to consider (e.g. "chr1", "chr2", "chr3", ...)
#' @param verbose (logical)
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' cnv <- validate.cnv(segdat_lung_ccle)
#' brk <- cnv.breaks(cnv)
#' chr.lim <- chromosome.limit.coords(cnv)
#' break.density(brk,chr.lim)


break.density <- function(brk, 
                         chr.lim, 
                         window.size = 10, 
                         slide.size=2,
                         chrlist=NULL,
                         verbose=TRUE){
  
  # make sure chr.lim format is correct
  chr.lim<- chr.lim[,1:2]
  colnames(chr.lim) <- c("begin","end")
  
  # make sure both chr.lim and breaks have same chromosome names 
  seqnames <- intersect(rownames(chr.lim),brk@breaks$chr)
  stopifnot(length(seqnames) > 0) 
  
  # a template vector to save breakpoint counts 
  templatevector <- brk@burden
  templatevector[]<-0
  
  WS <- window.size * 1e+6
  SS <- slide.size * 1e+6
  offset <- window.size/slide.size
  
  if(is.null(chrlist)) chrlist <- rownames(chr.lim)
  
  # count breaks for each chromosome for each fragment
  fragment <- list()
  for(chr in  chrlist){

    if(verbose) cat("\r",chr)

    chr_breaks <- brk@breaks[which(brk@breaks$chrom == chr),]
    frag <- seq(chr.lim[chr,"begin"],chr.lim[chr,"end"]+SS,SS)
    for(i in (1+offset):length(frag)){
      start <- frag[i - offset]
      stop <- frag[i]
      fragment[[paste(chr,start,stop)]] <- templatevector
      break.position <- chr_breaks$pos
      res_bp <- table(chr_breaks[intersect(which(break.position > start),which(break.position < stop)),"sample"])
      fragment[[paste(chr,start,stop)]][names(res_bp)] <- res_bp
    }
  }
  if(verbose) cat(" Done!\n")

  return( do.call(cbind,fragment))
  
}




#' Identify CNV breakpoints provided a segmentation file
#' 
#' @param brk1 (S4) an object of class breaks as returned by `svc.breaks` and `cnv.breaks`
#' @param brk2 (S4) an object of class breaks as returned by `svc.breaks` and `cnv.breaks` to compare against brk1
#' @param maxgap (numeric) distance (base pairs) limit for nreakpoints to be consider colocalized 
#' @param verbose (logical) 
#' @keywords CNV, SV, genomic breakpoints
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from segmentation data
#' cnv <- validate.cnv(segdat_lung_ccle)
#' brk1 <- cnv.breaks(cnv)
#' 
#' ## Obtain breakpoints from SV calls data
#' sv <- validate.svc(svdat_lung_ccle)
#' brk2 <- svc.breaks(svc)
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
                         maxgap=100000,
                         verbose=FALSE){
    
    common_samples <- intersect(names(brk1@burden),names(brk2@burden))
    stopifnot(length(common_samples) > 0, local = TRUE) 
    
    brk1_match <- brk2_match <- restab <- list()
    for(id in common_samples){
        
        brk1_i <- brk1@breaks[which(brk1@breaks$sample == id),]
        brk_ranges1 <- with(brk1_i, GRanges(chrom, IRanges(start=pos, end=pos)))
        
        brk2_i <- brk2@breaks[which(brk2@breaks$sample == id),]
        brk_ranges2 <- with(brk2_i, GRanges(chrom, IRanges(start=pos, end=pos)))
        
        
        options(warn=-1)
        seg_seg = GenomicAlignments::findOverlaps(brk_ranges1, brk_ranges2, maxgap=maxgap)
        options(warn=0)
        
        brk_match1 <- sort(unique(queryHits(seg_seg)))
        brk_match2 <- sort(unique(subjectHits(seg_seg)))
        
        restab[[id]] <- c(length(brk_match1), nrow(brk1_i), length(brk_match2), nrow(brk2_i))
        names(restab[[id]]) <- c("matched.brk1", "total.brk1", "matched.brk2", "total.brk2")
        
        brk1_match[[id]] <- brk1_i[brk_match1,]
        brk2_match[[id]] <- brk2_i[brk_match2,]
    }

    restab <- data.frame(do.call(rbind,restab))
    return(list(
        brk1_match = do.call(rbind,brk1_match),
        brk2_match = do.call(rbind,brk2_match),
        restab= restab))
}

