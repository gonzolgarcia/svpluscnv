#' 
#'
#' Caller for shattered genomic regions based on breakpoint densities
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents 20 percent fold change
#' @param min.seg.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @param num.breaks (numeric) size in megabases of the genmome bin to compute break density 
#' @param num.sd (numeric) size in megabases of the sliding genmome window
#' @param dist.iqm.cut (numeric) interquantile average of the distance between breakpoints within a shattered region
#' @param verbose (logical)
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' shattered.regions.cnv(seg)

shattered.regions.cnv <- function(seg,
                              fc.pct = 0.2,
                              min.seg.size = 0,
                              min.num.probes=0, 
                              low.cov = NULL,
                              clean.brk=NULL,
                              window.size = 10,
                              slide.size = 2,
                              num.breaks = 10,
                              num.sd = 5,
                              dist.iqm.cut = 1e+05,
                              verbose=FALSE
                              ){
  
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  
  segdat <- validate.seg(seg)
  
  chr.lim <- chromosome.limit.coords(segdat)
  
  breaks <- seg.breaks(seg = segdat, 
                       fc.pct = fc.pct, 
                       min.seg.size = min.seg.size, 
                       low.cov = low.cov, 
                       clean.brk=clean.brk,
                       verbose = verbose)
  
  if(verbose) message("Mapping CNV breakpoints across the genome:")
  seg.brk.dens <- break.density(breaks, 
                                chr.lim = chr.lim, 
                                window.size = window.size, 
                                slide.size = slide.size,
                                verbose = verbose)
  

  
  # calculate inter quantile mean and standard deviation per sample
  iqmdata1<- sddata<- breaks$brk.burden
  iqmdata1[] <- sddata[] <- 0
  
  iqmdata <- apply(seg.brk.dens,1,IQM,lowQ=0.2,upQ=0.8)
  sddata <- apply(seg.brk.dens,1,IQSD,lowQ=0.2,upQ=0.8)

  a <- sapply(rownames(seg.brk.dens),function(i) names(which(seg.brk.dens[i,] > iqmdata[i]+num.sd*sddata[i] )))
  b <- sapply(rownames(seg.brk.dens),function(i) names(which(seg.brk.dens[i,] >= num.breaks)))
  
  # condition for chromothripsis: at least n=breaks > 6 (sv SND seg)  AND n-breaks > u+2*sd (sv AND seg) 
  res <- sapply(rownames(seg.brk.dens),function(i) Reduce(intersect, list(b[[i]],a[[i]])) )
  highDensityRegions <- seg.brk.dens
  highDensityRegions[] <- 0
  for(cl in rownames(seg.brk.dens)) highDensityRegions[cl,res[[cl]]] <- 1
  
  res <- res[which(unlist(lapply(res,length)) >0)]

  if(verbose){
    message("Locating shattered regions by CNV only...")
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(res)
  }
  
  restab <- list()
  for(cl in names(res)){
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)

      if(length(res[[cl]]) > 0){
      
      tab <- as.data.frame(do.call(rbind,strsplit(res[[cl]]," ")))
      tab[,2]<-as.numeric(as.character(tab[,2]))
      tab[,3]<-as.numeric(as.character(tab[,3]))
      colnames(tab) <- c("chrom","start","end")
      tabgr = with(tab, GRanges(chrom, IRanges(start=start, end=end))) 
      hits = as.data.frame(GenomicAlignments::findOverlaps(tabgr,tabgr))
      
      agg <- aggregate(subjectHits ~ queryHits, hits, paste,simplify=FALSE)
      prev<-c(); cnum <- 0
      agglist <- list()
      for(x in agg$subjectHits){
        if(length(intersect(x,prev) > 0)){
          agglist[[cnum]] <- unique(c(x,prev))
          prev <- agglist[[cnum]]
        }else{
          cnum <- cnum+1
          agglist[[cnum]]<- x
          prev <-agglist[[cnum]]
        }
      }
      agglistUniq <- list()
      for(i in 1:length(agglist)){
        chr <- as.character(unique(tab[as.numeric(agglist[[i]]),"chrom"]))
        start <-min( tab[as.numeric(agglist[[i]]),"start"])
        end <- max( tab[as.numeric(agglist[[i]]),"end"])
        segNum <- length(agglist[[i]])
        agglistUniq[[i]] <-  c(chr,start,end,segNum)
      }
      tabmerged <- remove.factors(as.data.frame(do.call(rbind,agglistUniq)))
      for(i in 2:4) tabmerged[,i] <- as.numeric( tabmerged[,i] )
      colnames(tabmerged) <- c("chrom","start","end","nseg")
      restab[[cl]] <- tabmerged
    }
  }
  if(verbose) close(pb)
  
  if(verbose){
    message("Evaluating shattered regions by CNV data only...")
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(restab)
  }
  for(cl in names(restab)){
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)
    regions <-   restab[[cl]]
    br1 <- breaks$breaks[which(breaks$breaks$sample == cl),2:3]
    colnames(br1) <- c("chrom","pos")
    br1.gr <- with(br1, GRanges(chrom, IRanges(start=pos, end=pos)))
    regions_gr <- with(regions, GRanges(chrom, IRanges(start=start, end=end)))
    hits_1 = GenomicAlignments::findOverlaps(regions_gr,br1.gr)
    n.brk <- dist.iqm <- start <- end <- rep(0,nrow(regions))
    conf <- rep("HC",nrow(regions))
    for(i in 1:nrow(regions)){
      sites <- sort(unique(br1[subjectHits(hits_1)[which(queryHits(hits_1) == i)],"pos"]))
      dist.iqm[i]  <- IQM(sites[2:length(sites)] - sites[1:(length(sites)-1) ],lowQ = 0.2,upQ = 0.8)
      n.brk[i] <- length(sites)
      start[i] <- min(sites)
      end[i] <- max(sites)
    }
    conf[which(dist.iqm < dist.iqm.cut )] <-"lc"
    chrom <- regions$chrom
    nbins <- regions$nseg
    restab[[cl]] <- data.frame(chrom,start,end,nbins,dist.iqm,n.brk,conf)
  }
  if(verbose) close(pb)
  
  bins <- remove.factors(data.frame(do.call(rbind,strsplit(colnames(highDensityRegions)," "))))
  bins[,2] <- as.numeric( bins[,2])
  bins[,3] <- as.numeric( bins[,3])
  rownames(bins) <- colnames(highDensityRegions)
  binsGR <- with(bins, GRanges(X1, IRanges(start=X2, end=X3)))
  highDensityRegionsHC <- highDensityRegions
  for(cl in names(restab)){
    lc <- restab[[cl]][which(restab[[cl]]$conf == "lc"),]
    if(nrow(lc) > 0){
      lcGR<- with(lc, GRanges(chrom, IRanges(start=start, end=end)))
      hits = GenomicAlignments::findOverlaps(binsGR,lcGR)
      highDensityRegionsHC[cl,rownames(bins[unique(queryHits(hits)),])] <- 0
    }
  }
  
  results <- chromo.regs(
    regions.summary = restab,
    high.density.regions = highDensityRegions,
    high.density.regions.hc = highDensityRegionsHC,
    seg.brk.dens = seg.brk.dens,
    sv.brk.dens = matrix(),
    seg.brk.common.dens = matrix(),
    sv.brk.common.dens = matrix(),
    segbrk = breaks,
    svbrk = list(),
    common.brk = list(),
    segdat = segdat,
    svdat = data.frame(),
    param=list(
        fc.pct = fc.pct,
        min.seg.size = min.seg.size,
        min.num.probes=min.num.probes, 
        low.cov = low.cov,
        clean.brk=clean.brk,
        window.size = window.size,
        slide.size = slide.size,
        num.breaks = num.breaks,
        num.sd = num.sd,
        dist.iqm.cut = dist.iqm.cut)
  )
return(results)
}

