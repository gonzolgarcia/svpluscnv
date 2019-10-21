#' 
#'
#' Caller for shattered genomic regions based on breakpoint densities
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param min.seg.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @param num.breaks (numeric) size in megabases of the genmome bin to compute break density 
#' @param num.sd (numeric) size in megabases of the sliding genmome window
#' @keywords CNV, segmentation
#' @export
#' @examples
#' shattered.regions.cnv()



shattered.regions.cnv <- function(seg,
                              fc.pct = 0.2, min.seg.size = 0, min.num.probes=0, 
                              low.cov = NULL, clean.brk=clean.brk,
                              window.size = 10,slide.size = 2,
                              num.breaks = 10, num.sd = 5,
                              verbose=FALSE){
  
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
  
  seg.brk.dens <- break.density(breaks, 
                                chr.lim = chr.lim, 
                                window.size = window.size, 
                                slide.size = slide.size,
                                verbose = verbose)
  
  
  # calculate inter quantile mean and standard deviation per sample
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

  
  restab <- restab_bychr <- list()
  for(cl in names(res)){
    if(length(res[[cl]]) > 0){
      if(verbose == TRUE) message(cl)
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
      restab_bychr[[cl]]<-aggregate(nseg~chrom,tabmerged,sum) 
    }
  }
  
  if(verbose == TRUE) message("Calculating breakpoint dispersion")
  for(cl in names(restab)){
    regions <-   restab[[cl]]
    br1 <- breaks[which(breaks$sample == cl),2:3]
    colnames(br1) <- c("chrom","pos")
    br1.gr <- with(br1, GRanges(chrom, IRanges(start=pos, end=pos)))
    regions_gr <- with(regions, GRanges(chrom, IRanges(start=start, end=end)))
    hits_1 = GenomicAlignments::findOverlaps(regions_gr,br1.gr)
    density <- dist.iqm <- rep(0,nrow(regions))
    for(i in 1:nrow(regions)){
      sites <- sort(unique(br1[subjectHits(hits_1)[which(queryHits(hits_1) == i)],"pos"]))
      density[i] <- median(abs(sites - mean(sites)))/(regions[i,"end"]-regions[i,"start"])
      dist.iqm[i]  <- IQM(sites[2:length(sites)] - sites[1:(length(sites)-1) ],lowQ = 0.2,upQ = 0.8)
      n.brk <- length(sites)
    }
    restab[[cl]] <- data.frame(regions,density,dist.iqm,n.brk)
  }
  
  return(list(
    summary_chromo_bydensity = restab,
    summary_chromo_bydensity_bychr = restab_bychr,
    high.density.regions = highDensityRegions,
    seg.brk.dens=seg.brk.dens,
    segbrk = breaks,
    segdat = segdat
  ))
}

