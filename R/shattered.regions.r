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


shattered.regions <- function(seg, sv,
                              fc.pct = 0.2, min.seg.size = 0, min.num.probes=0, low.cov = NULL,
                              window.size = 10,slide.size = 2,
                              num.seg.breaks = 6, num.seg.sd = 5,
                              num.sv.breaks = 6, num.sv.sd = 5,
                              num.common.breaks = 3, num.common.sd = 3,
                              maxgap=10000, chrlist=NULL,
                              verbose=TRUE){
  

  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  
  segdat <- validate.seg(seg)
  svdat <- validate.sv(sv)
  
  chr.lim <- chromosome.limit.coords(segdat)
  
  segbrk <- seg.breaks(seg = segdat, 
                       fc.pct = fc.pct, 
                       min.seg.size = min.seg.size, 
                       low.cov = low.cov, 
                       verbose = verbose)
  svbrk <- sv.breaks(sv = svdat, 
                      low.cov = low.cov)
  
  common.brk <- match.variant.breaks(segdat,
                                     svdat,
                                     maxgap = maxgap, 
                                     low.cov=low.cov,
                                     verbose = verbose)
  segbrk.common <-  common.brk$seg_validated
  svbrk.common <-  sv.breaks(common.brk$sv_validated, 
                              low.cov = low.cov)
  
  segbrk.dens <- break.density(segbrk,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  svbrk.dens <- break.density(svbrk,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  segbrk.common.dens <- break.density(segbrk.common,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  svbrk.common.dens <- break.density(svbrk.common,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  
  commonSamples <- intersect(segbrk$sample,svbrk$sample)
  if(length(commonSamples) == 0) stop("There is no common samples between seg and sv input datasets.") 
  
  # calculate inter quantile mean and standard deviation per sample
  iqmdata1 <- apply(segbrk.dens,1,IQM,lowQ=0.2,upQ=0.8)
  sddata1  <- apply(segbrk.dens,1,IQSD,lowQ=0.2,upQ=0.8)
  iqmdata2 <- apply(svbrk.dens,1,IQM,lowQ=0.2,upQ=0.8)
  sddata2  <- apply(svbrk.dens,1,IQSD,lowQ=0.2,upQ=0.8)
  iqmdata3 <- apply(segbrk.common.dens,1,IQM,lowQ=0.2,upQ=0.8)
  sddata3  <- apply(segbrk.common.dens,1,IQSD,lowQ=0.2,upQ=0.8)
  
  a <- sapply(commonSamples, function(i) names(which(segbrk.dens[i,] > iqmdata1[i]+num.seg.sd*sddata1[i] )))
  b <- sapply(commonSamples, function(i) names(which(segbrk.dens[i,] >= num.seg.breaks)))
  c <- sapply(commonSamples, function(i) names(which(svbrk.dens[i,] > iqmdata2[i]+num.sv.sd*sddata2[i] )))
  d <- sapply(commonSamples, function(i) names(which(svbrk.dens[i,] >= num.sv.breaks)))
  e <- sapply(commonSamples, function(i) names(which(segbrk.common.dens[i,] > iqmdata3[i]+num.common.sd*sddata3[i] )))
  f <- sapply(commonSamples, function(i) names(which(segbrk.common.dens[i,] >= num.common.breaks)))
  
  
  # condition for chromothripsis: at least n=breaks > 6 (sv SND seg)  AND n-breaks > u+2*sd (sv AND seg) 
  res <- sapply(commonSamples,function(i) Reduce(intersect, list(a[[i]],b[[i]],c[[i]],d[[i]],e[[i]],f[[i]])))

  highDensityRegions <- segbrk.dens[commonSamples,]
  highDensityRegions[] <- 0
  for(cl in commonSamples) highDensityRegions[cl,res[[cl]]] <- 1
  
  #res <- res[which(unlist(lapply(res,length)) >0)]
  
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
        chrom <- as.character(unique(tab[as.numeric(agglist[[i]]),"chrom"]))
        start<-min( tab[as.numeric(agglist[[i]]),"start"])
        end<-max( tab[as.numeric(agglist[[i]]),"end"])
        segNum <- length(agglist[[i]])
        agglistUniq[[i]] <-  c(chrom,start,end,segNum)
      }
      tabmerged <- remove.factors(as.data.frame(do.call(rbind,agglistUniq)))
      for(i in 2:4) tabmerged[,i] <- as.numeric( tabmerged[,i] )
      colnames(tabmerged) <- c("chrom","start","end","nseg")
      restab[[cl]] <- tabmerged
      restab_bychr[[cl]]<-aggregate(nseg~chrom,tabmerged,sum) 
    }
  }
  out <- list(
    regions.summary = restab,
    regions.summary.bychr = restab_bychr,
    seg.brk.dens = segbrk.dens,
    sv.brk.dens = svbrk.dens,
    seg.brk.common.dens = segbrk.common.dens,
    sv.brk.common.dens = svbrk.common.dens,
    segbrk = segbrk,
    svbrk = svbrk,
    common.brk = common.brk,
    segdat = segdat,
    svdat = svdat
  )
  return(out)
}

