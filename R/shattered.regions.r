#' 
#'
#' Caller for shattered genomic regions based on breakpoint densities
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param min.seg.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param clean.brk (numeric) n of redundant breakpoints to filter out 
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @param num.seg.breaks (numeric) number of segmentation breakpoints per segments to be considered high-density break 
#' @param num.seg.sd (numeric) number of standard deviations above the sample average for num.seg.breaks
#' @param num.sv.breaks (numeric) number of SV breakpoints per segments to be considered high-density break  
#' @param num.sv.sd (numeric) number of standard deviations above the sample average for num.sv.breaks
#' @param num.common.breaks (numeric) number of common SV and segmentation breakpoints per segments to be considered high-density break
#' @param num.common.sd (numeric) number of standard deviations above the sample average for num.common.breaks
#' @param interleaved.cut (numeric) 0-1 value indicating percentage of interleaved SV breaks
#' @param dist.iqm.cut (numeric) interquantile average of the distance between breakpoints within a shattered region
#' @keywords chromothripsis, chromoplexy, chromosome shattering
#' @export
#' @examples
#' shattered.regions.cnv()


shattered.regions <- function(seg, sv,
                              fc.pct = 0.2, min.seg.size = 0, min.num.probes=0, 
                              low.cov = NULL,clean.brk=NULL,
                              window.size = 10,slide.size = 2,
                              num.seg.breaks = 6, num.seg.sd = 5,
                              num.sv.breaks = 6, num.sv.sd = 5,
                              num.common.breaks = 3, num.common.sd = 3,
                              maxgap=10000, chrlist=NULL, 
                              interleaved.cut = 0.5, dist.iqm.cut = 1e+05, 
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
                       clean.brk = clean.brk,
                       verbose = verbose)
  
  svbrk <- sv.breaks(sv = svdat, 
                      low.cov = low.cov)
  
  common.brk <- match.breaks(segbrk,svbrk,
                             maxgap = maxgap, 
                             verbose = verbose)
  
  if(verbose) message("Mapping CNV breakpoints across the genome:")
  segbrk.dens <- break.density(segbrk,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  if(verbose) message("Mapping SV breakpoints across the genome:")
  svbrk.dens <- break.density(svbrk,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  if(verbose) message("Mapping CNV validated breakpoints across the genome:")
  segbrk.common.dens <- break.density(common.brk$brk1_match,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  if(verbose) message("Mapping SV validated breakpoints across the genome:")
  svbrk.common.dens <- break.density(common.brk$brk2_match,chr.lim = chr.lim, window.size = window.size, slide.size = slide.size, chrlist = chrlist,verbose = verbose)
  
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
  if(verbose){
    message("Locating shattered regions...")
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(res)
  }
  
  restab  <- list()
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
    }
  }
  if(verbose) close(pb)
  
  results <- list(
    regions.summary = restab,
    high.density.regions = highDensityRegions,
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
  
  results_eval <- shattered.eval(results, interleaved.cut = interleaved.cut, dist.iqm.cut = dist.iqm.cut, verbose = verbose)
  
  bins <- remove.factors(data.frame(do.call(rbind,strsplit(colnames(highDensityRegions)," "))))
  bins[,2] <- as.numeric( bins[,2])
  bins[,3] <- as.numeric( bins[,3])
  rownames(bins) <- colnames(highDensityRegions)
  binsGR <- with(bins, GRanges(X1, IRanges(start=X2, end=X3)))
  highDensityRegionsHC <- highDensityRegions
  for(cl in names(results_eval$regions.summary)){
    lc <- results_eval$regions.summary[[cl]][which(results_eval$regions.summary[[cl]]$conf == "lc"),]
    if(nrow(lc) > 0){
      lcGR<- with(lc, GRanges(chrom, IRanges(start=start, end=end)))
      hits = GenomicAlignments::findOverlaps(binsGR,lcGR)
      highDensityRegionsHC[cl,rownames(bins[unique(queryHits(hits)),])] <- 0
    }
  }
  results_eval[["high.density.regions.hc"]] <- highDensityRegionsHC
  
  return(results_eval)
}

