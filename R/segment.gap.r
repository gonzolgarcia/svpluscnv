#' 
#' 
#' Filling gaps in a segmentation dat abased on chromosome limits defined withing the complete segmentation dataset
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @keywords CNV, segmentation
#' @export
#' @examples
#' segment.gap()

segment.gap <- function(seg,chrlist=NULL,verbose=TRUE){
  segdat <- validate.seg(seg)
  if(is.null(chrlist)) chrlist <- paste("chr",c(1:22,"X"), sep="" )
  
  chrlims <- chromosome.limit.coords(segdat)
  
  segdat <- segdat[order(segdat$start),]
  segdat <- segdat[order(match(segdat$chrom, chrlist)),]
  segdat <- segdat[order(segdat$sample),]
  
  
  if(verbose){
    message("Filling gaps is the segmentation data.frame")
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- nrow(segdat)
  }
  newsegments<-list()
  if(segdat[1,"start"] > chrlims[segdat[1,"chrom"],"begin"]){ 
    newsegments[["1"]] <- data.frame(segdat[1,c("sample","chrom")],chrlims[segdat[1,"chrom"],"begin"],segdat[1,"start"]-1,0,segdat[1,"segmean"])
  }
  
  for(i in 2:nrow(segdat)){
    if(segdat[i,"chrom"] == segdat[i-1,"chrom"] ){
      if( segdat[i,"start"] - segdat[i-1,"end"] > 100){
        newsegments[[as.character(i)]] <- data.frame(segdat[i,c("sample","chrom")],segdat[i-1,"end"]+1,segdat[i,"start"]-1,0,mean(segdat[c(i,i-1),"segmean"]) )
      }
    }else{
      if(segdat[i,"start"] > chrlims[segdat[i,"chrom"],"begin"]){ 
        newsegments[[as.character(i)]] <- data.frame(segdat[i,c("sample","chrom")],chrlims[segdat[i,"chrom"],"begin"],segdat[i,"start"]-1,0,segdat[i,"segmean"])
      }
      if(segdat[i-1,"end"] < chrlims[segdat[i-1,"chrom"],"end"]){
        newsegments[[as.character(i)]] <- data.frame(segdat[i-1,c("sample","chrom")],segdat[i-1,"end"]+1,chrlims[segdat[i-1,"chrom"],"end"],0,segdat[i,"segmean"])
      }
    }
    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)
  }
  if(segdat[i,"end"] < chrlims[segdat[i,"chrom"],"end"]) newsegments[[as.character(i)]] <- data.frame(segdat[i,c("sample","chrom")],segdat[i,"end"]+1,chrlims[segdat[i,"chrom"],"end"],0,segdat[i-1,"segmean"])
  if(verbose) close(pb)
  
  newsegments<- lapply(newsegments, setNames, colnames(segdat))
  
  segout <- rbind(segdat, do.call(rbind,newsegments))
  segout <- segout[order(segout$start),]
  segout <- segout[order(segout$sample),]
  
  return(segout)
}
