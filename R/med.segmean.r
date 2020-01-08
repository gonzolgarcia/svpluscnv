#' 
#'
#' Obtain the median segment mean from a segmentaton file
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' med.segmean(seg)


####################


med.segmean <- function(seg){

  require(taRifx)  # contains remove.factors
  
  segdat <- validate.seg(seg)
  
  glen <- segdat$end-segdat$start
  sample <- segdat$sample
  segmean <- segdat$segmean
  df <- remove.factors(data.frame(sample,glen,segmean))
  out <-rep(NA,length(unique(df$sample)))
  names(out) <- unique(df$sample)
  
  pb <- txtProgressBar(style=3)
  cc <-0
  tot <- length(unique(sample))
  
  
  for(i in unique(df$sample)){
    cc <- cc+1
    
    minidf <- df[which(df$sample == i),]
    miniord <-minidf[order(minidf$segmean),]
    medseg <- which(abs(cumsum(miniord$glen)/sum(miniord$glen) - 0.5) == min(abs(cumsum(miniord$glen)/sum(miniord$glen) - 0.5)))
    out[i] <- mean(miniord$segmean[medseg])
    
    setTxtProgressBar(pb, cc/tot)
    
  }
  close(pb)
  return(out)
}

