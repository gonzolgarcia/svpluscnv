#' 
#'
#' Obtain the average weighted segment mean from each sample in a segmentaton file
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation
#' @export
#' @examples
#' ave.segmean()


####################


ave.segmean <- function(seg){
  
  require(taRifx)  # contains remove.factors
  
  segdat <- validate.seg(seg)
  
  width <- segdat$end - segdat$start
  sample <- segdat[,c("sample")]
  segmean <- segdat[,c("segmean")]

  df <- aggregate(width~sample,data.frame(sample,width),sum)
  glen <- df$width
  names(glen) <- df$sample
  
  w.segmean <- segmean*width/glen[sample]
  df2 <- aggregate(w.segmean~sample,data.frame(sample,w.segmean),sum)
  ave <- df2$w.segmean
  names(ave) <- df2$sample
  return(ave)
  
  }

