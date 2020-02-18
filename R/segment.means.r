#' Obtain the average weighted segment mean log2 ratios from each sample within a CNV segmentaton data.frame
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @return (numeric) a vector containing the weighted average logR from segmented data
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' cnv <- validate.cnv(segdat_lung_ccle)
#' ave.segmean(cnv)


####################


ave.segmean <- function(cnv){
  
stopifnot(cnv@type == "cnv")
cnvdat <- cnv@data
    

  width <- as.numeric(cnvdat$end - cnvdat$start)
  sample <- cnvdat$sample
  segmean <- cnvdat$segmean

  df <- stats::aggregate(width~sample,data.table(sample,width),sum)
  glen <- df$width
  names(glen) <- df$sample
  
  w.segmean <- segmean*width/glen[sample]
  df2 <- stats::aggregate(w.segmean~sample,data.table(sample,w.segmean),sum)
  ave <- df2$w.segmean
  names(ave) <- df2$sample
  return(ave)
  
  }


#' Obtain the median segment mean from a segmentaton file; The median refers to the logR that occupies a center of all segments ordered by their log ratio
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @return (numeric) a vector containing the median logR value of a segmented data.frame
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' med.segmean(cnv)


####################


med.segmean <- function(cnv){
    
    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    glen <- as.numeric(cnvdat$end-cnvdat$start)
    sample <- cnvdat$sample
    segmean <- cnvdat$segmean
    dt <- data.table(sample,glen,segmean)
    out <-rep(NA,length(unique(dt$sample)))
    names(out) <- unique(dt$sample)
    
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(unique(sample))
    
    
    for(i in unique(dt$sample)){
        cc <- cc+1
        
        minidf <- dt[which(dt$sample == i)]
        miniord <-minidf[order(minidf$segmean)]
        medseg <- which(abs(cumsum(miniord$glen)/sum(miniord$glen) - 0.5) == min(abs(cumsum(miniord$glen)/sum(miniord$glen) - 0.5)))
        out[i] <- mean(miniord$segmean[medseg])
        
        setTxtProgressBar(pb, cc/tot)
        
    }
    close(pb)
    return(out)
}

