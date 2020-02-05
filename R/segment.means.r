#' 
#' Obtain the average weighted segment mean log2 ratios from each sample within a CNV segmentaton data.frame
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' cnv <- validate.cnv(segdat_lung_ccle)
#' ave.segmean(cnv)


####################


ave.segmean <- function(cnv){
  
  require(taRifx)  # contains remove.factors
  
    cnvdat <- validate.cnv(cnv)
  
  width <- cnvdat$end - cnvdat$start
  sample <- cnvdat[,c("sample")]
  segmean <- cnvdat[,c("segmean")]

  df <- aggregate(width~sample,data.frame(sample,width),sum)
  glen <- df$width
  names(glen) <- df$sample
  
  w.segmean <- segmean*width/glen[sample]
  df2 <- aggregate(w.segmean~sample,data.frame(sample,w.segmean),sum)
  ave <- df2$w.segmean
  names(ave) <- df2$sample
  return(ave)
  
  }


#' Obtain the median segment mean from a segmentaton file
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
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
    
    require(taRifx)  # contains remove.factors
    
    cnvdat <- validate.cnv(cnv)
    
    glen <- cnvdat$end-cnvdat$start
    sample <- cnvdat$sample
    segmean <- cnvdat$segmean
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

