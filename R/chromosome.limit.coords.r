#' 
#'
#' This function mapps defines mapping limmit regions based on a segmentation file 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation, mapping
#' @export
#' @examples
#' pcent.genome.changed()

  chromosome.limit.coords <- function(seg){
    
    require(taRifx)  # contains remove.factors
    
    segdat <- validate.seg(seg)
    
    chrlist <- unique(segdat$chrom)
    chrmin <- chrmax <- list()
    for(chr in paste("chr",c(1:22,"X","Y"),sep="") ){
      if(chr %in% segdat$chrom){
        chrmin[[chr]] <- min(segdat[which(segdat$chrom == chr),"start"]) 
        chrmax[[chr]] <- max(segdat[which(segdat$chrom == chr),"end"]) 
        }
      }
    begin <- unlist(chrmin)
    end <- unlist(chrmax)
    return(data.frame(begin,end))
  }
  