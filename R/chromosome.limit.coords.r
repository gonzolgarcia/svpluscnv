#' 
#'
#' This function mapps defines mapping limmit regions based on a segmentation file 
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation, mapping
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' chr.lim <- chromosome.limit.coords(cnv)

  chromosome.limit.coords <- function(cnv){
    
    chrlist <- unique(cnv$chrom)
    chrmin <- chrmax <- list()
    for(chr in paste("chr",c(1:22,"X","Y"),sep="") ){
      if(chr %in% cnv$chrom){
        chrmin[[chr]] <- min(cnv[which(cnv$chrom == chr),"start"]) 
        chrmax[[chr]] <- max(cnv[which(cnv$chrom == chr),"end"]) 
        }
      }
    begin <- unlist(chrmin)
    end <- unlist(chrmax)
    return(data.frame(begin,end))
  }
  