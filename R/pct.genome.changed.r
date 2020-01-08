#' 
#'
#' Calculate the percentage of genome changed by CNV
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) percentage CNV gain/loss for a segment to be considered changed (i.e. 0.2 = 20 percent change 0.8 < segmean && segmean > 1.2)
#' @param discard.sex (logical) whether sex chromosomes should be included or not
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' pct.genome.changed(seg)

pct.genome.changed <- function(seg, 
                               fc.pct=0.2, 
                               discard.sex=TRUE){

  require(taRifx)  # contains remove.factors

  segdat <- validate.seg(seg)
  
  if(discard.sex == TRUE) segdat <- segdat[which(!segdat$chrom %in% c("chrX","chrY")),]

  width <- segdat$end - segdat$start
  segmean <- segdat$segmean
  sample <- segdat$sample
  df <- remove.factors(data.frame(sample,width,segmean))
  idx_changed <- c(which(df$segmean < log2(1-fc.pct)),which(df$segmean >= log2(1+fc.pct)))
  idx_normal <- setdiff(1:nrow(df),idx_changed)
  df_normal <- df[idx_normal,]
  df_changed <-  df[idx_changed,]
  
  length_changed_df <- aggregate(width~sample ,df_changed,sum)
  length_normal_df <- aggregate(width~sample ,df_normal,sum)
  
  nochange <- setdiff(length_normal_df$sample,length_changed_df$sample)
  fullchange <- setdiff(length_changed_df$sample,length_normal_df$sample)
  nochange_x <- rep(0,length(nochange))
  names(nochange_x) <- nochange
  fullchange_x <- rep(0,length(fullchange))
  names(fullchange_x) <- fullchange
  
  length_changed <- c(length_changed_df[,2],nochange_x)
  names(length_changed)<- c(length_changed_df[,1],nochange)
  
  length_normal <- c(length_normal_df[,2],fullchange_x)
  names(length_normal)<- c(length_normal_df[,1],fullchange)
  
  pct.change<- length_changed/apply(cbind(length_normal[names(length_changed)],length_changed),1,sum)
  
  return(pct.change)
}


