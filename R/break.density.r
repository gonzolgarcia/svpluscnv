#' 
#'
#' Identify CNV breakpoints provided a segmentation file
#' @param breaks (data.frame) 3 column table (sample, chrom, position) i.e returned by read.depth.breaks
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) returned by chromosome.limit.coords
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @keywords CNV, segmentation
#' @export
#' @examples
#' break.density()


break.density <- function(brk, 
                         chr.lim, 
                         window.size = 10, 
                         slide.size=2,
                         chrlist=NULL,
                         verbose=TRUE){
  
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  
  # make sure breaks format is correct
  breaks <-brk$breaks
  
  # make sure chr.lim format is correct
  chr.lim<- chr.lim[,1:2]
  colnames(chr.lim) <- c("begin","end")
  
  # make sure both chr.lim and breaks have same chromosome names 
  seqnames <- intersect(rownames(chr.lim),breaks$chr)
  stopifnot(length(seqnames) > 0) 
  
  # a template vector to save breakpoint counts 
  templatevector <- brk$brk.burden
  templatevector[]<-0
  
  WS <- window.size * 1e+6
  SS <- slide.size * 1e+6
  offset <- window.size/slide.size
  
  if(is.null(chrlist)) chrlist <- rownames(chr.lim)
  
  # count breaks for each chromosome for each fragment
  fragment <- list()
  for(chr in  chrlist){

    if(verbose) cat("\r",chr)

    chr_breaks <- breaks[which(breaks$chr == chr),]
    frag <- seq(chr.lim[chr,"begin"],chr.lim[chr,"end"]+SS,SS)
    for(i in (1+offset):length(frag)){
      start <- frag[i - offset]
      stop <- frag[i]
      fragment[[paste(chr,start,stop)]] <- templatevector
      break.position <- chr_breaks$start + (chr_breaks$end - chr_breaks$start)/2
      res_bp <- table(chr_breaks[intersect(which(break.position > start),which(break.position < stop)),"sample"])
      fragment[[paste(chr,start,stop)]][names(res_bp)] <- res_bp
    }
  }
  if(verbose) cat(" Done!\n")

  return( do.call(cbind,fragment))
  
}
