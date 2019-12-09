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


break.density <- function(breaks, 
                         chr.lim, 
                         window.size = 10, 
                         slide.size=2,
                         chrlist=NULL,
                         verbose=TRUE){
  
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  

  stopifnot(window.size%%1 == 0, local = TRUE) 
  stopifnot(slide.size%%1 == 0, local = TRUE) 
  stopifnot((window.size/slide.size)%%1 == 0, local = TRUE) 

  stopifnot(ncol(breaks) >=3 , local = TRUE) 
  
  # make sure breaks format is correct
  breaks <- remove.factors(breaks[,1:3])
  colnames(breaks) <- c("sample","chr","pos")
  stopifnot(is.character(breaks$sample))
  stopifnot(is.character(breaks$chr))
  stopifnot(is.numeric(breaks$pos))
  
  # make sure chr.lim format is correct
  chr.lim<-chr.lim[,1:2]
  colnames(chr.lim) <- c("begin","end")
  stopifnot(is.numeric(chr.lim$begin))
  stopifnot(is.numeric(chr.lim$end))
  
  # make sure both chr.lim and breaks have same chromosome names 
  seqnames <- intersect(rownames(chr.lim),breaks$chr)
  stopifnot(length(seqnames) > 0) 
  
  # a template vector to save breakpoint counts 
  templatevector <- rep(0,length( unique(breaks$sample)))
  names(templatevector) <- unique(breaks$sample)
  
  WS <- window.size * 1e+6
  SS <- slide.size * 1e+6
  offset <- window.size/slide.size
  
  if(is.null(chrlist)) chrlist <- paste("chr",c(1:22,"X"),sep="")
  
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
      res_bp <-table(chr_breaks[intersect(which(chr_breaks$pos > start),which(chr_breaks$pos < stop)),"sample"])
      fragment[[paste(chr,start,stop)]][names(res_bp)] <- res_bp
    }
  }
  if(verbose) cat(" Done!\n")

  return( do.call(cbind,fragment))
  
}
