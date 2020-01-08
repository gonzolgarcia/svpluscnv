#' 
#'
#' Obtain a matrix with the weighted average CN per chromosome arm 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param genome.v (`hg19` or `hg38`) reference genome version to draw chromosome limits and centromeres
#' @param verbose (logical)
#' @keywords CNV, segmentation, chromosome arm
#' @export
#' @examples
#' chr.arm.cnv()


chr.arm.cnv <- function(seg,
                    genome.v="hg19",
                    verbose=TRUE){

  require(D3GB,quietly = TRUE, warn.conflicts = FALSE)
  require(taRifx,quietly = TRUE, warn.conflicts = FALSE)
  require(tidyr,quietly = TRUE, warn.conflicts = FALSE)
  require(GenomicRanges,quietly = TRUE, warn.conflicts = FALSE)
  
  segdat <- validate.seg(seg)
  
  if(genome.v %in% c("GRCh37","hg19")){ 
    bands <- GRCh37.bands
  }else if(genome.v %in% c("GRCh38","hg38")){ 
    bands <- GRCh38.bands
  }else{stop("Genome version not provided")}
  
  centromeres_start <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
  centromeres_end <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"end"]
  names(centromeres_start) <-  names(centromeres_end) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")
  
  chrends <- chromosome.limit.coords(segdat)
  chrarms <- rbind(cbind(chrends$begin,centromeres_start[rownames(chrends)]),cbind(centromeres_end[rownames(chrends)],chrends$end))
  chrarms <- data.frame(rownames(chrarms),chrarms,c(paste(rownames(chrends),"p",sep=""), paste(rownames(chrends),"q",sep="")))
  colnames(chrarms) <- c("chr","start","end","arm")
  
  chrarms <- remove.factors(chrarms[which(chrarms$end -chrarms$start > 0),])
  rownames(chrarms) <- NULL
  
  chrarmsGR <- with(chrarms,GRanges(chr, IRanges(start=start, end=end)))

  segdat_gr <- with(segdat, GRanges(chrom, IRanges(start=start, end=end)))
  hits <-GenomicAlignments::findOverlaps(chrarmsGR,segdat_gr)
  
  armcnvmat <- matrix(ncol=length(unique(segdat$sample)), nrow=nrow(chrarms) )
  colnames(armcnvmat) <- unique(segdat$sample)
  rownames(armcnvmat) <- chrarms$arm
  
  for(i in unique(queryHits(hits))){ 
    arm <- chrarms[i,"arm"]
    
    if(verbose) message(arm)
    
    armdf <- segdat[subjectHits(hits)[which(queryHits(hits) == i)],]
    armdf[which(armdf[,"start"] < chrarms[i,"start"]),"start"] <- chrarms[i,"start"]
    armdf[which(armdf[,"end"] > chrarms[i,"end"]),"end"] <- chrarms[i,"end"]

    arm.width <- armdf[,"end"] - armdf[,"start"]
    armdf <- data.frame(armdf,arm.width)
    armlength <- aggregate(arm.width~sample,armdf,sum)[,2]
    names(armlength) <- aggregate(arm.width~sample,armdf,sum)[,1]
    part <- armdf$segmean * armdf$arm.width / armlength[armdf$sample]
    
    armdf <- data.frame(armdf,arm.width,part,armlength[armdf$sample])
    
    meanArmSegment <- aggregate(part~sample,armdf,sum)
    
    num <-  as.numeric(meanArmSegment[,2])
    names(num) <- as.character(meanArmSegment[,1])
    armcnvmat[arm,names(num)] <- num
  }
  return(armcnvmat)
}



