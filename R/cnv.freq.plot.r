#' 
#'
#' Plot CNV frequency across the human genome from a seg file sontaining multiple samples 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param ch.pct (numeric) percentage CNV gain/loss for a segment to be considered changed (i.e. 0.2 = 20 percent change 0.8 < segmean && segmean > 1.2)
#' @param g.bin (numeric) size in megabases of the genmome bin to compute break density 
#' @param hg (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param cex.axis,cex.lab,axis.line,label.line plot parameters
#' @keywords CNV, segmentation
#' @export
#' @examples
#' cnv.freq.plot()

cnv.freq.plot <- function(seg=NULL,
                          ch.pct=0.2,
                          g.bin=1,
                          hg="hg19",
                          cex.axis=1,cex.lab=1,axis.line=-1.5,label.line=-1.2){
  
  require(D3GB,quietly = TRUE,warn.conflicts = FALSE)
  require(taRifx,quietly = TRUE,warn.conflicts = FALSE)
  require(tidyr,quietly = TRUE,warn.conflicts = FALSE)
  require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)
  
  segdat <- validate.seg(seg)
  
  if(hg == "hg19"){ bands <- GRCh37.bands
  }else if(hg=="h38"){ bands <- GRCh38.bands}
  
  centromeres <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
  names(centromeres) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")
  
  # define chromosome mapped limits and the global genome coordinates for each chromosome start
  chrlimits <-   chromosome.limit.coords(segdat)
  offset <- c(0,sapply(1:(nrow(chrlimits)-1), function(i) sum(chrlimits[1:i,"end"]) + length(1:i)*g.bin ))
  chrlabelpos <- offset + chrlimits$end/2
  chrlimits <-data.frame(offset,chrlimits,chrlabelpos)
  
  g.bin.mb <- g.bin*1e6
  
  message("Generating binned genome map ")
  chrbins <- list()
  for(chr in rownames(chrlimits)){
    seqpos <- seq(chrlimits[chr,"begin"],chrlimits[chr,"end"]+g.bin.mb,g.bin.mb)
    ranges <-  t(sapply(2:length(seqpos), function(i) c(seqpos[i-1],seqpos[i])))
    chrcol<- rep(chr,length(seqpos)-1)
    segcol_del <- segcol_gain <- rep("grey",length(chrcol))
    segcol_del[which(ranges[,2] <= centromeres[chr])] <- "lightblue"
    segcol_del[which(ranges[,2] > centromeres[chr])] <- "blue"
    segcol_gain[which(ranges[,2] <= centromeres[chr])] <- "salmon"
    segcol_gain[which(ranges[,2] > centromeres[chr])] <- "red"
    chrbins[[chr]] <- data.frame(chrcol,ranges,segcol_del,segcol_gain)
  }
  
  chrbins.df <- data.frame(do.call(rbind,unname(chrbins) ))
  colnames(chrbins.df) <- c("chr","start","end","segcol_del","segcol_gain")
  rownames(chrbins.df) <-  unite(chrbins.df[,1:3],paste)$paste
  chrbins.df<-remove.factors(chrbins.df)
  
  message("Calculating mean segmean per genomic bin")
  # find overlaps between bins and cnv segments
  binsGR <- with(chrbins.df, GRanges(chr, IRanges(start=start, end=end)))
  segGR <- with(segdat, GRanges(chrom, IRanges(start=start, end=end)))
  hits <-GenomicAlignments::findOverlaps(binsGR,segGR)
  
  outmat <- matrix(ncol=length(unique(segdat$sample)),nrow=nrow(chrbins.df))
  colnames(outmat) <- unique(segdat$sample)
  rownames(outmat) <- rownames(chrbins.df)

  for(i in 1:nrow(chrbins.df)){ 
    segtmp<- segdat[subjectHits(hits)[which(queryHits(hits) == i)],]
    a <- aggregate(segmean~sample,segtmp, sum)  
    outmat[i,a$sample]<- a$segmean
  }

  message("Calculating gain/loss frequencies per genomic bin")
  outmat[which(is.na(outmat),arr.ind=T)] <- 0
  
  outmat_gain<-outmat_loss<-outmat
  outmat_gain[]<-outmat_loss[]<-0
  
  outmat_gain[which(outmat > log2(1+ch.pct), arr.ind=T)] <-  1
  outmat_loss[which(outmat < log2(1-ch.pct), arr.ind=T)] <-  1
  allgains <- apply(outmat_gain,1,sum)/ncol(outmat_gain)
  allloss <- apply(outmat_loss,1,sum)/ncol(outmat_loss)
  
  plot.end<- chrlimits$offset[nrow(chrlimits)]+chrlimits$end[nrow(chrlimits)]
  bin.loc <- chrlimits[chrbins.df[names(allgains),"chr"],"offset"] + chrbins.df[names(allgains),"start"]
  
  message("Plotting ...")
  altcols <- rep(c(rgb(0.1,0.1,0.1,alpha=0.1),rgb(0.8,0.8,0.8,alpha=0.1)),12)
  altcols2<- rep(c(rgb(0.1,0.1,0.1,alpha=1),rgb(0.4,0.4,0.4,alpha=1)),12)
  
  plot(x=NULL,y=NULL,xlim=c(0,plot.end),ylim=c(-1,1),bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
  for(i in 1:length(chrlimits$offset) ) rect( chrlimits$offset[i],-1,chrlimits$offset[i]+chrlimits$end[i],1, col=altcols[i],border=NA )
  points(bin.loc,allgains,type='h',col=chrbins.df$segcol_gain)
  points(bin.loc,-allloss,type='h',col=chrbins.df$segcol_del)
  lines(c(0,plot.end),c(0,0),col="lightgrey")
  lines(c(0,plot.end),c(0.5,0.5),col="lightgrey",lty=3)
  lines(c(0,plot.end),c(-0.5,-0.5),col="lightgrey",lty=3)
  mtext(gsub("chr","",rownames(chrlimits))[seq(1,nrow(chrlimits),2)],side=1,at=chrlimits$chrlabelpos[seq(1,nrow(chrlimits),2)],las=1,col=altcols2[seq(1,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
  mtext(gsub("chr","",rownames(chrlimits))[seq(2,nrow(chrlimits),2)],side=3,at=chrlimits$chrlabelpos[seq(2,nrow(chrlimits),2)],las=1,col=altcols2[seq(2,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
  axis(2,c(100,50,0,50,100),at=c(-1,-0.5,0,0.5,1),las=1,line= axis.line,cex.axis=cex.axis)
}

