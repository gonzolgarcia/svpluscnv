#' 
#'
#' Plot CNV frequency across the human genome from a seg file sontaining multiple samples 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) percentage CNV gain/loss for a segment to be considered changed (i.e. 0.2 = 20 percent change 0.8 < segmean && segmean > 1.2)
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param g.bin (numeric) size in megabases of the genmome bin to compute break density 
#' @param sampleids (character) vector containing list of samples to include in plot. if set to NULL, all samples in the input will be used
#' @param cex.axis,cex.lab,label.line (numeric) plot parameters
#' @param plot (logical) whether to produce a graphical output
#' @param summary (logical)  whether to return an object with a the summary
#' @param verbose (logical) 
#' @keywords CNV, segmentation, plot
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' seg <- validate.seg(nbl_segdat)
#' 
#' cnv.freq(seg, genome.v = "hg19")

cnv.freq <- function(seg,
                     fc.pct= 0.2,
                     genome.v= "hg19",
                     g.bin= 1,
                     sampleids=NULL,
                     cex.axis= 1,
                     cex.lab= 1,
                     label.line= -1.2,
                     plot=TRUE,
                     summary=TRUE,
                     verbose=TRUE){
  
require(D3GB,quietly = TRUE,warn.conflicts = FALSE)
require(taRifx,quietly = TRUE,warn.conflicts = FALSE)
require(tidyr,quietly = TRUE,warn.conflicts = FALSE)
require(GenomicRanges,quietly = TRUE,warn.conflicts = FALSE)

segdat <- validate.seg(seg)
if(!is.null(sampleids)) segdat <- segdat[which(segdat$sample %in% sampleids),]
  
if(genome.v == "hg19"){ bands <- GRCh37.bands
}else if(genome.v=="h38"){ bands <- GRCh38.bands}

centromeres <- bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"start"]
names(centromeres) <- paste("chr",bands[intersect(which(bands$score == "acen"),grep("q",bands$name)),"chr"],sep="")

# define chromosome mapped limits and the global genome coordinates for each chromosome start
chrlimits <-   chromosome.limit.coords(segdat)
offset <- c(0,sapply(1:(nrow(chrlimits)-1), function(i) sum(chrlimits[1:i,"end"]) + length(1:i)*g.bin ))
chrlabelpos <- offset + chrlimits$end/2
chrlimits <- data.frame(offset,chrlimits,chrlabelpos)
  
g.bin.mb <- g.bin*1e6
  
if(verbose) message("Generating binned genome map ")

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
  
  if(verbose) message("Calculating mean segmean per genomic bin")
  # find overlaps between bins and cnv segments
  binsGR <- with(chrbins.df, GRanges(chr, IRanges(start=start, end=end)))
  segGR <- with(segdat, GRanges(chrom, IRanges(start=start, end=end)))
  hits <-GenomicAlignments::findOverlaps(binsGR,segGR)
  
  outmat <- matrix(ncol=length(unique(segdat$sample)),nrow=nrow(chrbins.df))
  colnames(outmat) <- unique(segdat$sample)
  rownames(outmat) <- rownames(chrbins.df)

  for(i in 1:nrow(chrbins.df)){
    segtmp<- segdat[subjectHits(hits)[which(queryHits(hits) == i)],]
    if(nrow(segtmp)>0){
      a <- aggregate(segmean~sample,segtmp, sum)  
      outmat[i,a$sample]<- a$segmean
    }else{
      outmat[i,a$sample]<- NA
    }
  }

  if(verbose) message("Calculating gain/loss frequencies per genomic bin")
  outmat[which(is.na(outmat),arr.ind=T)] <- 0
  
  outmat_gain<-outmat_loss<-outmat
  outmat_gain[]<-outmat_loss[]<-0
  nsamples <- ncol(outmat_gain)
  
  outmat_gain[which(outmat > log2(1+fc.pct), arr.ind=T)] <-  1
  outmat_loss[which(outmat < log2(1-fc.pct), arr.ind=T)] <-  1
  freq.gains <- apply(outmat_gain,1,sum)/nsamples
  freq.loss <- apply(outmat_loss,1,sum)/nsamples
  
if(plot == TRUE){
    plot.end<- chrlimits$offset[nrow(chrlimits)]+chrlimits$end[nrow(chrlimits)]
    bin.loc <- chrlimits[chrbins.df[names(freq.gains),"chr"],"offset"] + chrbins.df[names(freq.gains),"start"]

    if(verbose) message("Plotting ...")
    altcols <- rep(c(rgb(0.1,0.1,0.1,alpha=0.1),rgb(0.8,0.8,0.8,alpha=0.1)),12)
    altcols2<- rep(c(rgb(0.1,0.1,0.1,alpha=1),rgb(0.4,0.4,0.4,alpha=1)),12)
  
    plot(x=NULL,y=NULL,xlim=c(0,plot.end),ylim=c(-1,1),bty='n',xaxt='n',yaxt='n',xlab="",ylab="")
    for(i in 1:length(chrlimits$offset) ) rect( chrlimits$offset[i],-1,chrlimits$offset[i]+chrlimits$end[i],1, col=altcols[i],border=NA )
    points(bin.loc,freq.gains,type='h',col=chrbins.df$segcol_gain)
    points(bin.loc,-freq.loss,type='h',col=chrbins.df$segcol_del)
    lines(c(0,plot.end),c(0,0),col="lightgrey")
    lines(c(0,plot.end),c(0.5,0.5),col="lightgrey",lty=3)
    lines(c(0,plot.end),c(-0.5,-0.5),col="lightgrey",lty=3)
    mtext(gsub("chr","",rownames(chrlimits))[seq(1,nrow(chrlimits),2)],side=1,at=chrlimits$chrlabelpos[seq(1,nrow(chrlimits),2)],las=1,col=altcols2[seq(1,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
    mtext(gsub("chr","",rownames(chrlimits))[seq(2,nrow(chrlimits),2)],side=3,at=chrlimits$chrlabelpos[seq(2,nrow(chrlimits),2)],las=1,col=altcols2[seq(2,nrow(chrlimits),2)],line=label.line,cex=cex.lab)
    mtext("Frequency",side=4,line=1)
    mtext("#samples",side=2,line=1)
    axis(4,c(100,50,0,50,100),at=c(-1,-0.5,0,0.5,1),las=1,pos=plot.end, cex.axis=cex.axis)
    axis(2,c(nsamples,round(nsamples/2),0,round(nsamples/2),nsamples),at=c(-1,-0.5,0,0.5,1),las=1, pos=0, cex.axis=cex.axis)
}
  
if(summary == TRUE){
    summary <- data.frame(chrbins.df[,c("chr","start","end")],freq.gains,freq.loss)
    return(list(freqsum = summary,
                bin.mat = outmat,
                param = list(
                    fc.pct= fc.pct,
                    genome.v= genome.v,
                    g.bin= g.bin,
                    sampleids=sampleids,
                    cex.axis= cex.axis,
                    cex.lab= cex.lab,
                    label.line= label.line   
                )))
    }
}

