#' Circos plot combining segmentation and SV calls
#' @param cnv (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param svc (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param sample.id (character) the id of the sample to be plotted
#' @param lrr.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param lrr.max (numeric) maximum CNV to be plotted
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param chrlist (character) vector containing chromosomes to plot; by default all chromosomes plotted
#' @keywords CNV, segmentation, structural variant, visualization, circular plot
#' @export
#' @examples
#' 
#' require(circlize)
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' ## select a random sample id
#' id <-  sample(intersect(cnv$sample,svc$sample))[1]
#' 
#' circ.wg.plot(cnv, svc, sample.id=id)


circ.wg.plot <- function(cnv, 
                         svc, 
                         sample.id=NULL,
                         lrr.pct = 0.2,
                         lrr.max = 4,
                         genome.v = "hg19",
                         chrlist=NULL){
  
  require(circlize, quietly = F, warn.conflicts = FALSE)

  cnvdat <- validate.cnv(cnv)
  svcdat <- validate.svc(svc)
  
  if(!is.null(sample.id)){ 
    cnvdat <- cnvdat[which(cnvdat$sample == sample.id),]
    svcdat <- svcdat[which(svcdat$sample == sample.id),]
    }
  
  stopifnot(length(unique(cnvdat$sample)) == 1)
  stopifnot(length(unique(svcdat$sample)) == 1)
  stopifnot(unique(cnvdat$sample) == unique(svcdat$sample) )
  
  sample_id <- unique(cnvdat$sample)

  if(is.null(chrlist)) chrlist <- paste("chr",c(1:22,"X","Y"),sep="")

  alllinks1 <- data.frame(svcdat$chrom1,svcdat$pos1,svcdat$pos1 )
  alllinks2 <- data.frame(svcdat$chrom2,svcdat$pos2,svcdat$pos2 )
  colnames(alllinks1) <- colnames(alllinks2) <- c("chr","start","end")
  map = setNames(c("blue", "red", "orange","black","green"), c("DEL", "DUP","INV","TRA","INS"))
  alllinkcolors <- map[as.character(svcdat$svclass)]
  
  cnv <- cnvdat[,c("chrom","start","end","segmean")]
  colores <- rep("black",nrow(cnv))
  colores[which(cnv$segmean < log2(1 - lrr.pct)) ] <- "blue"
  colores[which(cnv$segmean > log2(1 + lrr.pct)) ] <- "red"
  cnv <- remove.factors(data.frame(cnv,colores))
  cnv[which(cnv$segmean < log2(1/lrr.max) ),"segmean"] <- log2(1/lrr.max) 
  cnv[which(cnv$segmean > log2(lrr.max)),"segmean"] <- log2(lrr.max)
  allcnvlist <- list()
  for(i in chrlist) allcnvlist[[i]] <- cnv[which(cnv$chrom == i),]
  
  circos.initializeWithIdeogram(species=genome.v, chromosome.index=chrlist, plotType=c("ideogram","labels"))
  text(0, 0,  gsub("_","\n",sample_id), cex = 1)
  circos.genomicTrackPlotRegion(allcnvlist, bg.lwd =0.2, bg.col=rainbow(length(allcnvlist),alpha=0.1),ylim=c(-2.4,2.4), track.height=0.2, panel.fun = function(region, value, ...) {
    circos.genomicLines(region, value, col=as.character(allcnvlist[[CELL_META$sector.index]][,"colores"]), numeric.column = c(1), type="segment")
  })
  circos.genomicLink(alllinks1, alllinks2, col = alllinkcolors, border = NA)
}



