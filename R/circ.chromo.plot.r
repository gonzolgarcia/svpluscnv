#' 
#'
#' Circos plot combining segmentation and SV calls
#' 
#' @param chromo.regs2 (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param sample.id (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param lrr.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param lrr.max (numeric) maximum fold change to be plotted (default = 4)
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @keywords CNV, segmentation, structural variant, visualization, circular plot
#' @export
#' @examples
#' read.depth.breaks()



circ.chromo.plot <- function(chromo.regs2, 
                             sample.id,
                             lrr.pct = 0.2,
                             lrr.max= 4,
                             genome.v="hg19",
                             chr.list=NULL
                             ){
  
  require(circlize)
  library(taRifx)  # contains remove.factors

  segdat <- chromo.regs2$segdat[which(chromo.regs2$segdat$sample == sample.id),]
  svdat <- chromo.regs2$svdat[which(chromo.regs2$svdat$sample == sample.id),]
  regions <- remove.factors(chromo.regs2$regions.summary[[sample.id]])

  stopifnot(nrow(segdat) >  1)
  stopifnot(nrow(svdat) >  1)

  if(is.null(chr.list)) chr.list <- unique(regions$chrom)
  
  alllinks1 <- data.frame(svdat$chrom1,svdat$pos1,svdat$pos1 )
  alllinks2 <- data.frame(svdat$chrom2,svdat$pos2,svdat$pos2 )
  colnames(alllinks1) <- colnames(alllinks2) <- c("chr","start","end")
  map = setNames(c("blue", "red", "orange", "orange","orange","black","green"), c("DEL", "DUP","INV","h2hINV","t2tINV","TRA","INS"))
  alllinkcolors <- map[as.character(svdat$svclass)]
  
  cnv <- segdat[,c("chrom","start","end","segmean")]
  colores <- rep("black",nrow(cnv))
  colores[which(cnv$segmean < log2(1 - lrr.pct)) ] <- "blue"
  colores[which(cnv$segmean > log2(1 + lrr.pct)) ] <- "red"
  cnv <- remove.factors(data.frame(cnv,colores))
  cnv[which(cnv$segmean < log2(1/lrr.max) ),"segmean"] <- log2(1/lrr.max) 
  cnv[which(cnv$segmean > log2(lrr.max)),"segmean"] <- log2(lrr.max)
  allcnvlist <- list()
  for(i in chr.list) allcnvlist[[i]] <- cnv[which(cnv$chrom == i),]

  zoomchr <- intersect(which(alllinks1$chr %in% chr.list),which(alllinks2$chr %in% chr.list))
  links1<-alllinks1[zoomchr,]
  links2<-alllinks2[zoomchr,]
  linkcolors<-alllinkcolors[zoomchr]
  cnvlist <- list()
  for(i in chr.list) cnvlist[[i]] <- cnv[which(cnv$chrom == i),]
  reg.map = setNames(c("pink", "purple"), c("lc", "HC"))
  reg.col <- unname(reg.map[regions$conf])
  value <- rep(0.1,nrow(regions))
  regions.plot <- remove.factors(data.frame(regions,reg.col,value))
  
  p.regions <- list()
  for(chr in chr.list){  
    p.regions[[chr]] <- remove.factors(regions.plot[which(regions$chrom == chr),c("chrom","start","end","value","reg.col")])
    colnames(p.regions[[chr]]) <- c("chrom","start","end","value","color")
  }
  
    circos.initializeWithIdeogram(species=genome.v,chromosome.index=chr.list,plotType=c("axis","labels"), track.height=0.05, axis.labels.cex=0.4,labels.cex=1.3)
    circos.genomicIdeogram(track.height = 0.03)
    text(0, 0,  gsub("_","\n",sample.id), cex = 1.3)
    circos.genomicTrack(p.regions, bg.lwd =0.01, ylim=c(0,0.02), track.height=0.05,
                        panel.fun = function(region, value, ...) {
                          circos.genomicRect(region, value, ytop = 0.02, ybottom = 0, col = p.regions[[CELL_META$sector.index]][,"color"],  border = NA, ...)
                          circos.lines(CELL_META$cell.xlim, c(0.01, 0.01), lty = 2, col = "#00000040")
                        })
    
    circos.genomicTrackPlotRegion(cnvlist, bg.lwd =0.2, bg.col=rainbow(length(cnvlist),alpha=0.1),ylim=c(-2.5,2.5), track.height=0.2, 
                                  panel.fun = function(region, value, ...) {
                                    circos.genomicLines(region, value, col=as.character(cnvlist[[CELL_META$sector.index]][,"colores"]), numeric.column = c(1), type="segment")
                                  })
    circos.genomicLink(links1, links2, col = linkcolors, border = NA)
}

