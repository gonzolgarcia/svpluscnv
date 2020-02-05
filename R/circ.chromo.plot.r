#' Circos plot combining segmentation and SV calls
#' @param chromo.regs.obj (chromo.regs) An object of class chromo.regs 
#' @param sample.id (character) the id of a sample to be plotted within names(chromo.regs.obj@regions.summ)
#' @param lrr.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents 20 percent fold change
#' @param lrr.max (numeric) CNV plot limit
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param chrlist (character) vector containing chromosomes to plot; by default only chromosomes with shattered regions are ploted
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
#' ## obtain shattered regions
#' shatt.regions <- shattered.regions(cnv,svc)
#' 
#' id <-  names(chromo.regs.obj@regions.summary)[1]
#' circ.chromo.plot(shatt.regions, sample.id = id)

circ.chromo.plot <- function(chromo.regs.obj, 
                             sample.id,
                             genome.v = "hg19",
                             lrr.pct = 0.2,
                             lrr.max = 4,
                             chrlist=NULL
                             ){

require(circlize)

cnvdat <- chromo.regs.obj@cnvdat[which(chromo.regs.obj@cnvdat$sample == sample.id),]
svcdat <- chromo.regs.obj@svcdat[which(chromo.regs.obj@svcdat$sample == sample.id),]
regions <- chromo.regs.obj@regions.summary[[sample.id]]

stopifnot(nrow(cnvdat) > 0 | nrow(svcdat) >  0)

if(is.null(chrlist)) chrlist <- unique(regions$chrom)
  
if(nrow(svcdat) >  0){
    alllinks1 <- data.frame(svcdat$chrom1,svcdat$pos1,svcdat$pos1 )
    alllinks2 <- data.frame(svcdat$chrom2,svcdat$pos2,svcdat$pos2 )
    colnames(alllinks1) <- colnames(alllinks2) <- c("chr","start","end")
    map = setNames(c("blue", "red", "orange","black","green"), c("DEL", "DUP","INV","TRA","INS"))
    alllinkcolors <- map[as.character(svcdat$svclass)]
    zoomchr <- intersect(which(alllinks1$chr %in% chrlist),which(alllinks2$chr %in% chrlist))
    links1<-alllinks1[zoomchr,]
    links2<-alllinks2[zoomchr,]
    linkcolors<-alllinkcolors[zoomchr]
}

if(nrow(cnvdat) >  0){
    cnv <- cnvdat[,c("chrom","start","end","segmean")]
    colores <- rep("black",nrow(cnv))
    colores[which(cnv$segmean < log2(1 - lrr.pct)) ] <- "blue"
    colores[which(cnv$segmean > log2(1 + lrr.pct)) ] <- "red"
    cnv <- remove.factors(data.frame(cnv,colores))
    cnv[which(cnv$segmean < log2(1/lrr.max) ),"segmean"] <- log2(1/lrr.max) 
    cnv[which(cnv$segmean > log2(lrr.max)),"segmean"] <- log2(lrr.max)
    allcnvlist <- list()
    for(i in chrlist) allcnvlist[[i]] <- cnv[which(cnv$chrom == i),]

    cnvlist <- list()
    for(i in chrlist) cnvlist[[i]] <- cnv[which(cnv$chrom == i),]
}

reg.map = setNames(c("pink", "purple"), c("lc", "HC"))
reg.col <- unname(reg.map[regions$conf])
value <- rep(0.1,nrow(regions))
regions.plot <- remove.factors(data.frame(regions,reg.col,value))
  
p.regions <- list()
for(chr in chrlist){
    p.regions[[chr]] <- regions.plot[which(regions$chrom == chr),c("chrom","start","end","value","reg.col")]
    colnames(p.regions[[chr]]) <- c("chrom","start","end","value","color")
}
  
circos.initializeWithIdeogram(species=genome.v,chromosome.index=chrlist,plotType=c("axis","labels"), track.height=0.05, axis.labels.cex=0.4,labels.cex=1.3)
circos.genomicIdeogram(track.height = 0.03)
circos.genomicTrack(p.regions, bg.lwd =0.01, ylim=c(0,0.02), track.height=0.05,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0.02, ybottom = 0, col = p.regions[[CELL_META$sector.index]][,"color"],  border = NA, ...)
                      circos.lines(CELL_META$cell.xlim, c(0.01, 0.01), lty = 2, col = "#00000040")
                    })
    
circos.genomicTrackPlotRegion(cnvlist, bg.lwd =0.2, bg.col=rainbow(length(cnvlist),alpha=0.1),ylim=c(-2.5,2.5), track.height=0.2, 
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, col=as.character(cnvlist[[CELL_META$sector.index]][,"colores"]), numeric.column = c(1), type="segment")
                              })
if(nrow(svcdat) >  0) circos.genomicLink(links1, links2, col = linkcolors, border = NA)
text(0, 0,  gsub("_","\n",sample.id), cex = 1.3)

}

