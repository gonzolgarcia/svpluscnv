#' Circos plot combining segmentation and SV calls
#' @param chromo.regs.obj (chromo.regs) An object of class chromo.regs 
#' @param sample.id (character) the id of a sample to be plotted within 
#' @param lrr.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents 20 percent fold change
#' @param lrr.max (numeric) CNV plot limit
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param chrlist (character) vector containing chromosomes to plot; by default only chromosomes with shattered regions are ploted
#' @return circos plot into open device
#' @keywords CNV, segmentation, structural variant, visualization, circular plot
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' cnv <- validate.cnv(segdat_lung_ccle)
#' svc <- validate.svc(svdat_lung_ccle)
#' 
#' ## obtain shattered regions
#' shatt.regions <- shattered.regions(cnv,svc)
#' 
#' # select a random sample from the 
#' id <- "SCLC21H_LUNG"
#' 
#' circ.chromo.plot(shatt.regions, sample.id = id)

circ.chromo.plot <- function(chromo.regs.obj, 
                             sample.id,
                             genome.v = "hg19",
                             lrr.pct = 0.2,
                             lrr.max = 4,
                             chrlist=NULL,
                             ...){


if(sample.id %in% chromo.regs.obj@cnv@data$sample){
    cnvdat <- chromo.regs.obj@cnv@data[which(chromo.regs.obj@cnv@data$sample == sample.id),]
}
if(sample.id %in% chromo.regs.obj@svc@data$sample){
    svcdat <- chromo.regs.obj@svc@data[which(chromo.regs.obj@svc@data$sample == sample.id),]
}else{
    svcdat <- data.table()
}
regions <- chromo.regs.obj@regions.summary[[sample.id]]

stopifnot(nrow(chromo.regs.obj@cnv@data) > 0 | nrow(chromo.regs.obj@svc@data) >  0)

if(is.null(chrlist)) chrlist <- unique(regions$chrom)
  
if(nrow(svcdat) >  0){
    alllinks1 <- data.table(svcdat$chrom1,svcdat$pos1,svcdat$pos1 )
    alllinks2 <- data.table(svcdat$chrom2,svcdat$pos2,svcdat$pos2 )
    colnames(alllinks1) <- colnames(alllinks2) <- c("chr","start","end")
    map = setNames(c("blue", "red", "orange","black","green","black"), c("DEL", "DUP","INV","TRA","INS","BND"))
    alllinkcolors <- map[svcdat$svclass]
    zoomchr <- intersect(which(alllinks1$chr %in% chrlist),which(alllinks2$chr %in% chrlist))
    links1<-alllinks1[zoomchr,]
    links2<-alllinks2[zoomchr,]
    linkcolors<-alllinkcolors[zoomchr]
}

if(nrow(cnvdat) >  0){
    colores <- rep("black",nrow(cnvdat))
    colores[which(cnvdat$segmean < log2(1 - lrr.pct)) ] <- "blue"
    colores[which(cnvdat$segmean > log2(1 + lrr.pct)) ] <- "red"
    cnv.df <- data.frame(cnvdat[,c("chrom","start","end","segmean")],colores)
    cnv.df[,"colores"] <- as.character(cnv.df[,"colores"])
    cnv.df[which(cnv.df$segmean < log2(1/lrr.max) ),"segmean"] <- log2(1/lrr.max) 
    cnv.df[which(cnv.df$segmean > log2(lrr.max)),"segmean"] <- log2(lrr.max)
    allcnvlist <- list()
    for(i in chrlist) allcnvlist[[i]] <- cnv.df[which(cnv.df$chrom == i),]

    cnvlist <- list()
    for(i in chrlist) cnvlist[[i]] <- cnv.df[which(cnv.df$chrom == i),]
}

reg.map = setNames(c("pink", "purple"), c("lc", "HC"))
reg.col <- unname(reg.map[regions$conf])
value <- rep(0.1,nrow(regions))
regions.plot <- as.data.frame(data.table(regions,reg.col,value))
  
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
text(0, 0,  gsub("_","\n",sample.id),...)

}




#' Circos plot combining segmentation and SV calls
#' @param cnv (S4) an object of class svcnvio containing data type 'cnv' validated by validate.cnv
#' @param svc (S4) an object of class svcnvio containing data type 'svc' validated by validate.svc
#' @param sample.id (character) the id of the sample to be plotted
#' @param lrr.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param lrr.max (numeric) maximum CNV to be plotted
#' @param genome.v (hg19 or h38) reference genome version to draw chromosome limits and centromeres
#' @param chrlist (character) vector containing chromosomes to plot; by default all chromosomes plotted
#' @return circos plot into open device
#' @keywords CNV, segmentation, structural variant, visualization, circular plot
#' @export
#' @examples
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
    
    stopifnot(cnv@type == "cnv")
    cnvdat <- cnv@data
    
    stopifnot(svc@type == "svc")
    svcdat <- svc@data
    
    if(is.null(sample.id)){ 
        sample.id <- intersect(cnvdat$sample,svcdat$sample)
        stopifnot(length(sample.id) == 1)
    }
    cnvdat <- cnvdat[which(cnvdat$sample == sample.id),]
    svcdat <- svcdat[which(svcdat$sample == sample.id),]
    
    if(is.null(chrlist)) chrlist <- chr.sort(unique(cnvdat$chrom))
    
    alllinks1 <- data.table(svcdat$chrom1,svcdat$pos1,svcdat$pos1 )
    alllinks2 <- data.table(svcdat$chrom2,svcdat$pos2,svcdat$pos2 )
    colnames(alllinks1) <- colnames(alllinks2) <- c("chr","start","end")
    map = setNames(c("blue", "red", "orange","black","green","black"), c("DEL", "DUP","INV","TRA","INS","BND"))
    alllinkcolors <- map[as.character(svcdat$svclass)]
    
    cnvcirc <- cnvdat[,c("chrom","start","end","segmean")]
    colores <- rep("black",nrow(cnvcirc))
    colores[which(cnvcirc$segmean < log2(1 - lrr.pct)) ] <- "blue"
    colores[which(cnvcirc$segmean > log2(1 + lrr.pct)) ] <- "red"
    cnvcirc <- data.table(cnvcirc,colores)
    cnvcirc[which(cnvcirc$segmean < log2(1/lrr.max) ),"segmean"] <- log2(1/lrr.max) 
    cnvcirc[which(cnvcirc$segmean > log2(lrr.max)),"segmean"] <- log2(lrr.max)
    allcnvlist <- list()
    for(i in chrlist) allcnvlist[[i]] <- as.data.frame(cnvcirc[which(cnvcirc$chrom == i),])
    
    circos.initializeWithIdeogram(species=genome.v, chromosome.index=chrlist, plotType=c("ideogram","labels"))
    text(0, 0,  gsub("_","\n",sample.id), cex = 1)
    circos.genomicTrackPlotRegion(allcnvlist, bg.lwd =0.2, bg.col=rainbow(length(allcnvlist),alpha=0.1),ylim=c(-2.4,2.4), track.height=0.2, panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value, col=as.character(allcnvlist[[CELL_META$sector.index]][,"colores"]), numeric.column = c(1), type="segment")
    })
    circos.genomicLink(alllinks1, alllinks2, col = alllinkcolors, border = NA)
}






