#' 
#' Integrated visualization of SVs and SNV in local genomic regions
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param chr (character) chromosome (e.g chr9)
#' @param start (numeric) genomic coordinate from specified chromosome to start plotting
#' @param end (numeric) genomic coordinate from specified chromosome to stop plotting
#' @param cnvlim (numeric) limits for color coding of background CNV log-ratios
#' @param addlegend (character) one of 'sv' (show SV type legend), 'cnv' (show CNV background color legend) or 'both'.
#' @param addtext (character) a vector indicating what SV types should include text labels indicating brakpoint partners genomic locations. 
#' The labeels are only added to point tbreakpoint locations outside the plot area. (e.g. c("TRA","INV") )
#' @param text.size (numeric) point size of the text 
#' @param ... additional plot parameters from graphics plot function 
#' @keywords CNV, segmentation
#' @export
#' @examples
#' plot()


sv.model.view <- function(sv, seg, chr, start, stop, 
                          cnvlim=c(-2,2), 
                          addlegend='both',
                          interval=NULL,
                          addtext=NULL,
                          text.size=.8,
                          ...){
    

    stopifnot(!is.null(chr) && !is.null(start) && !is.null(stop))
    
    svdat <- validate.sv(sv)
    segdat <- validate.seg(seg)

    
    sv1gr = with(svdat, GRanges(chrom1, IRanges(start=pos1, end=pos1))) 
    sv2gr = with(svdat, GRanges(chrom2, IRanges(start=pos2, end=pos2))) 
    genegr <- with(data.frame(chr,start,stop), GRanges(chr, IRanges(start=start, end=stop))) 
    
    hits1 = GenomicAlignments::findOverlaps(sv1gr,genegr)
    hits2 = GenomicAlignments::findOverlaps(sv2gr,genegr)
    svtab <- svdat[sort(unique(c(queryHits(hits1),queryHits(hits2)))),]
    svcolormap = setNames(c("blue", "red", "orange", "black", "green","grey20"), 
                   c("DEL", "DUP", "INV", "TRA", "INS","BND"))
    svcolor <- svcolormap[svtab$svclass]
    svtab_plot <- remove.factors(data.frame(svtab,svcolor))
    svtab_plot_seg <- svtab_plot[which(svtab_plot$svclass != "TRA"),]
    svtab_plot_tra <- svtab_plot[which(svtab_plot$svclass == "TRA"),]
    
    segwithsv <- segdat[which(segdat$sample %in% svtab$sample),]
    seggr <- with(segwithsv, GRanges(chrom, IRanges(start=start, end=end))) 
    hits_seg = GenomicAlignments::findOverlaps(seggr,genegr)
    seg_plot <- segwithsv[queryHits(hits_seg),]
    segcolor <- map2color(seg_plot$segmean,
              pal=colorRampPalette(c("lightblue","white","salmon"))(256),
              limit=cnvlim)
    seg_plot <- remove.factors(data.frame(seg_plot,segcolor))
        
    sample_order <- 1:length(unique(svtab[,"sample"]))
    names(sample_order) <- unique(svtab[,"sample"])

    if(!is.null(addlegend)){
        plot_ylim <- length(sample_order)*10/100+length(sample_order)
        legend_ypos <- plot_ylim - length(sample_order)*3/100 
        if(length(sample_order) < 10) plot_ylim <- length(sample_order) +1
    }else{
        plot_ylim <- length(sample_order)
    }
    
    plot(x=NULL,y=NULL,xlim=range(c(start,stop)),ylim=range(c(0,plot_ylim)),
         xaxt='n',yaxt='n',xlab='',ylab='',bty='n', ...)

    mtext(side=2,at=sample_order-0.5,text=names(sample_order),las=2,line = 0.5, ...)
    
    for(sid in names(sample_order)){
        ypos <- sample_order[sid]
        polygon(rbind(
            c(start-1e7,ypos+0.02),
            c(start-1e7,ypos-0.98),
            c(stop+1e7,ypos-0.98),
            c(stop+1e7,ypos+0.02)),
            col=rep(c("grey80","grey80"),length(sample_order))[ypos],border=NA)
    }
        
    for(sid in names(sample_order)){
        seg_sample_plot <- seg_plot[which(seg_plot$sample == sid),]
        ypos <- sample_order[sid]
        for(i in 1:nrow(seg_sample_plot)){
            polygon(rbind(
                c(seg_sample_plot[i,"start"],ypos),
                c(seg_sample_plot[i,"start"],ypos-1),
                c(seg_sample_plot[i,"end"],ypos-1),
                c(seg_sample_plot[i,"end"],ypos)
            ),col=seg_sample_plot[i,"segcolor"],border=NA)
        }
    }
    
    for(sid in unique(svtab_plot_seg$sample)){
        svtab_plot_seg_i <- svtab_plot_seg[which(svtab_plot_seg$sample == sid),]
        ypos <- sample_order[sid]
        addrnorm <- rep(c(0,0.2,-0.2,0.1,-0.1,0.3,-0.3),nrow(svtab_plot_seg_i))
        for(i in 1:nrow(svtab_plot_seg_i)){
            polygon(rbind(
                c(svtab_plot_seg_i[i,"pos1"],ypos-0.4-addrnorm[i]),
                c(svtab_plot_seg_i[i,"pos1"],ypos-0.6-addrnorm[i]),
                c(svtab_plot_seg_i[i,"pos2"],ypos-0.6-addrnorm[i]),
                c(svtab_plot_seg_i[i,"pos2"],ypos-0.4-addrnorm[i])
            ),col=NA,border=svtab_plot_seg_i[i,"svcolor"])

            if(svtab_plot_seg_i[i,"svclass"] %in% addtext){
                if(svtab_plot_seg_i[i,"pos1"] < start){
                    text(start,ypos-0.5-addrnorm[i],
                         paste("<--",svtab_plot_seg_i[i,"pos1"],sep=""),
                         pos=4,offset=0,cex=text.size)
                }
                if(svtab_plot_seg_i[i,"pos2"] > stop){
                    text(stop,ypos-0.5-addrnorm[i],
                         paste(svtab_plot_seg_i[i,"pos2"],"-->",sep=""),
                         pos=2,offset=0,cex=text.size)
                }
            }
        }
    }

    for(sid in unique(svtab_plot_tra$sample)){
        svtab_plot_tra_i <- svtab_plot_tra[which(svtab_plot_tra$sample == sid),]
        ypos <- sample_order[sid]
        addrnorm <- rep(c(0,0.3,-0.3,0.1,-0.1,0.2,-0.2),nrow(svtab_plot_seg_i))
        for(i in 1:nrow(svtab_plot_tra_i)){
            if(svtab_plot_tra_i[i,"chrom2"] == chr){ 
                points(svtab_plot_tra_i[i,"pos2"],ypos-0.5+addrnorm[i],pch=10)
                lines(c(svtab_plot_tra_i[i,"pos2"],svtab_plot_tra_i[i,"pos2"]),c(ypos,ypos-1),lwd=1,lty=3)
                if("TRA" %in% addtext){
                    text(svtab_plot_tra_i[i,"pos2"],ypos-0.5+addrnorm[i],
                         paste("-->",svtab_plot_tra_i[i,"chrom1"],":",svtab_plot_tra_i[i,"pos1"],sep=""),
                         pos=4,offset=0,cex=text.size)
                }
            }            
            if(svtab_plot_tra_i[i,"chrom1"] == chr){
                points(svtab_plot_tra_i[i,"pos1"],ypos-0.5+addrnorm[i],pch=10)
                lines(c(svtab_plot_tra_i[i,"pos1"],svtab_plot_tra_i[i,"pos1"]),c(ypos,ypos-1),lwd=1,lty=3)
                if("TRA" %in% addtext) {
                    text(svtab_plot_tra_i[i,"pos1"],ypos-0.5+addrnorm[i],
                         paste("-->",svtab_plot_tra_i[i,"chrom2"],":",svtab_plot_tra_i[i,"pos2"],sep=""),
                         pos=4,offset=0,cex=text.size)
                }
            }
        }
    }
    
    if(is.null(interval)) interval <- round((stop - start)/5000) * 1000
    xlabs <- seq(floor(start/10000)*10000, ceiling(stop/10000)*10000,interval)
    axis(1, at = xlabs,labels=TRUE, lwd.ticks=1.5, pos=0,...)
    #mtext(xlabs, at=xlabs, side=1, cex=cex.axis)
    
    if(addlegend %in% c("sv","both")) {
        legend(x= start, y =legend_ypos, legend = c("DEL", "DUP", "INV","BND", "TRA"), 
               bty = "n", fill = c("white", "white", "white", "white",NA), border=c("blue", "red","orange","grey20",NA), 
               pch = c(NA, NA,NA, NA,10), horiz = TRUE, x.intersp=0.2)
    }
    if(addlegend %in% c("cnv","both")) {
        legend(x=start + (stop-start)/2, y = legend_ypos,legend = c(paste("CNV= ",cnvlim[1],sep=""), "CNV= 0", paste("CNV= ",cnvlim[2],sep=""), "no-data"),
               bty = "n",fill=c("lightblue","white","salmon","grey80"), border=c(NA,"black",NA,NA), 
               horiz = TRUE, x.intersp=0.2)
    }
    #par()
}


    
