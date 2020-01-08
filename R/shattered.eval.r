#' 
#'
#' Evaluate shattered regions
#' @param chromo.regs (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param interleaved.cut (numeric) the percentage of non interleaved structural variant calls 
#' @param dist.iqm.cut (numeric) interquantile average of the distance between breakpoints within a shattered region
#' @param verbose (logical)
#' @keywords chromothripsis, segmentation, structural variants
#' @export
#' @examples
#' shattered.eval()



shattered.eval <- function(chromo.regs,
                           interleaved.cut=0.5,
                           dist.iqm.cut=100000,
                           verbose=TRUE){
  
  svdat <- chromo.regs$svdat
  segbrk  <- chromo.regs$segbrk
  svbrk <- chromo.regs$svbrk

  if(verbose){
    message(paste("Evaluating shattered regions in ",length(names(chromo.regs$regions.summary))," samples...",sep="") )
    pb <- txtProgressBar(style=3)
    cc <-0
    tot <- length(chromo.regs$regions.summary)
  }
  
  for(cl in names(chromo.regs$regions.summary)){

    if(verbose) cc <- cc+1
    if(verbose) setTxtProgressBar(pb, cc/tot)

    regions <-   chromo.regs$regions.summary[[cl]]
    br1 <- segbrk$breaks[which(segbrk$breaks$sample == cl),2:3]
    br2 <- svbrk$breaks[which(svbrk$breaks$sample == cl),2:3]
    colnames(br1) <- colnames(br2) <- c("chrom","pos")
    br1.gr <- with(br1, GRanges(chrom, IRanges(start=pos, end=pos)))
    br2.gr <- with(br2, GRanges(chrom, IRanges(start=pos, end=pos)))
    
    c.br1 <- chromo.regs$common.brk$brk1_match[which(chromo.regs$common.brk$brk1_match$sample == cl),2:3]
    c.br2 <- chromo.regs$common.brk$brk2_match[which(chromo.regs$common.brk$brk2_match$sample == cl),2:3]
    colnames(c.br1) <- colnames(c.br2) <- c("chrom","pos")
    c.br1.gr <- with(c.br1, GRanges(chrom, IRanges(start=pos, end=pos)))
    c.br2.gr <- with(c.br2, GRanges(chrom, IRanges(start=pos, end=pos)))
    
    regions_gr <- with(regions, GRanges(chrom, IRanges(start=start, end=end)))
    hits_1 = GenomicAlignments::findOverlaps(regions_gr,br1.gr)
    hits_2 = GenomicAlignments::findOverlaps(regions_gr,br2.gr)
    
    c.hits_1 = GenomicAlignments::findOverlaps(regions_gr,c.br1.gr)
    c.hits_2 = GenomicAlignments::findOverlaps(regions_gr,c.br2.gr)
    
    svdat_i <- svdat[which(svdat$sample == cl),]
    sv_ranges_ori <-   with(svdat_i, GRanges(chrom1, IRanges(start=pos1, end=pos1)))
    sv_ranges_dest <-   with(svdat_i, GRanges(chrom2, IRanges(start=pos2, end=pos2)))
    hits_ori = GenomicAlignments::findOverlaps(regions_gr,sv_ranges_ori)
    hits_dest = GenomicAlignments::findOverlaps(regions_gr,sv_ranges_dest)
    
    # calculate interleaved score, n.breaks, dispersion, and median distance between consecutive breaks
    start <- end <- reg.size <- dist.iqm.seg <- dist.iqm.sv <- n.brk.seg <- n.brk.sv <- n.orth.seg <- n.orth.sv <- interleaved <- rep(0,nrow(regions))
    conf <- rep("HC",nrow(regions))
    for(i in 1:nrow(regions)){
      sites1 <- sort(unique(br1[subjectHits(hits_1)[which(queryHits(hits_1) == i)],"pos"]))
      dist.iqm.seg[i]  <- IQM(sites1[2:length(sites1)] - sites1[1:(length(sites1)-1) ],lowQ = 0.2,upQ = 0.8)
      n.brk.seg[i] <- length(sites1)  
      
      sites2 <- sort(unique(br2[subjectHits(hits_2)[which(queryHits(hits_2) == i)],"pos"]))
      dist.iqm.sv[i]  <- IQM(sites2[2:length(sites2)] - sites2[1:(length(sites2)-1) ],lowQ = 0.2,upQ = 0.8)
      n.brk.sv[i] <- length(sites2)  
      
      c.sites1 <- sort(unique(c.br1[subjectHits(c.hits_1)[which(queryHits(c.hits_1) == i)],"pos"]))
      n.orth.seg[i] <- length(c.sites1)  

      c.sites2 <- sort(unique(c.br2[subjectHits(c.hits_2)[which(queryHits(c.hits_2) == i)],"pos"]))
      n.orth.sv[i] <- length(c.sites2)  
      
      start[i] <- min(c(sites1,sites2))
      end[i] <- max(c(sites1,sites2))
      
      a<-svdat_i[subjectHits(hits_ori)[which(queryHits(hits_ori) == i)],]
      b<-svdat_i[subjectHits(hits_dest)[which(queryHits(hits_dest) == i)],]
      apos<-a$pos1
      names(apos) <- rownames(a)
      bpos<-b$pos2
      names(bpos) <- rownames(b)
      brkids <- names(sort(c(apos,bpos)))
      reg.size <- regions[,"end"]-regions[,"start"]
      ct<-0
      for(x in 2:length(brkids)) if(brkids[x] == brkids[x-1]) ct<-ct+1
      interleaved[i] <- ct/length(unique(brkids))
    }    
    chrom <-regions[,"chrom"]
    nseg<-regions[,"nseg"]
    regions <- data.frame(chrom,start,end,nseg)
    
    # Find links between shattered regions   
    if(nrow(chromo.regs$regions.summary[[cl]]) > 1 ){
      record_mat<- matrix(nrow=0,ncol=2)
      links <- rep("",nrow(regions))
      for(i in 1:nrow(regions)){
        for(j in 1:nrow(regions)){
          region_a_hits <- subjectHits(hits_ori)[which(queryHits(hits_ori) == i)]
          region_b_hits <- subjectHits(hits_dest)[which(queryHits(hits_dest) == j)]
          if(length(intersect(region_a_hits,region_b_hits)) > 0 ){
            record_mat <- rbind(record_mat,c(i,j),c(j,i))
            }
          }
        }
      for(i in 1:nrow(regions)){
        links[i] <- paste(setdiff(as.character(sort(unique(c(i, record_mat[which(record_mat[,1] == i),2]) ))),as.character(i)),collapse=",")
        if(length(grep("[0-9]",links[i])) == 0) links[i]<-"-"
        }
    }else{
      links <- "-"
    }

  if(!is.null(interleaved.cut)) conf[which(interleaved >= interleaved.cut)] <- "lc"
  if(!is.null(dist.iqm.cut)) conf[which(apply(cbind(dist.iqm.seg,dist.iqm.sv),1,mean) < dist.iqm.cut)] <- "lc"
  conf[which(links != "-")] <- "HC"
  chromo.regs$regions.summary[[cl]] <- remove.factors(data.frame(regions,links,reg.size,dist.iqm.seg,dist.iqm.sv,n.brk.seg,n.brk.sv,n.orth.seg,n.orth.sv,interleaved,conf))
    
  }
  if(verbose) close(pb)
  
  return(chromo.regs)
}

