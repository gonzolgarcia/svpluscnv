#' Define peak regions and samples in shattered regions hot spots
#' @param chromo.regs.obj (chromo.regs) An object of class chromo.regs 
#' @param freq.cut (numeric) the hot spot threshold above which peaks are defined for sample ID retrieval
#' @return a list comprising two lists: peakRegions, peakRegionsSamples
#' @export
#' 

hot.spot.samples <- function(chromo.regs.obj, freq.cut){

freq.matrix <- apply(chromo.regs.obj@high.density.regions.hc,2,sum)
textRegions <- names(which(freq.matrix >= freq.cut))
hitRegions <- data.table(do.call(rbind,strsplit(textRegions," ")),textRegions)
colnames(hitRegions) <- c("chr","start","end","regid")
hitRegions$start <- as.numeric(hitRegions$start)
hitRegions$end <- as.numeric(hitRegions$end)


# collapes contiguous bins into unique regions
bins2remove <- c()
for(i in 2:nrow(hitRegions)){ 
    if(hitRegions[i]$chr == hitRegions[i-1]$chr){
        if(hitRegions[i]$start < (hitRegions[i-1]$end)){
            hitRegions[i]$start <- hitRegions[i-1]$start
            bins2remove <- c(bins2remove,textRegions[i-1])
        }
    }
}
hitRegionsPost<- hitRegions[which(hitRegions$regid %in% setdiff(hitRegions$regid,bins2remove))]

hitRegions_gr <- with(hitRegions, GRanges(chr, IRanges(start=start, end=end)))
hitRegionsPost_gr <- with(hitRegionsPost, GRanges(chr, IRanges(start=start, end=end)))
hits <-GenomicAlignments::findOverlaps(hitRegionsPost_gr,hitRegions_gr)

regList <- list()
for(i in unique(queryHits(hits))) regList[[hitRegionsPost[i]$regid]] <- textRegions[subjectHits(hits)[which(queryHits(hits) == i)]]

# obtain the genomic bins with maximum number of samples
peakRegions <- lapply(regList, function(x) 
    names(which(freq.matrix[x] == max(freq.matrix[x]))))

# collect samples with shattered region in the peaks 
peakRegionsSamples <- lapply(peakRegions, function(x) 
    names(which(apply(cbind(chromo.regs.obj@high.density.regions.hc[,x]),1,sum) > 0)))

return(list(peakRegions=peakRegions,peakRegionsSamples=peakRegionsSamples))

}

