#' 
#'
#' Define peak regions and samples in shattered regions hot spots
#' @param chromo.regs (list) an object returned by 'shattered.regions' and 'shattered.regions.cnv'
#' @param freq.cut (numeric) the hot spot threshold above which peaks are defined and sample ids collected 
#' @export
#' 

hot.spot.samples <- function(chromo.regs,freq.cut){

freq.matrix <- apply(chromo.regs$high.density.regions.hc,2,sum)
textRegions <- names(which(freq.matrix >= freq.cut))
hitRegions <- remove.factors((data.frame(do.call(rbind,strsplit(textRegions," ")))))
hitRegions[,2] <- as.numeric(hitRegions[,2])
hitRegions[,3] <- as.numeric(hitRegions[,3])
colnames(hitRegions) <- c("chr","start","end")
rownames(hitRegions) <-textRegions

# collapes contiguous bins into unique regions
bins2remove <- c()
for(i in 2:nrow(hitRegions)){ 
    if(hitRegions[i,"chr"] == hitRegions[i-1,"chr"] ){
        if(hitRegions[i,"start"] < (hitRegions[i-1,"end"])){
            hitRegions[i,"start"] <- hitRegions[i-1,"start"]
            bins2remove <- c(bins2remove,textRegions[i-1])
        }
    }
}
hitRegionsPost<- hitRegions[setdiff(rownames(hitRegions),bins2remove),]

require(GenomicRanges)
hitRegions_gr <- with(hitRegions, GRanges(chr, IRanges(start=start, end=end)))
hitRegionsPost_gr <- with(hitRegionsPost, GRanges(chr, IRanges(start=start, end=end)))
hits <-GenomicAlignments::findOverlaps(hitRegionsPost_gr,hitRegions_gr)

regList <- list()
for(i in unique(queryHits(hits))) regList[[paste(hitRegionsPost[i,],collapse=" ") ]] <- textRegions[subjectHits(hits)[which(queryHits(hits) == i)]]

# obtain the genomic bins with maximum number of samples
peakRegions <- lapply(regList, function(x) 
    names(which(freq.matrix[x] == max(freq.matrix[x]))))

# collect samples with shattered region in the peaks 
peakRegionsSamples <- lapply(peakRegions, function(x) 
    names(which(apply(cbind(shatt_lung_cnv$high.density.regions.hc[,x]),1,sum) > 0)))

return(list(peakRegions=peakRegions,peakRegionsSamples=peakRegionsSamples))

}

