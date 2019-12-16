#' 
#'
#' Identify recurrently altered genes
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) returned by chromosome.limit.coords
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @keywords CNV, segmentation
#' @export
#' @examples
#' cnv.gene.breaks()

cnv.gene.breaks <- function(seg, fc.pct = 0.2, 
                        genome.v="hg19",
                        geneid="Symbol",
                        maxgap=1000,
                        upstr=50000,
                        min.seg.size = NULL, 
                        min.num.probes = NULL, 
                        low.cov = NULL,
                        clean.brk = NULL,
					    verbose=TRUE){
	
    segdf <- validate.seg(seg)
	if(genome.v %in% c("hg19","GRCh37")){
		require(TxDb.Hsapiens.UCSC.hg19.knownGene)
		require(org.Hs.eg.db)
		genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
		if(geneid== "Symbol") genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID')
		genesgr <- genesgr[which(!is.na(genesgr$gene_id)),]
	}else if(genome.v %in% c("hg38","GRCh38")){
		require(TxDb.Hsapiens.UCSC.hg38.knownGene)
		require(org.Hs.eg.db)
		genesgr = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
		if(geneid== "Symbol") genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID')
		genesgr <- genesgr[which(!is.na(genesgr$gene_id)),]
	}else{stop("Unspecified, or non available genome")}
	
	genes_tab <- as.data.frame(genesgr)
	genes_ups <- genes_tab
	genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"]
	genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr
	genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"]
	genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr
	upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start=start, end=end), strand=strand, gene_id=gene_id))
	
	breaks <- seg.breaks(segdf, fc.pct = fc.pct, min.seg.size = min.seg.size, min.num.probes = min.num.probes, 
						 low.cov = low.cov, clean.brk = clean.brk,verbose=verbose)
	rownames(breaks) <- gsub(" ","",apply(breaks[,1:3],1,paste,collapse="-"))
	
	breaksgr = with(breaks[,c("chrom","start","end")], GRanges(chrom, IRanges(start=start, end=end)))
	
	overlapGenes <- GenomicAlignments::findOverlaps(genesgr,breaksgr,ignore.strand=TRUE,maxgap = maxgap)
	
	overlapUpstream <- GenomicAlignments::findOverlaps(upstreamgr,breaksgr,ignore.strand=TRUE)

	geneHits <- cbind(rownames(breaks)[subjectHits(overlapGenes)],genes_tab[queryHits(overlapGenes),"gene_id"])
	upstreamHits <- cbind(rownames(breaks)[subjectHits(overlapUpstream)],genes_ups[queryHits(overlapUpstream),"gene_id"])
	
	if(verbose){
	    pb <- txtProgressBar(style=3)
	    cc <-0
	    tot <- length(c(unique(geneHits[,2]),unique(geneHits[,2])))
	}
	
	geneBreaks<-list()
	for(i in unique(geneHits[,2])){
	    geneBreaks[[i]] <- unique(geneHits[which(geneHits[,2] == i),1])
	    if(verbose) cc <- cc+1
	    if(verbose) setTxtProgressBar(pb, cc/tot)
	}
	geneSamples <- lapply(geneBreaks,function(x) unique(breaks[x,"sample"]) )

	upstreamBreaks <-list()
	for(i in unique(upstreamHits[,2])){
	    upstreamBreaks[[i]] <- unique(upstreamHits[which(upstreamHits[,2] == i),1])
	    if(verbose) cc <- cc+1
	    if(verbose) setTxtProgressBar(pb, cc/tot)
	}
	upstreamSamples <- lapply(upstreamBreaks,function(x) unique(breaks[x,"sample"]) )

	if(verbose) close(pb)

	
	results<-list()
	results[["breaks"]] <-  breaks
	results[["feature_tab"]] <-  genesgr
	results[["upstreamBreaks"]] <-  upstreamBreaks
	results[["upstreamSamples"]] <-  upstreamSamples
	results[["geneBreaks"]] <-  geneBreaks
	results[["geneSamples"]] <-  geneSamples
	
	return(results)
}

