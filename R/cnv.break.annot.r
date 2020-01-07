#' 
#'
#' Identify recurrently altered genes
#' @param segdf (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param chr.lim (data.frame) 3 column table (chrom, begin, end) returned by chromosome.limit.coords
#' @param window.size (numeric) size in megabases of the genmome bin to compute break density 
#' @param slide.size (numeric) size in megabases of the sliding genmome window
#' @keywords CNV, segmentation
#' @export
#' @examples
#' cnv.break.annot()

cnv.break.annot <- function(segdf, fc.pct = 0.2, 
                        genome.v="hg19",
                        maxgap= 1000,
                        upstr=50000,
                        min.seg.size = NULL, 
                        min.num.probes = NULL, 
                        low.cov = NULL,
                        clean.brk = NULL,
					    verbose=TRUE){
	
    segdat <- validate.seg(segdf)
    chrlist <- unique(segdat$chrom)
    require(org.Hs.eg.db)
    
    if(genome.v %in% c("hg19","GRCh37")){
		require(TxDb.Hsapiens.UCSC.hg19.knownGene)
		genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
	}else if(genome.v %in% c("hg38","GRCh38")){
		require(TxDb.Hsapiens.UCSC.hg38.knownGene)
		genesgr = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
	}else{stop("Unspecified, or non available genome")}
	
    genesgr = genesgr[which(as.character(genesgr@seqnames) %in% chrlist)]
    genesgr@elementMetadata$SYMBOL <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID')
    
    genesgr <- genesgr[which(!is.na(genesgr$SYMBOL)),]
    genesgr <- genesgr[which(lapply(genesgr@elementMetadata$SYMBOL,length) > 0),]
    
    genes_tab <- as.data.frame(genesgr)
    rownames(genes_tab) <- genes_tab$SYMBOL
    
    genes_ups <- genes_tab
    genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - maxgap
    genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr - maxgap
    genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + maxgap
    genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr + maxgap
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id, SYMBOL=SYMBOL))
    upstr_tab <- as.data.frame(upstreamgr)
    rownames(upstr_tab) <- upstr_tab$SYMBOL
    
	cnv_breaks <- seg.breaks(segdat, fc.pct = fc.pct, min.seg.size = min.seg.size, min.num.probes = min.num.probes, 
						 low.cov = low.cov, clean.brk = clean.brk,verbose=verbose)
	breaks<-cnv_breaks$breaks
	rownames(breaks) <- gsub(" ","",apply(breaks[,1:3],1,paste,collapse="-"))
	
	breaksgr = with(breaks[,c("chrom","start","end")], GRanges(chrom, IRanges(start=start, end=end)))
	
	overlapGenes <- GenomicAlignments::findOverlaps(genesgr,breaksgr,ignore.strand=TRUE,maxgap = maxgap)
	
	overlapUpstream <- GenomicAlignments::findOverlaps(upstreamgr,breaksgr,ignore.strand=TRUE)

	geneHits <- cbind(rownames(breaks)[subjectHits(overlapGenes)],genes_tab[queryHits(overlapGenes),"SYMBOL"])
	upstreamHits <- cbind(rownames(breaks)[subjectHits(overlapUpstream)],genes_ups[queryHits(overlapUpstream),"SYMBOL"])
	
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

	return(list(
	    breaks= breaks,
	    genes = genes_tab,
	    disruptSamples = geneSamples,
	    disruptBreaks = geneBreaks,
	    upstreamSamples = upstreamSamples,
	    upstreamBreaks = upstreamBreaks,
	    param = list(
	        genome.v = genome.v,
	        maxgap = maxgap, 
	        upstr = upstr, 
	        verbose = verbose)
	))
	
}

