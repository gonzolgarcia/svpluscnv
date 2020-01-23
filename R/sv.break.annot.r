#' Class to store breakpoint annotations
#' @param input (data.frame): chromosome for the first breakpoint
#' @param geneinfo (data.frame): position for the first breakpoint
#' @param disruptSamples (list): chromosome for the second breakpoint
#' @param disruptBreaks (list): position for the second breakpoint
#' @param upstreamSamples (list):
#' @param upstreamBreaks (list):
#' @param param (list):
#' @return an instance of the class 'break.annot' containing breakpoint mapping onto genes
#' @export
break.annot <- setClass("break.annot",
                   representation(
                       input  = "data.frame",
                       geneinfo = "data.frame",
                       disruptSamples = "list",
                       disruptBreaks = "list",
                       upstreamSamples = "list",
                       upstreamBreaks = "list",
                       param = "list"
                   ))


setMethod("show","break.annot",function(object){
    writeLines(paste("An object of class annotsv from svcnvplus containing the following stats:
                \nNumber of samples=",length(unique(object@input$sample)),
                "\nAltered genes=",length(object@disruptSamples)))
})



#' Identify recurrently altered genes by strutural variants
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param maxgap (numeric) offset region size (base pairs) beyond transcription start and end sites to extend for overlaps
#' @param upstr (numeric) upstream region size (base pairs) to identify breakpoints in gene upstream regulatory regions 
#' @param sv.seg.size (numeric) segmental variants (DEL, DUP, INV or INS) size limit; bigger sigments are considered translocations
#' @return an instance of the class 'break.annot' containing breakpoint mapping onto genes
#' @keywords Structural variants, annotation
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from SV calls data
#' sv <- validate.sv(svdat_lung_ccle)
#' 
#' sv.break.annot(sv, genome.v="hg19")

sv.break.annot <- function(sv, 
                     genome.v="hg19",
                     maxgap=0,
                     upstr=50000, 
                     sv.seg.size=200000,
                     verbose=TRUE){
    
    svdat <- validate.sv(sv)
    chrlist <- unique(c(svdat$chrom1,svdat$chrom2))
    
    rownames(svdat) <- createRandomString(nrow(svdat),8)
    
    if(verbose) message("# Generating all ranges objects")

    if(genome.v %in% c("hg19","GRCh37")){
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene, columns=c("gene_id"))
    }else if(genome.v %in% c("hg38","GRCh38")){
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        genesgr = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns=c("gene_id"))
    }else{stop("Unspecified, or non-available genome version; use hg19 or hg38")}
    
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
    
    # We will create genomic ranges objects with left and right junction sides 
    
    sv_df  <-  svdat[intersect(which(svdat$pos2 - svdat$pos1 < sv.seg.size), which(svdat$svclass %in% c("DUP","DEL","BND","INV","INS"))),]
    segmentEventGR = with(sv_df, GRanges(chrom1, IRanges(start=pos1, end=pos2),strand1=strand1,strand2=strand2,rowid = rownames(sv_df),svclass=svclass))
    non_seg_df <- svdat[setdiff(rownames(svdat),rownames(sv_df)),]
    leftJunctionGR = with(non_seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos1),Rle(strand1),rowid = rownames(non_seg_df),svclass=svclass))
    rightJunctionGR = with(non_seg_df, GRanges(chrom2, IRanges(start=pos2, end=pos2),Rle(strand2),rowid = rownames(non_seg_df),svclass=svclass))


    ###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs
    
    if(verbose) message("Finding disrupting overlaps with segmental junctions")
    segOverlaps = GenomicAlignments::findOverlaps(genesgr, segmentEventGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    segResults <- as.data.frame(segOverlaps)
    segHits <- sort(table(queryHits(segOverlaps)))[which( sort(table(queryHits(segOverlaps))) > 0)]
    segSamples <- sapply(names(segHits), function(i) svdat[segResults$subjectHits[which(segResults$queryHits == i)],"sample"],simplify=FALSE)
    segJunctions <- sapply(names(segHits), function(i) rownames(svdat)[segResults$subjectHits[which(segResults$queryHits == i)]],simplify=FALSE)
    names(segJunctions) <- names(segSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(segJunctions))]))
    
    if(verbose) message("Finding disrupting overlaps with left and right junction from translocation and lare SVs")
    leftOverlaps = GenomicAlignments::findOverlaps(genesgr,leftJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    leftResults <- as.data.frame(leftOverlaps)
    leftHits <- sort(table(queryHits(leftOverlaps)))[which( sort(table(queryHits(leftOverlaps))) > 0)]
    leftSamples <- sapply(names(leftHits), function(i) svdat[leftResults$subjectHits[which(leftResults$queryHits == i)],"sample"],simplify=FALSE)
    leftJunctions <- sapply(names(leftHits),function(i) rownames(svdat)[leftResults$subjectHits[which(leftResults$queryHits == i)]],simplify=FALSE)
    names(leftJunctions) <- names(leftSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(leftJunctions))]))

    rightOverlaps = GenomicAlignments::findOverlaps(genesgr,rightJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    rightResults <- as.data.frame(rightOverlaps)
    rightHits <- sort(table(queryHits(rightOverlaps)))[which( sort(table(queryHits(rightOverlaps))) > 0)]
    rightSamples <- sapply(names(rightHits),function(i) svdat[rightResults$subjectHits[which(rightResults$queryHits == i)],"sample"],simplify=FALSE)
    rightJunctions <- sapply(names(rightHits),function(i) rownames(svdat)[rightResults$subjectHits[which(rightResults$queryHits == i)]],simplify=FALSE)
    names(rightJunctions) <- names(rightSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(rightJunctions))]))
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
    
    if(verbose) message("Finding upstream gene overlaps with segmental junctions")
    segOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr, segmentEventGR,ignore.strand=TRUE,type="any")
    segResultsUp <- as.data.frame(segOverlapsUp)
    segHitsUp <- sort(table(queryHits(segOverlapsUp)))[which( sort(table(queryHits(segOverlapsUp))) > 0)]
    segSamplesUp <- sapply(names(segHitsUp), function(i) svdat[segResultsUp$subjectHits[which(segResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    segJunctionsUp <- sapply(names(segHitsUp), function(i) rownames(svdat)[segResultsUp$subjectHits[which(segResultsUp$queryHits == i)]],simplify=FALSE)
    names(segJunctionsUp) <- names(segSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(segJunctionsUp))]))
    
    if(verbose) message("Finding  upstream gene overlaps with left and right junction from translocation and lare SVs")
    leftOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,leftJunctionGR,ignore.strand=TRUE,type="any")
    leftResultsUp <- as.data.frame(leftOverlapsUp)
    leftHitsUp <- sort(table(queryHits(leftOverlapsUp)))[which( sort(table(queryHits(leftOverlapsUp))) > 0)]
    leftSamplesUp <- sapply(names(leftHitsUp), function(i) svdat[leftResultsUp$subjectHits[which(leftResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    leftJunctionsUp <- sapply(names(leftHitsUp),function(i) rownames(svdat)[leftResultsUp$subjectHits[which(leftResultsUp$queryHits == i)]],simplify=FALSE)
    names(leftJunctionsUp) <- names(leftSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(leftJunctionsUp))]))

    if(verbose) message("Finding  upstream gene overlaps with left and right junction from translocation and lare SVs")
    rightOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,rightJunctionGR,ignore.strand=TRUE,type="any")
    rightResultsUp <- as.data.frame(rightOverlapsUp)
    rightHitsUp <- sort(table(queryHits(rightOverlapsUp)))[which( sort(table(queryHits(rightOverlapsUp))) > 0)]
    rightSamplesUp <- sapply(names(rightHitsUp), function(i) svdat[rightResultsUp$subjectHits[which(rightResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    rightJunctionsUp <- sapply(names(rightHitsUp), function(i) rownames(svdat)[rightResultsUp$subjectHits[which(rightResultsUp$queryHits == i)]],simplify=FALSE)
    names(rightJunctionsUp) <- names(rightSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(rightJunctionsUp))]))
    
    if(verbose) message("Wrapping up!")
    disruptSamples <- merge2lists(merge2lists(leftSamples,rightSamples),segSamples)
    disruptJunctions <- merge2lists(merge2lists(leftJunctions,rightJunctions),segJunctions)

    upstreamSamples <- merge2lists(merge2lists(leftSamplesUp,rightSamplesUp),segSamplesUp)
    upstreamJunctions <- merge2lists(merge2lists(leftJunctionsUp,rightJunctionsUp),segJunctionsUp)
    
    return(break.annot(
    input = svdat,
    geneinfo = genes_tab,
    disruptSamples = disruptSamples,
    disruptBreaks = disruptJunctions,
    upstreamSamples = upstreamSamples,
    upstreamBreaks = upstreamJunctions,
    param = list(
        genome.v = genome.v,
        maxgap = maxgap, 
        upstr = upstr, 
        sv.seg.size = sv.seg.size, 
        verbose = verbose)
    ))
    
}


#' Identify recurrently altered genes
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param maxgap (numeric) offset region size (base pairs) beyond transcription start and end sites to extend for overlaps
#' @param upstr (numeric) upstream region size (base pairs) to identify breakpoints in gene upstream regulatory regions 
#' @param min.seg.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param clean.brk (numeric) identical CNV segments across multiple samples are removed and replaced by neighbours average   
#' @param verbose (logical) 
#' @return an instance of the class 'break.annot' containing breakpoint mapping onto genes
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' cnv.break.annot(seg)

cnv.break.annot <- function(seg, 
                            fc.pct = 0.2, 
                            genome.v="hg19",
                            maxgap= 1000,
                            upstr=50000,
                            min.seg.size = NULL, 
                            min.num.probes = NULL, 
                            low.cov = NULL,
                            clean.brk = NULL,
                            verbose=TRUE){
    
    segdat <- validate.seg(seg)
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
                             low.cov = low.cov, clean.brk = clean.brk, verbose=verbose)
    
    breaks<- cnv_breaks$breaks
    rownames(breaks) <- createRandomString(nrow(breaks),8)
    
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
    
    return(break.annot(
        input= breaks,
        geneinfo = genes_tab,
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




