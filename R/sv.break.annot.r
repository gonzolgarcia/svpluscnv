#' 
#'
#' Identify recurrently altered genes by strutural variants
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param maxgap (numeric) offset region size (base pairs) beyond transcription start and end sites to extend for overlaps
#' @param upstr (numeric) upstream region size (base pairs) to identify breakpoints in gene upstream regulatory regions 
#' @param sv.seg.size (numeric) segmental variants (DEL, DUP, INV or INS) size limit; bigger sigments are considered translocations
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
    
    rownames(svdat) <- createRandomString(nrow(svdat),10)
    
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
    
    seg_df  <-  svdat[intersect(which(svdat$pos2 - svdat$pos1 < sv.seg.size), which(svdat$svclass %in% c("DUP","DEL","BND","INV","INS"))),]
    segmentEventGR = with(seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos2),strand1=strand1,strand2=strand2,rowid = rownames(seg_df),svclass=svclass))
    non_seg_df <- svdat[setdiff(rownames(svdat),rownames(seg_df)),]
    leftJunctionGR = with(non_seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos1),Rle(strand1),rowid = rownames(non_seg_df),svclass=svclass))
    rightJunctionGR = with(non_seg_df, GRanges(chrom2, IRanges(start=pos2, end=pos2),Rle(strand2),rowid = rownames(non_seg_df),svclass=svclass))


    ###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs
    
    # find disrupting overlaps with segmental junctions
    segOverlaps = GenomicAlignments::findOverlaps(genesgr, segmentEventGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    segResults <- as.data.frame(segOverlaps)
    segHits <- sort(table(queryHits(segOverlaps)))[which( sort(table(queryHits(segOverlaps))) > 0)]
    segSamples <- sapply(names(segHits), function(i) svdat[segResults$subjectHits[which(segResults$queryHits == i)],"sample"],simplify=FALSE)
    segJunctions <- sapply(names(segHits), function(i) rownames(svdat)[segResults$subjectHits[which(segResults$queryHits == i)]],simplify=FALSE)
    names(segJunctions) <- names(segSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(segJunctions))]))
    
    # find disrupting overlaps with left and right junction from translocation and lare SVs
    leftOverlaps = GenomicAlignments::findOverlaps(genesgr,leftJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    leftResults <- as.data.frame(leftOverlaps)
    leftHits <- sort(table(queryHits(leftOverlaps)))[which( sort(table(queryHits(leftOverlaps))) > 0)]
    leftSamples <- sapply(names(leftHits), function(i) svdat[leftResults$subjectHits[which(leftResults$queryHits == i)],"sample"],simplify=FALSE)
    leftJunctions <- sapply(names(leftHits),function(i) rownames(svdat)[leftResults$subjectHits[which(leftResults$queryHits == i)]],simplify=FALSE)
    names(leftJunctions) <- names(leftSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(leftJunctions))]))
    # sort(unlist(lapply(leftSamples,function(x) length(unique(x)))),decreasing=T)[1:100]

    rightOverlaps = GenomicAlignments::findOverlaps(genesgr,rightJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    rightResults <- as.data.frame(rightOverlaps)
    rightHits <- sort(table(queryHits(rightOverlaps)))[which( sort(table(queryHits(rightOverlaps))) > 0)]
    rightSamples <- sapply(names(rightHits),function(i) svdat[rightResults$subjectHits[which(rightResults$queryHits == i)],"sample"],simplify=FALSE)
    rightJunctions <- sapply(names(rightHits),function(i) rownames(svdat)[rightResults$subjectHits[which(rightResults$queryHits == i)]],simplify=FALSE)
    names(rightJunctions) <- names(rightSamples) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(rightJunctions))]))
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
    
    # find disrupting overlaps with segmental junctions
    segOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr, segmentEventGR,ignore.strand=TRUE,type="any")
    segResultsUp <- as.data.frame(segOverlapsUp)
    segHitsUp <- sort(table(queryHits(segOverlapsUp)))[which( sort(table(queryHits(segOverlapsUp))) > 0)]
    segSamplesUp <- sapply(names(segHitsUp), function(i) svdat[segResultsUp$subjectHits[which(segResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    segJunctionsUp <- sapply(names(segHitsUp), function(i) rownames(svdat)[segResultsUp$subjectHits[which(segResultsUp$queryHits == i)]],simplify=FALSE)
    names(segJunctionsUp) <- names(segSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(segJunctionsUp))]))
    
    # find disrupting overlaps with left and right junction from translocation and lare SVs
    leftOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,leftJunctionGR,ignore.strand=TRUE,type="any")
    leftResultsUp <- as.data.frame(leftOverlapsUp)
    leftHitsUp <- sort(table(queryHits(leftOverlapsUp)))[which( sort(table(queryHits(leftOverlapsUp))) > 0)]
    leftSamplesUp <- sapply(names(leftHitsUp), function(i) svdat[leftResultsUp$subjectHits[which(leftResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    leftJunctionsUp <- sapply(names(leftHitsUp),function(i) rownames(svdat)[leftResultsUp$subjectHits[which(leftResultsUp$queryHits == i)]],simplify=FALSE)
    names(leftJunctionsUp) <- names(leftSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(leftJunctionsUp))]))
    # sort(unlist(lapply(leftSamples,function(x) length(unique(x)))),decreasing=T)[1:100]
    
    # find disrupting overlaps with left and right junction from translocation and lare SVs
    rightOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,rightJunctionGR,ignore.strand=TRUE,type="any")
    rightResultsUp <- as.data.frame(rightOverlapsUp)
    rightHitsUp <- sort(table(queryHits(rightOverlapsUp)))[which( sort(table(queryHits(rightOverlapsUp))) > 0)]
    rightSamplesUp <- sapply(names(rightHitsUp), function(i) svdat[rightResultsUp$subjectHits[which(rightResultsUp$queryHits == i)],"sample"],simplify=FALSE)
    rightJunctionsUp <- sapply(names(rightHitsUp), function(i) rownames(svdat)[rightResultsUp$subjectHits[which(rightResultsUp$queryHits == i)]],simplify=FALSE)
    names(rightJunctionsUp) <- names(rightSamplesUp) <- unname(unlist(genes_tab$SYMBOL[as.numeric(names(rightJunctionsUp))]))
    
    # combine left anbd right hits 
    disruptSamples <- merge2lists(merge2lists(leftSamples,rightSamples),segSamples)
    disruptJunctions <- merge2lists(merge2lists(leftJunctions,rightJunctions),segJunctions)

    upstreamSamples <- merge2lists(merge2lists(leftSamplesUp,rightSamplesUp),segSamplesUp)
    upstreamJunctions <- merge2lists(merge2lists(leftJunctionsUp,rightJunctionsUp),segJunctionsUp)
    
    return(list(
    svdat= svdat,
    genes = genes_tab,
    disruptSamples = disruptSamples,
    disruptJunctions = disruptJunctions,
    upstreamSamples = upstreamSamples,
    upstreamJunctions = upstreamJunctions,
    param = list(
        genome.v = genome.v,
        maxgap = maxgap, 
        upstr = upstr, 
        sv.seg.size = sv.seg.size, 
        verbose = verbose)
    ))
    
}

