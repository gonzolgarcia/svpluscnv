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
    writeLines(paste("An object of class break.annot from svcnvplus containing the following stats:
                \nNumber of samples=",length(unique(object@input$sample)),
                "\nAltered genes=",length(object@disruptSamples)))
})



#' Identify recurrently altered genes by strutural variants
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param genome.v (hg19 or hg38) reference genome version to retrieve gene annotations 
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
                     genetab=NULL,
                     maxgap=0,
                     upstr=50000, 
                     sv.seg.size=200000,
                     verbose=TRUE){
    
    svdat <- validate.sv(sv)

    chrlist <- unique(c(svdat$chrom1,svdat$chrom2))
    chrlist <- chr.sort(chrlist)
    
    rownames(svdat) <- createRandomString(nrow(svdat),8)
    
    if(verbose) message("# Generating GRanges objects")

    genesgr <- get.genesgr(genome.v=genome.v)
    
    genes_df <- as.data.frame(genesgr)
    rownames(genes_df) <- genes_df$gene_id
    genes_ups <- genes_df
    genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - maxgap
    genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr - maxgap
    genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + maxgap
    genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr + maxgap
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))

    
    # We will create genomic ranges objects with left and right junction sides 
    
    seg_df  <-  svdat[intersect(which(svdat$pos2 - svdat$pos1 < sv.seg.size), which(svdat$svclass %in% c("DUP","DEL","BND","INV","INS"))),]
    segmentEventGR = with(seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos2),strand1=strand1,strand2=strand2,rowid = rownames(seg_df),svclass=svclass))
    non_seg_df <- svdat[setdiff(rownames(svdat),rownames(seg_df)),]
    leftJunctionGR = with(non_seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos1),Rle(strand1),rowid = rownames(non_seg_df),svclass=svclass))
    rightJunctionGR = with(non_seg_df, GRanges(chrom2, IRanges(start=pos2, end=pos2),Rle(strand2),rowid = rownames(non_seg_df),svclass=svclass))

    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs
    
    if(verbose) message("Finding disrupting overlaps with segmental junctions")
    segOverlaps = GenomicAlignments::findOverlaps(genesgr, segmentEventGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(segOverlaps)]
    svHits <- segmentEventGR@elementMetadata$rowid[subjectHits(segOverlaps)]
    segResults <- remove.factors(data.frame(geneHits,svHits))
    aggList <- aggregate(svHits~geneHits,segResults,unique,simplify=FALSE)
    segJunctions <- aggList$svHits
    segSamples <- lapply(segJunctions,function(x) unique(svdat[x,"sample"]))
    names(segSamples) <- names(segJunctions) <- aggList$geneHits

    if(verbose) message("Finding disrupting overlaps with left and right junction from translocation and lare SVs")
    leftOverlaps = GenomicAlignments::findOverlaps(genesgr,leftJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(leftOverlaps)]
    svHits <- leftJunctionGR@elementMetadata$rowid[subjectHits(leftOverlaps)]
    leftResults <- remove.factors(data.frame(geneHits,svHits))
    leftJunctions <- sapply(unique(leftResults$geneHits), function(i) leftResults$svHits[which(leftResults$geneHits == i)],simplify=FALSE)
    leftSamples <- lapply(leftJunctions,function(x) unique(svdat[x,"sample"]))

    rightOverlaps = GenomicAlignments::findOverlaps(genesgr,rightJunctionGR,ignore.strand=TRUE,type="any",maxgap=maxgap)
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(rightOverlaps)]
    svHits <- rightJunctionGR@elementMetadata$rowid[subjectHits(rightOverlaps)]
    rightResults <- remove.factors(data.frame(geneHits,svHits))
    rightJunctions <- sapply(unique(rightResults$geneHits), function(i) rightResults$svHits[which(rightResults$geneHits == i)],simplify=FALSE)
    rightSamples <- lapply(rightJunctions,function(x) unique(svdat[x,"sample"]))
    

    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
    
    if(verbose) message("Finding upstream gene overlaps with segmental junctions")
    segOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr, segmentEventGR,ignore.strand=TRUE,type="any")
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(segOverlapsUp)]
    svHits <- segmentEventGR@elementMetadata$rowid[subjectHits(segOverlapsUp)]
    segResultsUp <- remove.factors(data.frame(geneHits,svHits))
    segJunctionsUp <- sapply(unique(segResultsUp$geneHits), function(i) segResultsUp$svHits[which(segResultsUp$geneHits == i)],simplify=FALSE)
    segSamplesUp <- lapply(segJunctionsUp,function(x) unique(svdat[x,"sample"]))
    

    if(verbose) message("Finding  upstream gene overlaps with left and right junction from translocation and lare SVs")
    leftOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,leftJunctionGR,ignore.strand=TRUE,type="any")
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(leftOverlapsUp)]
    svHits <- leftJunctionGR@elementMetadata$rowid[subjectHits(leftOverlapsUp)]
    leftResultsUp <- remove.factors(data.frame(geneHits,svHits))
    leftJunctionsUp <- sapply(unique(leftResultsUp$geneHits), function(i) leftResultsUp$svHits[which(leftResultsUp$geneHits == i)],simplify=FALSE)
    leftSamplesUp <- lapply(leftJunctionsUp,function(x) unique(svdat[x,"sample"]))

    rightOverlapsUp = GenomicAlignments::findOverlaps(upstreamgr,rightJunctionGR,ignore.strand=TRUE,type="any")
    geneHits <- genesgr@elementMetadata$gene_id[queryHits(rightOverlapsUp)]
    svHits <- rightJunctionGR@elementMetadata$rowid[subjectHits(rightOverlapsUp)]
    rightResultsUp <- remove.factors(data.frame(geneHits,svHits))
    rightJunctionsUp <- sapply(unique(rightResultsUp$geneHits), function(i) rightResultsUp$svHits[which(rightResultsUp$geneHits == i)],simplify=FALSE)
    rightSamplesUp <- lapply(rightJunctionsUp,function(x) unique(svdat[x,"sample"]))
    
    if(verbose) message("Wrapping up!")
    disruptSamples <- merge2lists(merge2lists(leftSamples,rightSamples,fun = "unique"),segSamples, fun="unique")
    disruptJunctions <- merge2lists(merge2lists(leftJunctions,rightJunctions,fun = "unique"),segJunctions,fun = "unique")

    upstreamSamples <- merge2lists(merge2lists(leftSamplesUp,rightSamplesUp,fun = "unique"),segSamplesUp,fun = "unique")
    upstreamJunctions <- merge2lists(merge2lists(leftJunctionsUp,rightJunctionsUp,fun = "unique"),segJunctionsUp,fun = "unique")
    
    return(break.annot(
    input = svdat,
    geneinfo = genes_df,
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
                            genetab=NULL,
                            maxgap= 1000,
                            upstr=50000,
                            min.seg.size = NULL, 
                            min.num.probes = NULL, 
                            low.cov = NULL,
                            clean.brk = NULL,
                            verbose=TRUE){
    
    segdat <- validate.seg(seg)
    chrlist <- unique(segdat$chrom)
    chrlist <- chr.sort(chrlist)
    
    if(verbose) message("# Generating GRanges objects")
    
    genesgr <- get.genesgr(genome.v=genome.v)
    
    genes_df <- as.data.frame(genesgr)
    rownames(genes_df) <- genes_df$gene_id
    genes_ups <- genes_df
    genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - maxgap
    genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr - maxgap
    genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + maxgap
    genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr + maxgap
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    
    
    cnv_breaks <- seg.breaks(segdat, fc.pct = fc.pct, min.seg.size = min.seg.size, min.num.probes = min.num.probes, 
                             low.cov = low.cov, clean.brk = clean.brk, verbose=verbose)
    
    breaks<- cnv_breaks$breaks
    rownames(breaks) <- createRandomString(nrow(breaks),8)
    
    breaksgr = with(breaks[,c("chrom","start","end")], GRanges(chrom, IRanges(start=start, end=end)))
    
    overlapGenes <- GenomicAlignments::findOverlaps(genesgr,breaksgr,ignore.strand=TRUE,maxgap = maxgap)
    
    overlapUpstream <- GenomicAlignments::findOverlaps(upstreamgr,breaksgr,ignore.strand=TRUE)

    if(verbose) message("Finding gene disrupting CNV breakpoint")
    break_id <- rownames(breaks)[subjectHits(overlapGenes)]
    gene_id <- genes_df[queryHits(overlapGenes),"gene_id"]
    geneHits <- remove.factors(data.frame(break_id,gene_id))
    
    aggList <- aggregate(break_id~gene_id,geneHits,unique,simplify=FALSE)
    geneBreaks <- aggList$break_id
    geneSamples <- lapply(geneBreaks,function(x) unique(breaks[x,"sample"]))
    names(geneSamples) <- names(geneBreaks) <- aggList$gene_id
    
    if(verbose) message("Finding upstream gene disrupting CNV breakpoint")
    break_id <- rownames(breaks)[subjectHits(overlapUpstream)]
    gene_id <- genes_df[queryHits(overlapUpstream),"gene_id"]
    upstreamHits <- remove.factors(data.frame(break_id,gene_id))

    aggList <- aggregate(break_id~gene_id,upstreamHits,unique,simplify=FALSE)
    upstreamBreaks <- aggList$break_id
    upstreamSamples <- lapply(upstreamBreaks,function(x) unique(breaks[x,"sample"]))
    names(upstreamSamples) <- names(upstreamBreaks) <- aggList$gene_id
    
     
    return(break.annot(
        input= breaks,
        geneinfo = genes_df,
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




