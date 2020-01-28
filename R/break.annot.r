#' Class to store breakpoint annotations
#' @param input (data.frame): chromosome for the first breakpoint
#' @param genesgr (GRanges): a GRanges object with genomic segments (e.g. genes) 
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
                       genesgr = "GRanges",
                       disruptSamples = "list",
                       disruptBreaks = "list",
                       upstreamSamples = "list",
                       upstreamBreaks = "list",
                       dnstreamSamples = "list",
                       dnstreamBreaks = "list",
                       param = "list"
                   ))


setMethod("show","break.annot",function(object){
    writeLines(paste("An object of class break.annot from svcnvplus containing the following stats:
                \nNumber of samples=",length(unique(object@input$sample)),
                "\nDirrupted genes=",length(object@disruptSamples),
                "\nUpstream disrupted genes=",length(object@upstreamSamples),
                "\nDownstream disrupted genes=",length(object@dnstreamSamples)))
})



#' Identify recurrently altered genes by strutural variants
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param genome.v (character): either 'hg19' or 'hg38' accepted; reference genome version to retrieve gene annotations 
#' @param genesgr (S4) a GenomicRanges object containing gene annotations (if not NULL overides genome.v). It must containg 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...) 
#' @param upstr (numeric) upstream region size (base pairs) to identify breakpoints in gene upstream  regions 
#' @param dnstr (numeric) upstream region size (base pairs) to identify breakpoints in gene downstream  regions 
#' @param sv.seg.size (numeric) segmental variants (DEL, DUP, INV or INS) size limit; bigger sigments are considered translocations, in that case, only the breakpoint 
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
                     genesgr=NULL,
                     upstr=50000,
                     dnstr=50000,
                     sv.seg.size=200000,
                     verbose=TRUE){
    
svdat <- validate.sv(sv)

chrlist <- unique(c(svdat$chrom1,svdat$chrom2))
chrlist <- chr.sort(chrlist)
    
rownames(svdat) <- createRandomString(nrow(svdat),8)

if(!is.null(genesgr)){
    if(anyDuplicated(genesgr@elementMetadata$gene_id > 0)) stop("The genesgr provided object contains duplicated gene_id values")
}else{
    genesgr <- get.genesgr(genome.v=genome.v)
}

genes_df <- as.data.frame(genesgr)
rownames(genes_df) <- genes_df$gene_id

# We will create genomic ranges objects with left and right junction sides 
if(verbose) message(paste("Generating GRanges objects for SV breaks and SV segments < ",sv.seg.size,"bp",sep=""))
seg_df  <-  svdat[intersect(which(svdat$pos2 - svdat$pos1 < sv.seg.size), which(svdat$svclass %in% c("DUP","DEL","BND","INV","INS"))),]
segmentEventGR = with(seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos2),rowid = rownames(seg_df)))
non_seg_df <- svdat[setdiff(rownames(svdat),rownames(seg_df)),]
leftJunctionGR = with(non_seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos1),rowid = rownames(non_seg_df)))
rightJunctionGR = with(non_seg_df, GRanges(chrom2, IRanges(start=pos2, end=pos2),rowid = rownames(non_seg_df)))
svRanges <- c(segmentEventGR,leftJunctionGR,rightJunctionGR)
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs

if(verbose) message("Finding disrupting breakpoint overlaps")
segOverlaps = GenomicAlignments::findOverlaps(genesgr, svRanges,ignore.strand=TRUE,type="any")
gene_id <- genesgr@elementMetadata$gene_id[queryHits(segOverlaps)]
sv_id <- svRanges@elementMetadata$rowid[subjectHits(segOverlaps)]
sample_id <- svdat[sv_id,"sample"]
svResults <- data.table(gene_id,sv_id,sample_id)
ans <- svResults[, .(sv_id = list(sv_id), sample_id=list(unique(sample_id))), by = "gene_id"]
disruptBreaks <- ans$sv_id
disruptSamples <- ans$sample_id
names(disruptBreaks) <- names(disruptSamples) <- ans$gene_id

    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs

if(!is.null(upstr)){
    if(verbose) message("Finding upstream breakpoint overlaps")

    genes_ups <- genes_df
    genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"]  
    genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr 
    genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"] 
    genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    
    
    overlapsUp = GenomicAlignments::findOverlaps(upstreamgr, svRanges,ignore.strand=TRUE,type="any")
    gene_id <- genesgr@elementMetadata$gene_id[queryHits(overlapsUp)]
    sv_id <- svRanges@elementMetadata$rowid[subjectHits(overlapsUp)]
    sample_id <- svdat[sv_id,"sample"]
    svResultsUp <- data.table(gene_id,sv_id,sample_id)
    ans <- svResultsUp[, .(sv_id = list(sv_id), sample_id=list(unique(sample_id))), by = "gene_id"]
    upstreamBreaks <- ans$sv_id
    upstreamSamples <- ans$sample_id
    names(upstreamBreaks) <- names(upstreamSamples) <- ans$gene_id
}else{
    upstreamBreaks <- upstreamSamples <- list()
}

    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
if(!is.null(dnstr)){
    
    if(verbose) message("Finding upstream breakpoint overlaps")
    
    genes_dns <- genes_df
    genes_dns[which(genes_dns$strand == "+"),"end"] <- genes_dns[which(genes_dns$strand == "+"),"end"]  + dnstr
    genes_dns[which(genes_dns$strand == "+"),"start"] <- genes_dns[which(genes_dns$strand == "+"),"end"] 
    genes_dns[which(genes_dns$strand == "-"),"start"] <- genes_dns[which(genes_dns$strand == "-"),"start"]  - dnstr
    genes_dns[which(genes_dns$strand == "-"),"end"] <- genes_dns[which(genes_dns$strand == "-"),"start"] 
    dnstreamgr <- with(genes_dns, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    
    overlapsDn = GenomicAlignments::findOverlaps(dnstreamgr, svRanges,ignore.strand=TRUE,type="any")
    gene_id <- genesgr@elementMetadata$gene_id[queryHits(overlapsDn)]
    sv_id <- svRanges@elementMetadata$rowid[subjectHits(overlapsDn)]
    sample_id <- seg_df[sv_id,"sample"]
    svResultsDn <- data.table(gene_id,sv_id,sample_id)
    ans <- svResultsDn[, .(sv_id = list(sv_id), sample_id=list(unique(sample_id))), by = "gene_id"]
    dnstreamBreaks <- ans$sv_id
    dnstreamSamples <- ans$sample_id
    names(dnstreamBreaks) <- names(dnstreamSamples) <- ans$gene_id
}else{
    dnstreamBreaks <- dnstreamSamples <- list()
}

    return(break.annot(
    input = svdat,
    genesgr = genesgr,
    disruptSamples = disruptSamples,
    disruptBreaks = disruptBreaks,
    upstreamSamples = upstreamSamples,
    upstreamBreaks = upstreamBreaks,
    dnstreamSamples = dnstreamSamples,
    dnstreamBreaks = dnstreamBreaks,
    param = list(
        genome.v = genome.v,
        upstr = upstr, 
        sv.seg.size = sv.seg.size, 
        verbose = verbose)
    ))
    
}


#' Identify recurrently altered genes
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param genegr (GRanges) a GRanges object containing gene annotations (if not NULL overides genom.v). It must containg 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...) 
#' @param upstr (numeric) genomic region size (base pairs) to identify breakpoints overlapping with gene upstream  regions 
#' @param dnstr (numeric) genomic region size (base pairs) to identify breakpoints overlapping with gene downstream  regions 
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
                            genesgr=NULL,
                            upstr=50000,
                            dnstr=50000,
                            min.seg.size = NULL, 
                            min.num.probes = NULL, 
                            low.cov = NULL,
                            clean.brk = NULL,
                            verbose=TRUE){
    
segdat <- validate.seg(seg)
chrlist <- unique(segdat$chrom)
chrlist <- chr.sort(chrlist)
    

if(!is.null(genesgr)){
    if(anyDuplicated(genesgr@elementMetadata$gene_id > 0)) stop("The genesgr provided object contains duplicated gene_id values")
}else{
    genesgr <- get.genesgr(genome.v=genome.v)
}
    
genes_df <- as.data.frame(genesgr)
rownames(genes_df) <- genes_df$gene_id
    
    
    
cnv_breaks <- seg.breaks(segdat, fc.pct = fc.pct, min.seg.size = min.seg.size, min.num.probes = min.num.probes, 
                         low.cov = low.cov, clean.brk = clean.brk, verbose=verbose)
    
breaks<- cnv_breaks$breaks
rownames(breaks) <- createRandomString(nrow(breaks),8)
    
breaksgr = with(breaks[,c("chrom","start","end")], GRanges(chrom, IRanges(start=start, end=end)))

    

if(verbose) message("Finding gene disrupting CNV breakpoint")
overlapGenes <- GenomicAlignments::findOverlaps(genesgr,breaksgr,ignore.strand=TRUE)
break_id <- rownames(breaks)[subjectHits(overlapGenes)]
gene_id <- genes_df[queryHits(overlapGenes),"gene_id"]

sample_id <- breaks[break_id,"sample"]
brResults <- data.table(gene_id,break_id,sample_id)
ans <- brResults[, .(break_id = list(break_id), sample_id=list(unique(sample_id))), by = "gene_id"]
disruptBreaks <- ans$break_id
disruptSamples <- ans$sample_id
names(disruptBreaks) <- names(disruptSamples) <- ans$gene_id

    
if(!is.null(upstr)){
    if(verbose) message("Finding upstream CNV breakpoint overlaps")
    genes_ups <- genes_df
    genes_ups[which(genes_ups$strand == "+"),"end"] <- genes_ups[which(genes_ups$strand == "+"),"start"]
    genes_ups[which(genes_ups$strand == "+"),"start"] <- genes_ups[which(genes_ups$strand == "+"),"start"] - upstr
    genes_ups[which(genes_ups$strand == "-"),"start"] <- genes_ups[which(genes_ups$strand == "-"),"end"]
    genes_ups[which(genes_ups$strand == "-"),"end"] <- genes_ups[which(genes_ups$strand == "-"),"end"] + upstr
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    overlapUpstream <- GenomicAlignments::findOverlaps(upstreamgr,breaksgr,ignore.strand=TRUE)
    
    break_id <- rownames(breaks)[subjectHits(overlapUpstream)]
    gene_id <- genes_df[queryHits(overlapUpstream),"gene_id"]
    sample_id <- breaks[break_id,"sample"]
    brResultsUp <- data.table(gene_id,break_id,sample_id)
    ans <- brResultsUp[, .(break_id = list(break_id), sample_id=list(unique(sample_id))), by = "gene_id"]
    upstreamBreaks <- ans$break_id
    upstreamSamples <- ans$sample_id
    names(upstreamBreaks) <- names(upstreamSamples) <- ans$gene_id
}else{
    upstreamBreaks <- upstreamSamples <- list()
}

if(!is.null(dnstr)){
    if(verbose) message("Finding downstream CNV breakpoint overlaps")
    genes_dns <- genes_df
    genes_dns[which(genes_dns$strand == "+"),"end"] <- genes_dns[which(genes_dns$strand == "+"),"end"]  + dnstr
    genes_dns[which(genes_dns$strand == "+"),"start"] <- genes_dns[which(genes_dns$strand == "+"),"end"] 
    genes_dns[which(genes_dns$strand == "-"),"start"] <- genes_dns[which(genes_dns$strand == "-"),"start"] - dnstr
    genes_dns[which(genes_dns$strand == "-"),"end"] <- genes_dns[which(genes_dns$strand == "-"),"start"] 
    dnstreamgr <- with(genes_dns, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    overlapDnstream <- GenomicAlignments::findOverlaps(dnstreamgr,breaksgr,ignore.strand=TRUE)
    
    break_id <- rownames(breaks)[subjectHits(overlapDnstream)]
    gene_id <- genes_df[queryHits(overlapDnstream),"gene_id"]
    sample_id <- breaks[break_id,"sample"]
    brResultsDn <- data.table(gene_id,break_id,sample_id)
    ans <- brResultsDn[, .(break_id = list(break_id), sample_id=list(unique(sample_id))), by = "gene_id"]
    dnstreamBreaks <- ans$break_id
    dnstreamSamples <- ans$sample_id
    names(dnstreamBreaks) <- names(dnstreamSamples) <- ans$gene_id
}else{
    dnstreamBreaks <- dnstreamSamples <- list()
}

return(break.annot(
    input= breaks,
    genesgr = genesgr,
    disruptSamples = disruptSamples,
    disruptBreaks = disruptBreaks,
    upstreamSamples = upstreamSamples,
    upstreamBreaks = upstreamBreaks,
    dnstreamSamples = dnstreamSamples,
    dnstreamBreaks = dnstreamBreaks,
    param = list(
        fc.pct=fc.pct,
        genome.v = genome.v,
        upstr = upstr, 
        dnstr = dnstr,
        min.seg.size=min.seg.size,
        min.num.probes=min.num.probes,
        low.cov=low.cov,
        clean.brk=clean.brk)
    ))
    
}

