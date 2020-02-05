#' Class to store breakpoint annotations in association with genomic features (e.g. gene loci)
#' @param input (data.frame): the breakpoint info containing data.frame, this will be occupied by the CNV segmentation data in the case of cnv.break.annot or SV for sv.break.annot. Unique random string rownames are added to the provided data.frame.
#' @param genesgr (GRanges): a GRanges object with genomic features (e.g. genes) to which breakpoints are mapped 
#' @param disruptSamples (list): a list which names correspond to genomic features and values correspond to sample ids harboring breakpoints overlapping with said features
#' @param disruptBreaks (list): a list which names correspond to genomic features and values correspond to the ids of breakpount mapped onto them. Break ids are linked to the 'input' data.frame rownames 
#' @param upstreamSamples (list): a list which names correspond to genomic features and values correspond to sample ids harboring breakpoints overlapping with upstream region of said features
#' @param upstreamBreaks (list): a list which names correspond to genomic features and values correspond to the ids of breakpount mapped onto upstream regions Break ids are linked to the 'input' data.frame rownames 
#' @param param (list): a list of parametres provided for the annotation function
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
    writeLines(paste("break.annot class object from svcnvplus containing the following stats:
                \nNumber of samples=",length(unique(object@input$sample)),
                "\nDirrupted genes=",length(object@disruptSamples),
                "\nUpstream disrupted genes=",length(object@upstreamSamples),
                "\nDownstream disrupted genes=",length(object@dnstreamSamples)))
})



#' Finf overlaps between genomic features and SV breakpoints; 
#' @param ggr (S4) a GenomicRanges object containing gene annotations. It is crutial that the genome version 'genesgr' and the input 'sv' are the same. The GRanges object must contain 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...). 
#' @param svgr (S4) a GenomicRanges object containing SV breakpoint ends. Metadata must contain 'rowid' and 'sampleid' fields. Seqnames are expected in the format (chr1, chr2, ...). Used by 'svc.break.annot' and 'cnv.break.annot'
#' @export
#' @keywords internal

geneBreakOverlap <- function(ggr,svgr){
    segOverlaps = GenomicAlignments::findOverlaps(ggr, svgr, ignore.strand=TRUE, type="any")
    gene_id <- ggr@elementMetadata$gene_id[queryHits(segOverlaps)]
    sv_id <- svgr@elementMetadata$rowid[subjectHits(segOverlaps)]
    sample_id <- svgr@elementMetadata$sampleid[subjectHits(segOverlaps)]
    svResults <- data.table(gene_id,sv_id,sample_id)
    ans <- svResults[, .(sv_id = list(sv_id), sample_id=list(unique(sample_id))), by = "gene_id"]
    geneBreaks <- ans$sv_id
    geneSamples <- ans$sample_id
    names(geneBreaks) <- names(geneSamples) <- ans$gene_id
    return(list(
        geneBreaks=geneBreaks,
        geneSamples=geneSamples
    ))
}

#' generate GRanges of upstream regions 
#' @param ggr (S4) a GenomicRanges object containing gene annotations. It is crutial that the genome version 'genesgr' and the input 'sv' are the same. The GRanges object must contain 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...). 
#' @param upstr (numeric) size in base pairs to define gene upstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of upstream regions.
#' @export


upgr <- function(ggr,upstr=50000){
    genes_df <- genes_ups <- as.data.frame(ggr)
    genes_ups[which(genes_df$strand == "+"),"end"] <- genes_df[which(genes_df$strand == "+"),"start"] - 1
    genes_ups[which(genes_df$strand == "+"),"start"] <- genes_df[which(genes_df$strand == "+"),"start"] - upstr
    genes_ups[which(genes_df$strand == "-"),"start"] <- genes_df[which(genes_df$strand == "-"),"end"] + 1
    genes_ups[which(genes_df$strand == "-"),"end"] <- genes_df[which(genes_df$strand == "-"),"end"] + upstr
    upstreamgr <- with(genes_ups, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    return(upstreamgr)
}

#' generate GRanges of downstream regions 
#' @param ggr (S4) a GenomicRanges object containing gene annotations. It is crutial that the genome version 'genesgr' and the input 'sv' are the same. The GRanges object must contain 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...). 
#' @param dnstr (numeric) size in base pairs to define gene downstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of downstream regions.
#' @export


dngr <- function(ggr,dnstr=50000){
    genes_df <- genes_dns <- as.data.frame(ggr)
    genes_dns[which(genes_df$strand == "+"),"start"] <- genes_df[which(genes_df$strand == "+"),"end"] +1
    genes_dns[which(genes_df$strand == "+"),"end"] <- genes_df[which(genes_df$strand == "+"),"end"]  + dnstr
    genes_dns[which(genes_df$strand == "-"),"start"] <- genes_df[which(genes_df$strand == "-"),"start"] - dnstr
    genes_dns[which(genes_df$strand == "-"),"end"] <- genes_df[which(genes_df$strand == "-"),"start"] -1
    dnstreamgr <- with(genes_dns, GRanges(seqnames, IRanges(start = start, end = end), strand = strand, gene_id=gene_id))
    return(dnstreamgr)
}

#' Identify recurrently altered genes by strutural variants. The function will identify overlaps between genomic features (e.g. genes) and SVs breakpoints. 
#' 
#' @param svc (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @param genome.v (character): either 'hg19' or 'hg38' accepted; reference genome version to retrieve gene annotations including genomic coordinates and strand
#' @param genesgr (S4) a GenomicRanges object containing gene annotations (if not NULL overides genome.v). It is crutial that the genome version 'genesgr' and the input 'sv' are the same. The GRanges object must contain 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...). 
#' @param upstr (numeric) size in base pairs to define gene upstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of upstream regions.
#' @param dnstr (numeric) size in base pairs to define gene downstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of downstream regions.
#' @param svc.seg.size (numeric) base pairs for maximum allowed segmental variants (DEL, DUP, INV or INS) size. Larger segmental SVs are treated as translocations and only the breakpoint position will be overlapped with genomic features.
#' @return an instance of the class 'break.annot' containing breakpoint mapping onto genes
#' @keywords Structural variants, annotation
#' @export
#' @examples
#' 
#' ## Obtain breakpoints from SV calls data
#' svc <- validate.svc(svcdat_lung_ccle)
#' 
#' svc.break.annot(svc, genome.v="hg19")

svc.break.annot <- function(svc, 
                            genome.v="hg19",
                            genesgr=NULL,
                            upstr=50000,
                            dnstr=50000,
                            svc.seg.size=200000,
                            verbose=TRUE){
    
    svcdat <- validate.svc(svc)
    
    chrlist <- unique(c(svcdat$chrom1,svcdat$chrom2))
    chrlist <- chr.sort(chrlist)
    
    rownames(svcdat) <- createRandomString(nrow(svcdat),8)
    
    if(!is.null(genesgr)){
        if(anyDuplicated(genesgr@elementMetadata$gene_id > 0)) stop("The genesgr provided object contains duplicated gene_id values")
        genome.v<-"custom"
    }else{
        genesgr <- get.genesgr(genome.v=genome.v)
    }
    
    # We will create genomic ranges objects with small segmental SVs and left and right breakpoint ends of large SVs and translocations
    if(verbose) message(paste("Generating GRanges objects for SV breaks and SV segments < ",svc.seg.size,"bp",sep=""))
    seg_df  <-  svcdat[intersect(which(svcdat$pos2 - svcdat$pos1 < svc.seg.size), which(svcdat$svclass %in% c("DUP","DEL","BND","INV","INS"))),]
    segmentEventGR = with(seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos2),
                                          rowid = rownames(seg_df),sampleid=seg_df$sample))
    
    non_seg_df <- svcdat[setdiff(rownames(svcdat),rownames(seg_df)),]
    leftJunctionGR = with(non_seg_df, GRanges(chrom1, IRanges(start=pos1, end=pos1),
                                              rowid = rownames(non_seg_df),sampleid=non_seg_df$sample))
    rightJunctionGR = with(non_seg_df, GRanges(chrom2, IRanges(start=pos2, end=pos2),
                                               rowid = rownames(non_seg_df),sampleid=non_seg_df$sample))
    # combine all three GRs 
    svRanges <- c(segmentEventGR,leftJunctionGR,rightJunctionGR)
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN FEATURES AND SVs
    
    if(verbose) message("Finding disrupting breakpoint overlaps")
    disrupt <- geneBreakOverlap(genesgr, svRanges)
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
    
    if(!is.null(upstr)){
        if(verbose) message("Finding upstream breakpoint overlaps")
        upstreamgr <- upgr(genesgr,upstr)
        upstream <- geneBreakOverlap(upstreamgr, svRanges)
    }else{
        upstream <- list(geneSamples=list(),geneBreaks=list())
    }
    
    ###  IDENTIFICATION OF OVERLAPS BETWEEN UPSTREAM REGIONS AND SVs
    if(!is.null(dnstr)){
        if(verbose) message("Finding upstream breakpoint overlaps")
        dnstreamgr <- dngr(genesgr,dnstr)
        dnstream <- geneBreakOverlap(dnstreamgr, svRanges)
    }else{
        dnstream <- list(geneSamples=list(),geneBreaks=list())
    }
    
    return(break.annot(
        input = svcdat,
        genesgr = genesgr,
        disruptSamples = disrupt$geneSamples,
        disruptBreaks = disrupt$geneBreaks,
        upstreamSamples = upstream$geneSamples,
        upstreamBreaks = upstream$geneBreaks,
        dnstreamSamples = dnstream$geneSamples,
        dnstreamBreaks = dnstream$geneBreaks,
        param = list(
            genome.v = genome.v,
            upstr = upstr, 
            dnstr = dnstr, 
            svc.seg.size = svc.seg.size, 
            verbose = verbose)
    ))
    
}


#' Identify recurrently altered genes by CNV. The function will identify overlaps between genomic features (e.g. genes) and CNV  breakpoints. 
#' 
#' As opposed to 'gene.cnv' function that returns the overal CNV of each gene, this function allows identifying sub-genic events and may help detecting other rearrangements.
#' 
#' @param cnv (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean.
#' @param fc.pct (numeric) copy number change between 2 consecutive segments: i.e (default) cutoff = 0.2 represents a fold change of 0.8 or 1.2.
#' @param genome.v (character): either 'hg19' or 'hg38' accepted; reference genome version to retrieve gene annotations including genomic coordinates and strand
#' @param genesgr (S4) a GenomicRanges object containing gene annotations (if not NULL overides genome.v). It is crutial that the genome version 'genesgr' and the input 'sv' are the same. The GRanges object must contain 'strand' and a metadata field 'gene_id' with unique values. Seqnames are expected in the format (chr1, chr2, ...). 
#' @param upstr (numeric) size in base pairs to define gene upstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of upstream regions.
#' @param dnstr (numeric) size in base pairs to define gene downstream region onto which breakpoint overlaps will be identified. The strand value, start and stop positions defined in genesgr will be used to create a GRanges object of downstream regions.
#' @param min.cnv.size (numeric) The minimun segment size (in base pairs) to include in the analysis 
#' @param min.num.probes (numeric) The minimun number of probes per segment to include in the analysis 
#' @param low.cov (data.frame) a data.frame (chr, start, end) indicating low coverage regions to exclude from the analysis
#' @param clean.brk (numeric) Identical segments removal when present in above a given number. Identical CNV segments across multiple samples may represent artifact of common germline variants, this is particularly relevant when the segmentation data was generated with a non-paired reference. For paired datasets (e.g. tumor vs. normal) better leave as NULL.
#' @param verbose (logical) 
#' @return an instance of the class 'break.annot' containing breakpoint mapping onto genes
#' @keywords CNV, segmentation
#' @export
#' @examples
#' 
#' ## validate input data.frame
#' cnv <- validate.cnv(segdat_lung_ccle)
#' 
#' cnv.break.annot(cnv)

cnv.break.annot <- function(cnv, 
                            fc.pct = 0.2, 
                            genome.v="hg19",
                            genesgr=NULL,
                            upstr=50000,
                            dnstr=50000,
                            min.cnv.size = NULL, 
                            min.num.probes = NULL, 
                            low.cov = NULL,
                            clean.brk = NULL,
                            verbose=TRUE){
    
cnvdat <- validate.cnv(cnv)

chrlist <- chr.sort(unique(cnvdat$chrom))
    

if(!is.null(genesgr)){
    if(anyDuplicated(genesgr@elementMetadata$gene_id) > 0) stop("The genesgr provided object contains duplicated gene_id values")
}else{
    genesgr <- get.genesgr(genome.v=genome.v)
}

cnv_breaks <- cnv.breaks(cnvdat, fc.pct = fc.pct, min.cnv.size = min.cnv.size, min.num.probes = min.num.probes, 
                         low.cov = low.cov, clean.brk = clean.brk, verbose=verbose)
    
cnvbrk<- cnv_breaks@breaks
rownames(cnvbrk) <- createRandomString(nrow(cnvbrk),8)
    
breaksgr = with(cnvbrk, GRanges(chrom, IRanges(start=pos, end=pos)))

if(verbose) message("Finding gene disrupting CNV breakpoint")
disrupt <- geneBreakOverlap(genesgr, breaksgr)

    
if(!is.null(upstr)){
    if(verbose) message("Finding upstream CNV breakpoint overlaps")
    upstreamgr <- upgr(genesgr,upstr)
    upstream <- geneBreakOverlap(upstreamgr, breaksgr)
    
}else{
    upstream <- list(geneSamples=list(),geneBreaks=list())
}

if(!is.null(dnstr)){
    if(verbose) message("Finding downstream CNV breakpoint overlaps")
    dnstreamgr <- upgr(genesgr,dnstr)
    dnstream <- geneBreakOverlap(dnstreamgr, breaksgr)
}else{
    dnstream <- list(geneSamples=list(),geneBreaks=list())
}

return(break.annot(
    input= cnvbrk,
    genesgr = genesgr,
    disruptSamples = disrupt$geneSamples,
    disruptBreaks = disrupt$geneBreaks,
    upstreamSamples = upstream$geneSamples,
    upstreamBreaks = upstream$geneBreaks,
    dnstreamSamples = dnstream$geneSamples,
    dnstreamBreaks = dnstream$geneBreaks,
    param = list(
        fc.pct=fc.pct,
        genome.v = genome.v,
        upstr = upstr, 
        dnstr = dnstr,
        min.cnv.size=min.cnv.size,
        min.num.probes=min.num.probes,
        low.cov=low.cov,
        clean.brk=clean.brk)
    ))
    
}

