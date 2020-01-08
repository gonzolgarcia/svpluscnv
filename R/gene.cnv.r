#' 
#'
#' Obtain a matrix with the weighted average CN per chromosome arm 
#' @param seg (data.frame) segmentation data with 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @param genome.v (hg19 or hg38) reference genome version to draw chromosome limits and centromeres
#' @param chrlist (character) list of chromosomes to include chr1, chr2, etc...
#' @param geneid (character) only "Symbol" accepted; if NULL, entrez ID will be used
#' @param fill.gaps (logical) whether to fill the gaps in the segmentation file using gap neighbour segmean average as log ratio
#' @param verbose (logical) 
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' gene.cnv(seg)


gene.cnv <- function(seg, 
                     genome.v="hg19",
                     chrlist=NULL, 
                     geneid="Symbol",
                     fill.gaps=TRUE,
                     verbose=TRUE){

if(is.null(chrlist)) chrlist <- paste("chr",c(1:22,"X"), sep="" )

segdat <- validate.seg(seg)

if(fill.gaps) segdat <- segment.gap(seg,chrlist=chrlist,verbose=verbose)

if(genome.v %in% c("hg19","GRCh37")){
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db)
  genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  if(geneid== "Symbol") genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID')
  genesgr <- genesgr[which(!is.na(genesgr$gene_id)),]
}else if(genome.v %in% c("hg38","GRCh38")){
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(org.Hs.eg.db)
  genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  if(geneid== "Symbol") genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID')
  genesgr <- genesgr[which(!is.na(genesgr$gene_id)),]
}else{stop("Unspecified, or non available genome")}

genes <- remove.factors(as.data.frame(genesgr))
genes <- genes[which(genes$seqnames %in% chrlist),]
genes <- genes[which(!is.na(genes$gene_id)),]
genes <- genes[order(genes$start),]
genes <- genes[order(genes$seqnames),]
rownames(genes) <- genes$gene_id


geneLimits_gr <- with(genes, GRanges(seqnames, IRanges(start=start, end=end)))
segdat_gr <- with(segdat, GRanges(chrom, IRanges(start=start, end=end)))

hits <-GenomicAlignments::findOverlaps(geneLimits_gr,segdat_gr)

overlaps_all <- pintersect(geneLimits_gr[queryHits(hits),], segdat_gr[subjectHits(hits),])
width_overlap <- width(overlaps_all)


df<-data.frame(segdat[subjectHits(hits),c("sample","segmean")],genes[queryHits(hits),"gene_id"],width_overlap)
colnames(df) <- c("sample","segmean","gene_id","width")

cnvmat <- matrix(ncol=length(unique(segdat$sample)), nrow=nrow(genes) )
colnames(cnvmat) <- unique(segdat$sample)
rownames(cnvmat) <- genes$gene_id

if(verbose){
  message("Calculating gene level CNV")
  pb <- txtProgressBar(style=3)
  cc <-0
  tot <- ncol(cnvmat)
}
for(i in unique(segdat$sample)){ 
  dfi <- df[which(df$sample == i),]
  num <- aggregate(segmean~gene_id,dfi,mean)
  
  gene_cn <- as.numeric(num[,2])
  names(gene_cn) <- as.character(num[,1])
  cnvmat[names(gene_cn),i] <- gene_cn
  
  if(verbose) cc <- cc+1
  if(verbose) setTxtProgressBar(pb, cc/tot)
  }
if(verbose) close(pb)
cnvmat<- na.omit(cnvmat)
return(list(
  cnvmat=cnvmat,
  genesdf=genes[rownames(cnvmat),],
  seg=segdat))
}



