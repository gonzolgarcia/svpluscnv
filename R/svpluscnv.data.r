#' Class to store refseq  data from UCSC containing exon level info
#' @param data (data.table): transcript information
#' @param exonStarts (list): every transcript exonic start position
#' @param exonStarts (list): every transcript exonic end position
#' @param genome.v (character): the genome version encoding transcript data
#' @return an instance of the class 'refSeqDat' containing transcript exonic coordinates
#' @export

refSeqDat <- setClass("refSeqDat", representation(
    data  = "data.table",
    exonStarts = "list",
    exonEnds= "list",
    genome.v="character"
))

setMethod("show","refSeqDat",function(object){
    writeLines(paste("An object of class refSeqDat from svpluscnv with ",nrow(object@data),"transcipts from",object@genome.v,"genome version"))
})

#' refSeq annotations for hg19 and hg38 versions from UCSC (http://genome.ucsc.edu/cgi-bin/hgTables)
#'
#' @docType data
#' @keywords genes, transcripts, exons, hg19, hg38
#' 
"refseq_hg19"
"refseq_hg38"


#' CCLE degmentation data from breast tissue cell lines
#'
#' Affy SNP6.0 array based segmentation data for lung cancer (DepMap): https://depmap.org/portal/download/
#' TARGET CGI structural variants and CNV segmentation: https://target-data.nci.nih.gov/
#'
#' @docType data
#'
#' @usage data()
#'
#' @keywords CNV segmentation, SVs
"segdat_lung_ccle"
"svdat_lung_ccle"
"cnv_blacklist_regions"
"nbl_segdat"
"nbl_svdat"

