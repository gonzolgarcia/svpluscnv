#' 
#'
#' This function validates and reformats the SV (structural variant) calls input. It is used internaly by 'svcnvplus' functions that require this type of data.
#' A few formatting rules are enforced:
#' 1) The input must obtain 8 columns in the following order(sample ID, chromosome of origin, strand of origin, position of origin,, chromosome of destination, strand of destination, position of destination, SV class)
#' 2) SV classes accepted: DEL(deletion), DUP(duplication), INS(insertion), TRA(translocation), INV(inversion) and BND(break end)
#' 3) Any variant in which chromosome of origin and destination differ are encoded as TRA (translocation)
#' 4) pos1 < pos2 is enforced for all variants in which chromosome of origin and destination are the same
#' 5) The class BND can be used to operate with complex events as long as both break ends are the same chromosome
#' 
#' @param svc (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @keywords SV, structural variants
#' @export
#' @examples
#' 
#' validate.svc(SVDATA)


validate.svc <- function(svc){
    
    stopifnot(ncol(svc) >= 8)
    svc <- remove.factors(data.frame(svc[,1:8]))
    colnames(svc) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass")
    if(length(grep("chr",svc[1,"chrom1"])) == 0) svc[,"chrom1"] <- paste("chr",svc$chrom1,sep="")
    if(length(grep("chr",svc[1,"chrom2"])) == 0) svc[,"chrom2"] <- paste("chr",svc$chrom2,sep="")
    
    stopifnot(is.numeric(svc$pos1))
    stopifnot(is.numeric(svc$pos2))
    stopifnot(is.character(svc$chrom1))
    stopifnot(is.character(svc$chrom2))
    stopifnot(is.character(svc$strand1))
    stopifnot(is.character(svc$strand2))
    
    svc[grep("INV",svc$svclass),"svclass"] <- "INV"
    svc[grep("DUP",svc$svclass),"svclass"] <- "DUP"
    
    extrachr <- which(unlist(lapply(apply(svc[,c("chrom1","chrom2")],1,unique),length)) == 2) 
    svc[extrachr,"svclass"] <- "TRA"
    
    wrong_class <- setdiff(unique(svc$svclass),c("DEL","DUP","TRA","INV","INS","BND"))
    try(if(length(wrong_class) > 0) message(paste("SV classes not accepted:", paste(wrong_class,collapse=","), "will be set as BND") ))
    svc[which(!svc$svclass %in% c("DEL","DUP","TRA","INV","INS","BND")),"svclass"] <- "BND"
    
    # ensure that pos1 is upstream pos2
    intrachr <- which(unlist(lapply(apply(svc[,c("chrom1","chrom2")],1,unique),length)) == 1) 
    intrachr_rev <- intersect(which(svc$pos2 -svc$pos1 < 0),intrachr)
    
    
    if(length(intrachr_rev) > 0){
        svcrev <- svc[intrachr_rev,c(1,2,6,7,5,3,4,8)]
        colnames(svcrev) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass")
        svc <- rbind(svcrev,svc[setdiff(1:nrow(svc),intrachr_rev),])
    }
    
    stopifnot(nrow(svc) > 0)
    
    return(svc)
    
}




#' 
#'
#' This function calculates the percent genome changed by CNV
#' @param cnv (data.frame) segmentation data with at least 6 columns: sample, chromosome, start, end, probes, segment_mean
#' @keywords CNV, segmentation
#' @export
#' @examples
#' validate.cnv()


validate.cnv <- function(cnv){
    
    stopifnot(ncol(cnv) >= 6)
    cnv <- remove.factors(data.frame(cnv[,1:6]))
    colnames(cnv) <- c("sample","chrom","start","end","probes","segmean")
    if(length(grep("chr",cnv[1,2])) == 0) cnv[,"chrom"] <- paste("chr",cnv$chrom,sep="")
    stopifnot(is.numeric(cnv$start))
    stopifnot(is.numeric(cnv$end))
    stopifnot(is.numeric(cnv$segmean))
    stopifnot(is.character(cnv$sample))
    stopifnot(is.character(cnv$chrom))
    
    cnv<- cnv[order(cnv$start),]
    options(warn=-1)
    cnv<- cnv[order(as.numeric(gsub("chr","",cnv$chrom))),]
    options(warn=0)
    cnv<- cnv[order(cnv$sample),]
    
    return(cnv)
    
}
