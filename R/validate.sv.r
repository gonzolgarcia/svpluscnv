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
#' @param sv (data.frame) structural variant table including  8 columns: sample, chrom1, pos1, strand1, chrom2, pos2, strand2, svclass
#' @keywords SV, structural variants
#' @export
#' @examples
#' validate.sv()


validate.sv <- function(sv){
  
    require(taRifx,quietly = TRUE,warn.conflicts = FALSE)  # contains remove.factors
  
    stopifnot(ncol(sv) >= 8)
    sv <- remove.factors(data.frame(sv[,1:8]))
    colnames(sv) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass")
    if(length(grep("chr",sv[1,"chrom1"])) == 0) sv[,"chrom1"] <- paste("chr",sv$chrom1,sep="")
    if(length(grep("chr",sv[1,"chrom2"])) == 0) sv[,"chrom2"] <- paste("chr",sv$chrom2,sep="")
  
    stopifnot(is.numeric(sv$pos1))
    stopifnot(is.numeric(sv$pos2))
    stopifnot(is.character(sv$chrom1))
    stopifnot(is.character(sv$chrom2))
    stopifnot(is.character(sv$strand1))
    stopifnot(is.character(sv$strand2))
  
    sv[grep("INV",sv$svclass),"svclass"] <- "INV"
    sv[grep("DUP",sv$svclass),"svclass"] <- "DUP"

    extrachr <- which(unlist(lapply(apply(sv[,c("chrom1","chrom2")],1,unique),length)) == 2) 
    sv[extrachr,"svclass"] <- "TRA"
  
    wrong_class <- setdiff(unique(sv$svclass),c("DEL","DUP","TRA","INV","INS","BND"))
    try(if(length(wrong_class) > 0) message(paste("SV classes not accepted:", paste(wrong_class,collapse=","), "will be set as BND") ))
    sv[which(!sv$svclass %in% c("DEL","DUP","TRA","INV","INS","BND")),"svclass"] <- "BND"

# ensure that pos1 is upstream pos2
intrachr <- which(unlist(lapply(apply(sv[,c("chrom1","chrom2")],1,unique),length)) == 1) 
intrachr_rev <- intersect(which(sv$pos2 -sv$pos1 < 0),intrachr)


if(length(intrachr_rev) > 0){
    svrev <- sv[intrachr_rev,c(1,2,6,7,5,3,4,8)]
    colnames(svrev) <- c("sample","chrom1","pos1","strand1","chrom2","pos2","strand2","svclass")
    sv <- rbind(svrev,sv[setdiff(1:nrow(sv),intrachr_rev),])
    }

stopifnot(nrow(sv) > 0)

  return(sv)

}
