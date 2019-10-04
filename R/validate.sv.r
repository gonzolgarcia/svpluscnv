#' 
#'
#' This function calculates the percent genome changed by CNV
#' sv classes accepted: DEL(deletion), DUP(duplication), h2hINV(head to head inversion), t2tINV(tail to tail inversion), TRA(translocation), INV(inversion)
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
  
wrong_class <- setdiff(unique(sv$svclass),c("DEL","DUP","h2hINV","t2tINV","TRA","INV"))
try(if(length(wrong_class) > 0) message(paste("SV classes not accepted:", paste(wrong_class,collapse=",")) ))
sv <- sv[which(sv$svclass %in% c("DEL","DUP","h2hINV","t2tINV","TRA","INV")),]

stopifnot(nrow(sv) > 0)

  return(sv)

}
