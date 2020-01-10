#' 
#'
#' Obtains significance cutoff for the frequency of binary events encoded in a matrix
#' @param mat (numeric matrix) a binary matrix where columns will be tested for their sum value compared to a permutated matrix
#' @param method (character) the method to pass to p.adjust function
#' @param p.cut (numeric) the cutoff for multiple hypothesis corrected p.value  
#' @param zerofreq (logical) whether to remove bins with observed frequency = 0; It is recommended to set to TRUE when the bins span genomic regions of low coverage   
#' @param iter (numeric) Number of iterations to produce null distribution (note that null size will be iter*ncol(mat))
#' @keywords empirical p.value, p.adjust  
#' @export
#' @examples
#' 
#' ## validate input data.frames
#' seg <- validate.seg(segdat_lung_ccle)
#' 
#' ## obtain a matrix of genomic bins vs samples indicating high density of breaks
#' shatt.regions <- shattered.regions.cnv(seg)
#' mat <- shatt.regions$high.density.regions.hc
#' 
#' freq.p.test(mat)



freq.p.test <- function(mat, 
                        method="fdr", 
                        p.cut= 0.05,
                        iter=100,
                        zerofreq=TRUE,
                        verbose=TRUE){

stopifnot(is.numeric(mat))

# obtain a frequency vector
highDensitiBinsFreq <- apply(mat,2,sum)

if(zerofreq == TRUE){
    bins.nozero <- names(which(highDensitiBinsFreq > 0))
    mat <- mat[,bins.nozero]
    highDensitiBinsFreq <- highDensitiBinsFreq[bins.nozero]
    if(verbose) message( paste("Testing ",dim(mat)[2],"non-zero bins in ",dim(mat)[1], "samples") )
}else{
    if(verbose) message( paste("Testing ",dim(mat)[2],"bins in ",dim(mat)[1], "samples") )
}


# create null distribution by sample shuffling
highDensitiBinsFreqRandomFreq<-list()
for(i in 1:iter){
    highDensitiBinsRandom<- t(apply(mat,1,sample))
    highDensitiBinsFreqRandomFreq[[i]] <- apply(highDensitiBinsRandom,2,sum)
    }
highDensitiBinsFreqRandomFreqNull <- unlist(highDensitiBinsFreqRandomFreq)

# obtain the frequency cutoff for statistical significance (e.g. FDR < 0.01)
pvalues <- highDensitiBinsFreq
for(i in 0:max(highDensitiBinsFreq)){
    pvalues[which(highDensitiBinsFreq == i)] <- length(which(highDensitiBinsFreqRandomFreqNull >i))/ length(highDensitiBinsFreqRandomFreqNull)
    }

freq.cut <- min(highDensitiBinsFreq[names(which(p.adjust(pvalues, method=method) < p.cut))])

return(list(
    freq.cut = freq.cut,
    pvalues = pvalues,
    observed = highDensitiBinsFreq,
    null = highDensitiBinsFreqRandomFreqNull,
    param = list(method=method, p.cut= p.cut, iter=iter)
    ))
}


