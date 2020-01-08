#' 
#'
#' Obtains significance cutoff for the frequency of binary events encoded in a matrix
#' @param mat (numeric matrix) a binary matrix where columns will be tested for their sum value compared to a permutated matrix
#' @param method (character) the method to pass to p.adjust function
#' @param p.cut (numeric) the cutoff for multiple hypothesis corrected p.value  
#' @param iter (numeric) Number of iterations to produce null distribution (note that null size will be iter*ncol(mat))
#' @keywords empirical p.value, p.adjust  
#' @export
#' @examples
#' freq.p.test()



freq.p.test <- function(mat, 
                             method="fdr", 
                             p.cut= 0.05, 
                             iter=100){

stopifnot(is.numeric(mat))
    
# obtain a frequency vector
highDensitiBinsFreq <- apply(mat,2,sum)

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
    param=list(method=method, p.cut= p.cut, iter=iter)
    ))
}


