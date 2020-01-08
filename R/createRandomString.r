#' 
#'
#' Generates n unique random character strings of a given length
#' @param n character vector length to return 
#' @param strlen random string length
#' @param seed set.seed value
#' @keywords random string 
#' @export
#' @examples
#' 
#' createRandomString(n=10 strlen=12)


createRandomString <- function(n=1, strlen=12, seed=123456789){

strlenchain <- strlen*n*2

set.seed(seed)
chain <- paste(sample(c(letters, LETTERS),strlenchain, replace=TRUE),collapse="")
idresult <- strsplit(gsub(paste("(.{",strlen,"})",sep=""), "\\1 ", chain)," ")

if(anyDuplicated(idresult[[1]]) != 0) stop("Repeated strings were produced; try modifying the 'seed' or increasing 'strlen'")

return(idresult[[1]][1:n])
}

