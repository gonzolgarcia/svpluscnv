#' 
#'
#' Obtains inter quantile average for a defined 'x' vector and both lower and upper quantiles
#' @param x numeric vector 
#' @param lowQ lower quantile
#' @param upQ upper quantile
#' @keywords statistics, interquartile 
#' @export
#' @examples
#' IQM()


IQM <- function(x,lowQ=0.1,upQ=0.9){

  stopifnot(is.numeric(x))

  rx <- rank(x,ties.method ='random')
  qt1<-quantile(rx,lowQ)
  qt2<-quantile(rx,upQ)
  inter_quantile_mean <- mean(x[intersect(which(rx > qt1),which(rx < qt2))])
  return(inter_quantile_mean)
}

