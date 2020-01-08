#' 
#'
#' Obtains inter quantile standard deviation for a defined 'x' vector and both lower and upper quantiles
#' @param x numeric vector 
#' @param lowQ lower quantile
#' @param upQ upper quantile
#' @keywords statistics, interquartile 
#' @export
#' @examples
#' 
#' x <- rnorm(100)
#' IQSD(x)


IQSD <- function(x,lowQ=0.1,upQ=0.9){
  stopifnot(is.numeric(x))
  
  rx <- rank(x,ties.method ='random')
  qt1<-quantile(rx,lowQ)
  qt2<-quantile(rx,upQ)
  
  inter_quantile_mean <- sd(x[intersect(which(rx > qt1),which(rx < qt2))])
  return(inter_quantile_mean)

  }

