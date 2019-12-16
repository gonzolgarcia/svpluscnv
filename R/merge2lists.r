#' 
#'
#' Given two lists with (or without common names) returns a combined list; if common names their values are returned unique
#' @param x list 1
#' @param y list 2
#' @return a combined list 
#' @keywords list 
#' @export
#' @examples
#' merge2lists(x,y)

merge2lists <- function(x,y){
    mergedList <- list()
    for(i in unique(c(names(x),names(y)))){
        if(length(y[[i]]) == 0 & length(x[[i]]) > 0){			mergedList[[i]] <- x[[i]]
        }else if(length(y[[i]]) > 0 & length(x[[i]]) == 0){	mergedList[[i]] <- y[[i]]
        }else if(length(y[[i]]) > 0 & length(x[[i]]) > 0){	mergedList[[i]] <- c(x[[i]],y[[i]])
        }
    }
    mergedList<-lapply(mergedList,unique)
    return(mergedList)
}

