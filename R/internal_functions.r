#' 
#'
#' Obtains inter quantile average for a defined 'x' vector and both lower and upper quantiles
#' @param x numeric vector 
#' @param lowQ lower quantile
#' @param upQ upper quantile
#' @keywords statistics, interquartile 
#' @export
#' @examples
#' 
#' x <- rnorm(100)
#' IQM(x)


IQM <- function(x, lowQ=0.1, upQ=0.9){
    
    stopifnot(is.numeric(x))
    
    rx <- rank(x,ties.method ='random')
    qt1<-quantile(rx,lowQ)
    qt2<-quantile(rx,upQ)
    
    inter_quantile_mean <- mean(x[intersect(which(rx > qt1),which(rx < qt2))])
    
    return(inter_quantile_mean)
}

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

#' Obtains inter quantile average for a defined 'x' vector and both lower and upper quantiles
#' @param x numeric vector 
#' @param pal color palette
#' @param limits numeric limit fr color mapping
#' @keywords color, number 
#' @export
#' @examples
#' 
#' x <- rnorm(100)
#' pal <- colorRampPalette(c("lightblue","white","salmon"))(256)
#' limits<-c(-1,1)
#' map2color(x,pal,limits)

map2color <- function(x, pal, limits=NULL){
    if(is.null(limits)) limits = range(x)
    return(pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal)+1), all.inside=TRUE)])
}



#' Given two lists with (or without common names) returns a combined list; if common names their values merged and returned as unique
#' @param x (list) input list 1
#' @param y (list) input list 2
#' @param func (character) Either 'unique' or 'intersect' are accepted
#' @return a combined list 
#' @keywords list 
#' @export
#' @examples
#' 
#' x <- sapply(letters[1:10], function(i) sample(1:10)[1:sample(2:10)[1]], simplify=FALSE )
#' y <- sapply(letters[5:15], function(i) sample(1:10)[1:sample(2:10)[1]], simplify=FALSE )
#' merge2lists(x,y)

merge2lists <- function(x,y,fun="unique"){
    mergedList <- list()
    if(fun == "unique"){
        for(i in unique(c(names(x),names(y)))){
            if(length(y[[i]]) == 0 & length(x[[i]]) > 0){
                mergedList[[i]] <- x[[i]]
            }else if(length(y[[i]]) > 0 & length(x[[i]]) == 0){	
                mergedList[[i]] <- y[[i]]
            }else if(length(y[[i]]) > 0 & length(x[[i]]) > 0){
                mergedList[[i]] <- unique(c(x[[i]],y[[i]]))
            }
        }
    }else if(fun == "intersect"){
        for(i in intersect(names(x),names(y)) ){
            commonElements <- intersect(x[[i]],y[[i]])
            if(length(commonElements) > 0){
                mergedList[[i]] <- commonElements
            }
        }
    }else{
        stop(paste("Unknown function:",func) )
    }
    return(mergedList)
}


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


#' retrieves a GRanges object containinng gene annotations for the specified genome version 
#' @param genome.v (hg19 or GRCh37 and hg38 or GRCh38) reference genome version to retrieve gene annotations 
#' @keywords CNV, segmentation, genes
#' @export
#' @examples
#' 
#' get.genesgr(genome.v="hg19")
#' 

get.genesgr<- function(genome.v="hg19",chrlist=NULL){

if(genome.v %in% c("hg19","GRCh37")){
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    genesgr = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
}else if(genome.v %in% c("hg38","GRCh38")){
    require(TxDb.Hsapiens.UCSC.hg38.knownGene)
    genesgr = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
}else{stop("Unspecified, or non available genome")}

genesgr = genesgr[which(as.character(genesgr@seqnames) %in% chrlist)]

err <- capture.output(
    genesgr@elementMetadata$gene_id <- mapIds(org.Hs.eg.db, genesgr@elementMetadata$gene_id,  'SYMBOL','ENTREZID'),
    type="message")

genesgr <- genesgr[which(!is.na(genesgr$gene_id))]
genesgr <- genesgr[which(lapply(genesgr@elementMetadata$gene_id,length) > 0)]
if(is.null(chrlist)) chrlist <- paste("chr",c(1:22,"X","Y"),sep="")
genesgr <- genesgr[which(genesgr@seqnames %in% chrlist)]

return(genesgr)
}


#' A function to order a list of chromosomes 
#' @param chrlist (character): a vector containing chromosome names (chr1, chr2...chrX,chrY  ) 
#' @export
#' @examples
#' 
#' chr.sort(chrlist)
#' 

chr.sort <- function(chrlist){ 
    chrunique <- gsub("chr","",unique(chrlist))
    chrsort <- chrlist[suppressWarnings(order(as.numeric(chrunique) ))]
    return(chrsort)
}




