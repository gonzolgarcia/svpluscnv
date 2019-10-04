## ----include=FALSE-------------------------------------------------------
require(taRifx)
require(data.table)
require(tidyr)
require(circlize)
require(GenomicRanges)
require(D3GB)

## ----message=FALSE-------------------------------------------------------
library(svcnvplus)
dim(brca.ccle.seg)
head(brca.ccle.seg)

dim(brca.ccle.sv)
head(brca.ccle.sv)

## ------------------------------------------------------------------------
segdf <- validate.seg(brca.ccle.seg)

## ------------------------------------------------------------------------
svdf <- validate.sv(brca.ccle.sv)

## ----message=FALSE-------------------------------------------------------
common.breaks <- match.variant.breaks(segdf,svdf,maxgap=50000,low.cov = NULL,verbose=FALSE)

dev.new()
par(mfrow=c(1,2),mar=c(3,9,2,1))
restab <- data.frame(common.breaks$restab)[order(common.breaks$restab$total.seg),]

m1 <- mean(restab$matched.seg/restab$total.seg)
barplot(rbind(restab$matched.seg,restab$total.seg - restab$matched.seg),
        border=NA, las=1, xlab="", horiz=T, cex.main=.7,cex.names=.7, names=rownames(restab),
        main=paste("CNV breaks matched by SV breaks\n","Average =",m1))
grid(ny=NA,nx=NULL)
m2 <- mean(restab$matched.sv/restab$total.sv)
barplot(rbind(restab$matched.sv,restab$total.sv - restab$matched.sv),
        border=NA,las=1,xlab="",horiz=T,cex.main=.7,cex.names=.7, names=rownames(restab),
        main=paste("SV breaks matched by CNV breaks\n","Average =",m2))
grid(ny=NA,nx=NULL)

