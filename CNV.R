#!/bin/Rscript
.libPaths("/home/sbowler/biocore_lts/scott/R")
.libPaths()
library(DNAcopy)
args <- commandArgs()

setwd( "/home/sbowler/biocore_lts/scott/" )
setwd( paste0( args[6] ) )
setwd( "loh" )
set.seed(0xcafe) # cafe is a hex number

fileName = paste0( args[7], ".copynumber.called" )

cn <- read.table(fileName,header=TRUE)
CNA.object <- CNA( genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio' )
smoothed.CNA.object <- smooth.CNA(CNA.object)
segs <- segment( smoothed.CNA.object, verbose=0, min.width = 2 )

seg.pvalue <- segments.p( segs, ngrid = 100, tol = 1e-6, alpha = 0.05, search.range = 100, nperm = 1000 )
seg.pvalue2 <- seg.pvalue[which(seg.pvalue[,2] %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")),]
seg.pvalue3 <- na.omit(seg.pvalue2)
seg.pvalue4 <- seg.pvalue3[,-1]
seg.pvalue5 <- cbind(ID=seq(1,nrow(seg.pvalue4)),seg.pvalue4)
newFileName = paste0( fileName, ".output" )
write.table(seg.pvalue5, file = newFileName, sep = "\t", quote=F, row.names=F, col.names = T)

detach(package:DNAcopy)
