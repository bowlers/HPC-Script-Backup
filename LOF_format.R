#!/bin/Rscript

library(VariantAnnotation)

args <- commandArgs()

print( getwd() )
setwd( "~/biocore_lts/scott")
setwd( paste0( args[6] ) )
setwd( "varscan/lof/filtered" )
print( getwd() )

print( args[7] )
data <- readVcf( paste0( args[7] ), "hg38" )
LOF <- info(data)[, 8:8]
LOF <- gsub( '^.|.$', '', LOF )
x <- sapply(LOF, function(x) strsplit(x, "\\|")[[1]], USE.NAMES=FALSE)
x <- t(x)
write.table(x, paste0( "./lof/", args[7], ".LOF.txt" ), row.names=F, quote=F, col.names=F)
