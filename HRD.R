#!/bin/Rscript
.libPaths("/home/sbowler/biocore_lts/scott/R_CNV")
.libPaths()
library(scarHRD)
args <- commandArgs()

setwd( "/home/sbowler/biocore_lts/scott/" )
setwd( paste0( args[6] ) )
setwd( "HRD" )

fileName = paste0( args[7], ".small.seqz.gz" )

scar_score(fileName, reference="grch38", seqz=TRUE)
