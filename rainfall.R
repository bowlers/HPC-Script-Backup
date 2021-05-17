#!/bin/Rscript
library(maftools)
setwd( "/home/sbowler/biocore_lts/scott/U54" )

data <- read.maf(maf="NH.maf")
fam <- subsetMaf(maf=data,genes="FAM186A")
rainfallPlot(maf=fam,detectChangePoints=TRUE,pointSize=0.4,width=12,height=4,savePlot=True)
