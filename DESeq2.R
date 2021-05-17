#!/bin/Rscript
library(DESeq2)
args <- commandArgs()

setwd( ".." )
setwd( paste0( args[6] ) )
setwd( "RNA/htseq" )

sampleFiles <- grep("tsv",list.files("."), value=TRUE)
sampleCondition <- gsub(".*\\.(.*)\\..*", "\\1", sampleFiles )
sampleTable <- data.frame(sampleName = sampleFiles,
		fileName = sampleFiles,
		condition = sampleCondition)
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount( sampleTable = sampleTable,
			directory = ".",
			design=~ condition)
dds <- DESeq(ddsHTSeq)

res <- results(dds)

write.csv( res, file="Results_DESeq2.csv" )
