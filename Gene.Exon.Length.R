#!/bin/Rscript
# Adapted from Irsan 2014
# Used to calculate non-overlapping exon lengths for downstream TMB based on 315 genes in FoundationMedicine panel

library(GenomicFeatures)
txdb <- makeTxDbFromGFF("yourFile.gtf",format="gtf")

# then collect the exons per gene id
exons.list.per.gene <- exonsBy(txdb,by="gene")

# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- sum(width(reduce(exons.list.per.gene)))
write.table(exonic.gene.sizes,"exon.sizes", sep="\t")

# Use python's myGene to convert Ensembl IDs to symbols.
# TMB = sum(mutations over each of the 315 genes) / sum( exonic lengths of 315 genes)
