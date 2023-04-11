#!/usr/bin/env Rscript

#Rscript rearrangeFeaturesCounts.R ./results/featureCounts_genes.txt ./input/sampleInfo.txt /results/featureCounts_genes_mod.txt
args = commandArgs(trailingOnly=TRUE)

# Load the libraries
suppressPackageStartupMessages({
	library(data.table)
	library(stringr)
})

dt <- fread(args[1], sep='\t', quote='', head=T, stringsAsFactors=F, skip=1)
sampleInfo <- fread(args[2], sep='\t', quote='', head=T, stringsAsFactors=F)

cols <- names(dt[,7:length(names(dt))])
cols_ord <- order(cols) + 6
dt_ord <- c(1, cols_ord)
#setcolorder(dt, dt_ord)
dt <- dt[, ..dt_ord]

fname <- basename(colnames(dt))
sname <- str_split_fixed(fname, "[.]", 2)[,1]
colnames(dt) <- sname

colinfo <- c("Geneid", sampleInfo$name)
setcolorder(dt, colinfo)

write.table(dt, file=args[3], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
