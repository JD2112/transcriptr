#!/usr/bin/env Rscript

# Rscript t2g.R /path/to/gtf path/to/results/t2g.tsv
args = commandArgs(trailingOnly=TRUE)


# Load the libraries
suppressPackageStartupMessages({
	library(data.table)
	library(rtracklayer)
})


# Read arguments
gtf <- args[1]
out <- args[2]

# Match trnascripts stable version
gtf <- as.data.table(readGFF(gtf))
gtf <- gtf[!is.na(gtf$transcript_version),]
gtf$transcript_id <- paste(gtf$transcript_id, gtf$transcript_version, sep=".")

# Print out conversion table
t2g <- gtf[, c(14, 9)]
t2g <- unique(t2g)
fwrite(t2g, file=out, quote=FALSE, sep='\t')
