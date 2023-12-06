#!/usr/bin/env Rscript

#Rscript rearrangeSalmon.R /path/to/input ./input/sampleInfo.txt path/to/tx2g.tsv gene_id path/to/output
args = commandArgs(trailingOnly=TRUE)

# Load the libraries
suppressPackageStartupMessages({
	library(tximport)
	library(data.table)
	library(stringr)
})


# Read arguments
dir <- args[1]
sampleInfo <- fread(args[2], sep='\t', quote='', head=T, stringsAsFactors=F)
tx2g <- fread(args[3])
attribute <- args[4]
#> list.files(dir)
#[1] "1_S1" "2_S2" "3_S3" "4_S4" "5_S5" "6_S6"

# Find all salmon results
samples <- list.files(dir)
files <- file.path(dir, samples, "quant.sf")
names(files) <- samples
ns <- length(files)

print(files)

if(attribute == 'gene_id'){
	txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2g)
} else {
	txi.salmon <- tximport(files, type = "salmon", txOut = TRUE)
}


#tab <- txi.salmon$counts[,1:ns]
#dt <- as.data.table(txi.salmon$counts[,1:ns], keep.rownames="Geneid")
dt <- as.data.table(round(txi.salmon$counts[,1:ns]), keep.rownames="Geneid")

cols <- names(dt[,2:length(names(dt))])
cols_ord <- order(cols) + 1
dt_ord <- c(1, cols_ord)
#setcolorder(dt, dt_ord)
dt <- dt[, ..dt_ord]

#fname <- basename(colnames(dt))
#sname <- str_split_fixed(fname, "[.]", 2)[,1]
#colnames(dt) <- sname

colinfo <- c("Geneid", sampleInfo$name)
setcolorder(dt, colinfo)

write.table(dt, file=args[5], row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
