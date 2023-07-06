#!/usr/bin/env Rscript

#Rscript edger.R ./Results/featureCounts_genes_mod.txt ./Results/sampleInfo.txt ./Results/compsTab.txt mmusculus_gene_ensembl 1 2 1.5 0.1 0 /path/to/results
args = commandArgs(trailingOnly=TRUE)

# Load the libraries
suppressPackageStartupMessages({
	library(edgeR)
	library(data.table)
	#library(biomaRt)
	library(openxlsx)
	#library(biomartr)
	library(rtracklayer)
	#library(xlsx)
	#library(stringr)
	#library(tidyr)
	#library(ggplot2)
	#library(RColorBrewer)
})


# Functions
test_output <- function(out) {
    if (dir.exists(out)) {
    write("Output directory already exists!", stdout())
    } else {
    dir.create(out, recursive=TRUE)
    }
}

plot_MDS <- function(countsPerMillion, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
	plotMDS(countsPerMillion, col=as.numeric(dgeFull$samples$group)+1)
	invisible(dev.off())
}

plot_box <- function(countsPerMillion, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
	boxplot(countsPerMillion, las=2, col=as.numeric(dgeFull$samples$group)+1, ylab="log-CPM")
	invisible(dev.off())
}

plot_BCV <- function(y, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
	plotBCV(y)
	invisible(dev.off())
}

dt_list <- function(fun, comps) {
    v <- mapply(FUN=fun, comps$comp, comps$cond1, comps$cond2, SIMPLIFY=FALSE)
    #return(v)
}

do_eT <- function(comp, cond1, cond2) {
	#print(paste(cond1, cond2, sep="_"))
	dt <- exactTest(y, pair=c(cond2,cond1))
	dt <- as.data.table(as.data.frame(topTags(dt, n=Inf)))
	#setnames(dt, "genes", "ensembl_gene_id")
	setnames(dt, "genes", "id")
	dt <- merge(dt, anno, all.x=T)
	dt <- dt[order(FDR, decreasing=FALSE),]
	return(dt)
}

edgeR2anno <- function(etObj, annoObj, fileName) {
  t <- as.data.table(as.data.frame(topTags(etObj, n=Inf)))
  #setnames(t, "genes", "ensembl_gene_id")
  setnames(t, "genes", "id")
  t <- merge(t, annoObj, all.x=T)
  t <- t[order(FDR, decreasing=FALSE),]
  write.table(t, file=fileName, row.names=FALSE, sep="\t", quote=FALSE)
  return(t)
}

print_tab <- function(comps, compT) {
	col <- ifelse(fpval == 0, quote(FDR), quote(PValue))
    for (i in comps$comp){
        dt <- compT[[i]]
        write.table(dt, file=file.path(out, paste(i, "_unfiltered.tsv", sep="")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
		write.xlsx(dt, file=file.path(out, paste(i, "_unfiltered.xlsx", sep="")), overwrite=TRUE, asTable=TRUE)
        dt <- dt[(logFC >= flogFC & eval(col) <= fFDR) | (logFC <= -flogFC & eval(col) <= fFDR) ,]
		write.table(dt, file=file.path(out, paste(i, "_filtered.tsv", sep="")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
		write.xlsx(dt, file=file.path(out, paste(i, "_filtered.xlsx", sep="")), overwrite=TRUE, asTable=TRUE)
    }
}


# Create the results folder
#out <- file.path(getwd(), "edgeR_results")
out <- file.path(args[10], "edgeR_results")
test_output(out)


# Load the data
data_raw <- read.table(args[1], sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F, row.names=1)
#print(as.data.table(data_raw))
sampleInfo <- read.table(args[2], sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)
#print(as.data.table(sampleInfo))
comps <- fread(args[3])
#print(comps)
species <- args[4]
nCPM <- as.numeric(args[5])
nSamples <- as.numeric(args[6])
flogFC <- as.numeric(args[7])
fFDR <- as.numeric(args[8])
fpval <- as.numeric(args[9])
attribute <- args[11]
gtf_path <- args[12]


# Create the DGE list
dgeFull <- DGEList(data_raw, genes=rownames(data_raw), group=sampleInfo$condition)
#dgeFull$samples


# DMS before normalization
countsPerMillion <- cpm(dgeFull, log=TRUE)
plot_MDS(countsPerMillion, out, 'MDS_Pre-Norm.pdf')
plot_box(countsPerMillion, out, 'Boxplot_Pre-Norm.pdf')


# Normalization
#countCheck <- countsPerMillion >= 1
#keep <- which(rowSums(countCheck) >= 3)
#dgList <- dgeFull[keep,]
#dim(dgList)

###Shorter code for same Result###
keep <- rowSums(cpm(dgeFull) >= nCPM) >= nSamples
dgList <- dgeFull[keep,]
#dim(dgList)

dgList$samples$lib.size <- colSums(dgList$counts)
#dgList$samples
fdgList <- calcNormFactors(dgList, method="TMM")
robj <- setDT(fdgList$samples, keep.rownames = "sample")
write.xlsx(robj, file.path(out, "edgeR_obj.xlsx"), overwrite=TRUE, asTable=TRUE)
#fdgList$samples


# DMS after normalization
countsPerMillion <- cpm(fdgList, log=TRUE)
plot_MDS(countsPerMillion, out, 'MDS_Post-Norm.pdf')
plot_box(countsPerMillion, out, 'Boxplot_Post-Norm.pdf')


# Print CPM for Heatmaps
countsPerMillion <- cpm(fdgList, log=FALSE)
write.table(countsPerMillion, file=file.path(out, 'raw_counts_cpm_after_norm.txt'), row.names=TRUE, sep="\t", quote=FALSE)
countsPerMillionLog2 <- cpm(fdgList, log=TRUE)
write.table(countsPerMillionLog2, file=file.path(out, 'raw_counts_log2cpm_after_norm.txt'), row.names=TRUE, sep="\t", quote=FALSE)

# Print .xlsx file format
tCPM <- setDT(as.data.frame(countsPerMillion), keep.rownames = "gene")
write.xlsx(tCPM, file.path(out, "raw_counts_cpm_after_norm.xlsx"), overwrite=TRUE, asTable=FALSE)
tCPMlog2 <- setDT(as.data.frame(countsPerMillionLog2), keep.rownames = "gene")
write.xlsx(tCPMlog2, file.path(out, "raw_counts_log2cpm_after_norm.xlsx"), overwrite=TRUE, asTable=FALSE)


# Estimate dispersion
y <- estimateDisp(fdgList)
plot_BCV(y, out, 'BCV_plot.pdf')


# Start biomaRt database
#mart = useMart('ensembl')
# list all the ensembl database of organisms
#listDatasets(mart)  
#choose database of your interest ; in this case its "cfamiliaris_gene_ensembl" I guess
#ensembl = useMart( "ensembl", dataset=species)  
# choose attributes of your interest
#listAttributes(ensembl)
#gene <- getBM( attributes = c("ensembl_gene_id","external_gene_name"), values=as.data.table(topTags(CvsRT, n=Inf))$table.genes, mart=ensembl)  
#anno <- as.data.table(getBM( attributes = c("ensembl_gene_id","external_gene_name","chromosome_name", "start_position", "end_position", "strand", "description"), mart=ensembl))
#anno

#gtf <- as.data.table(readGFF("/mnt/WD1/transcriptR/test/index/anno/anno.gtf"))
#print(attribute)
#print(gtf_path)

gtf <- as.data.table(readGFF(gtf_path))
#print(gtf)

#gtf_g <- gtf[type == "gene",]
#gtf_g <- gtf_g[, c(9, 11, 1, 4, 5, 7, 13)]

#gtf_t <- gtf[type == "transcript",]
#gtf_t <- gtf_t[, c(14, 16, 1, 4, 5, 7, 18)]

if (attribute == "gene_id") {
	gtf_g <- gtf[type == "gene",]
	gtf_g <- gtf_g[, c(9, 11, 1, 4, 5, 7, 13)]
	anno <- gtf_g
	colnames(anno) <- c("id", "name", "chr", "start", "end", "strand", "biotype")
} else {
	gtf_t <- gtf[type == "transcript",]
	gtf_t <- gtf_t[, c(14, 16, 1, 4, 5, 7, 18)]
	anno <- gtf_t
	colnames(anno) <- c("id", "name", "chr", "start", "end", "strand", "biotype")
}

#print(anno)

# Comparisons and Annotation
compT <- dt_list(do_eT, comps)
#compT

# Print complete and filtered tables
print_tab(comps, compT)

sink(file=file.path(out, "sessionInfo.txt"))
sessionInfo()
sink()

write("Done!", stdout())
