#!/usr/bin/env Rscript

#Rscript deseq2.R ./Results/featureCounts_genes_mod.txt ./Results/sampleInfo.txt ./Results/compsTab.txt mmusculus_gene_ensembl 1 2 1.5 0.1 0 /path/to/results
args = commandArgs(trailingOnly=TRUE)

# Load the libraries
suppressPackageStartupMessages({
	library(DESeq2)
    library(edgeR)
	library(data.table)
	library(openxlsx)
	library(rtracklayer)
	library(stringr)
	library(ggplot2)
})


# Functions
test_output <- function(out) {
    if (dir.exists(out)) {
    write("Output directory already exists!", stdout())
    } else {
    dir.create(out, recursive=TRUE)
    }
}

plot_PCA <- function(obj, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
    plotPCA(obj, intgroup=c("condition"))
	invisible(dev.off())
}

plot_MDS <- function(obj, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
    plotMDS(obj, col=as.numeric(colData(dds)$condition)+1)
	invisible(dev.off())
}

plot_box <- function(obj, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
	#boxplot(countsPerMillion, las=2, col=as.numeric(dgeFull$samples$group)+1, ylab="log-CPM")
    #boxplot(log10(assays(obj)[["cooks"]]), range=0, las=2)
    boxplot(obj, las=2, col=as.numeric(colData(dds)$condition)+1)
	invisible(dev.off())
}

plot_MA <- function(obj, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
    plotMA(obj, ylim=c(-2,2))
	invisible(dev.off())
}

plot_disp <- function(obj, out, fileName) {
	pdf(file.path(out, fileName), paper = "a4r", width = 0, height = 0)
    plotDispEsts(obj)
	invisible(dev.off())
}

dt_list <- function(fun, comps) {
    v <- mapply(FUN=fun, comps$comp, comps$cond1, comps$cond2, SIMPLIFY=FALSE)
    #return(v)
}

do_eT <- function(comp, cond1, cond2) {
    res <- results(dds, contrast=c("condition", cond1, cond2))
    #plot_MA(res, out, paste(cond1, "_vs_", cond2, ".pdf", sep=""))
	dt <- as.data.table(res, keep.rownames="id")
	dt <- merge(dt, anno, all.x=T)
	dt <- dt[order(padj, decreasing=FALSE),]
	return(dt)
}

print_tab <- function(comps, compT) {
	col <- ifelse(fpval == 0, quote(padj), quote(pvalue))
    for (i in comps$comp){
        dt <- compT[[i]]
        write.table(dt, file=file.path(out, paste(i, "_unfiltered.tsv", sep="")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
		write.xlsx(dt, file=file.path(out, paste(i, "_unfiltered.xlsx", sep="")), overwrite=TRUE, asTable=TRUE)
        dt <- dt[(log2FoldChange >= flogFC & eval(col) <= fFDR) | (log2FoldChange <= -flogFC & eval(col) <= fFDR) ,]
		write.table(dt, file=file.path(out, paste(i, "_filtered.tsv", sep="")), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
		write.xlsx(dt, file=file.path(out, paste(i, "_filtered.xlsx", sep="")), overwrite=TRUE, asTable=TRUE)
    }
}

# Create the results folder
out <- file.path(args[10], "deseq2_results")
test_output(out)


# Load the data
data_raw <- read.table(args[1], sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F, row.names=1)
#data_raw <- read.table("/mnt/WD1/test/biomartr/results/subread/featureCounts_sorted.txt", sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F, row.names=1)
#print(as.data.table(data_raw))
sampleInfo <- read.table(args[2], sep='\t', quote='', comment.char='', head=T, stringsAsFactors=TRUE)
#ampleInfo <- read.table("/mnt/WD2/CFFMHS-MH-2324/sampleInfo.txt", sep='\t', quote='', comment.char='', head=T, stringsAsFactors=F)
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
#dgeFull <- DGEList(data_raw, genes=rownames(data_raw), group=sampleInfo$condition)
#dgeFull$samples

#rownames(sampleInfo) <- sampleInfo$name
dds <- DESeqDataSetFromMatrix(countData = data_raw,
                              colData = sampleInfo,
                              design = ~ condition)
#dds
ldds <- log2(assay(dds)+1)
plot_box(ldds, out, 'Boxplot_Pre-Norm.pdf')

#smallestGroupSize <- 3
#dds <- estimateSizeFactors(dds)
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds, normalized=TRUE) >= nCPM) >= nSamples
dds <- dds[keep,]


# Annotation from GTF file
gtf <- as.data.table(readGFF(gtf_path))

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


# Comparisons and Annotation
dds <- DESeq(dds)
#resultsNames(dds)
#res <- results(dds, name="condition_CTR_vs_COND")
#res <- results(dds, contrast=c("condition","CTR","COND"))
compT <- dt_list(do_eT, comps)

# Print complete and filtered tables
print_tab(comps, compT)

#################
# Visualization #
#################
#plot_MA(res, out, 'MA_plot.pdf')
#vsd <- vst(dds, blind=FALSE)
#plot_PCA(vsd, out, 'PCA_plot.pdf')
rld <- rlog(dds, blind = FALSE)
plot_box(assay(rld), out, 'Boxplot_Post-Norm.pdf')
plot_MDS(dds, out, 'MDS_plot.pdf')
plot_disp(dds, out, 'Dispersion_plot.pdf')


sink(file=file.path(out, "sessionInfo.txt"))
sessionInfo()
sink()

write("Done!", stdout())
