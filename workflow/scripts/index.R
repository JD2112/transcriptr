#!/usr/bin/env Rscript

#Rscript edger.R ./Results/featureCounts_genes_mod.txt ./Results/sampleInfo.txt ./Results/compsTab.txt mmusculus_gene_ensembl 1 2 1.5 0.1 0 /path/to/results
args = commandArgs(trailingOnly=TRUE)


# Load the libraries
suppressPackageStartupMessages({
	library(data.table)
	library(biomartr)
})


# Functions
test_output <- function(out) {
    if (dir.exists(out)) {
    write("Output directory already exists!", stdout())
    } else {
    dir.create(out, recursive=TRUE)
    }
}


# Create the results folder
#out <- file.path(getwd(), "edgeR_results")
out <- file.path(args[10], "index")
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


# Download biomartr reference files
getGTF(
	db = "ensembl",
	organism = "homo_sapiens",
	release = NULL,
	path = file.path("/mnt/WD1/test/biomartr"),
	assembly_type = "primary_assembly"
)



getGenome(
	db = "ensembl",
	organism = "homo_sapiens",
	reference = FALSE,
	release = NULL,
	gunzip = TRUE,
	path = file.path("/mnt/WD1/test"),
	assembly_type = "primary_assembly"
)



mart = useMart('ensembl')
# list all the ensembl database of organisms
#listDatasets(mart)  
#choose database of your interest ; in this case its "cfamiliaris_gene_ensembl" I guess
ensembl = useMart( "ensembl", dataset=species)  
# choose attributes of your interest
#listAttributes(ensembl)
#gene <- getBM( attributes = c("ensembl_gene_id","external_gene_name"), values=as.data.table(topTags(CvsRT, n=Inf))$table.genes, mart=ensembl)  
anno <- as.data.table(getBM( attributes = c("ensembl_gene_id","external_gene_name","chromosome_name", "start_position", "end_position", "strand", "description"), mart=ensembl))
#anno


# Comparisons and Annotation
compT <- dt_list(do_eT, comps)
#compT


# Print complete and filtered tables
print_tab(comps, compT)

sink(file=file.path(out, "sessionInfo.txt"))
sessionInfo()
sink()

write("Done!", stdout())

