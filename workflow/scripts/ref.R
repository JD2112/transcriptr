#!/usr/bin/env Rscript

# Rscript ref.R homo_sapiens latest /mnt/WD1/test/biomartr
args = commandArgs(trailingOnly=TRUE)


# Load the libraries
suppressPackageStartupMessages({
	library(data.table)
	library(biomartr)
	library(R.utils)
})


# Functions
test_output <- function(out) {
    if (dir.exists(out)) {
    	#write("Output directory already exists!", stdout())
    } else {
    	dir.create(out, recursive=TRUE)
    }
}


rename_file <- function(i_name, f_name) {
	if(file.exists(i_name)){
		# Rename file name
  		file.rename(i_name, f_name)
	} else {
  		print('File Not found :')
	}
}

make_release <- function(ver) {
	if (ver == "latest") {
		NULL
	} else {
		as.numeric(ver)
	}
}


# Read arguments
species <- args[1]
version <- args[2]
workdir <- args[3]

# Make release
#release <- NULL if version == "latest" else release <- as.numeric(version)
release <- make_release(version)

# Create the results folder
#out <- file.path(getwd(), "edgeR_results")
index <- file.path(args[3], "index")
ref_out <- file.path(index, "ref")
gtf_out <- file.path(index, "gtf")
gff_out <- file.path(index, "gff")

test_output(index)
test_output(ref_out)
test_output(gtf_out)
test_output(gff_out)

# Download biomartr reference files
#getReleases(db="ensembl")
ref <- getGenome(
	db = "ensembl",
	organism = species,
	reference = FALSE,
	release = release,
	gunzip = TRUE,
	path = file.path(ref_out),
	assembly_type = "primary_assembly"
)

gtf <- getGTF(
	db = "ensembl",
	organism = species,
	release = release,
	path = file.path(gtf_out)
	#assembly_type = "primary_assembly"
)

gff <- getGFF(
  db = "ensembl",
  organism = species,
  #reference = FALSE,
  release = release,
  gunzip = FALSE,
  remove_annotation_outliers = FALSE,
  path = file.path(gff_out)
)

# Rename files
gunzip(gtf, remove=TRUE)
un_gtf <- tools::file_path_sans_ext(gtf)

rename_file(ref, file.path(ref_out, "ref.fa"))
rename_file(un_gtf, file.path(gtf_out, "anno.gtf"))
rename_file(gff, file.path(gff_out, "ref.gff3"))
