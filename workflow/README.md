# STARlight
RNA-seq STAR Pipeline

## How to run:
Load the required modules:

```bash
module load snakemake-7.12.1 STAR-2.7.10a subread-2.0.1 fastqc-0.11.9 R-4.2.1 multiqc-1.9 samtools-1.11 deeptools-3.5.0

## Removed
picard-2.23.8
cutadapt-3.3
```

Run the Snakemake command pointing to the location of the .smk file and the config file:

```bash
snakemake -p --cores 12 -s /path/to/STARlight/STARlight.smk --configfile /path/to/config.json
```

If you need to merge fastq files coming from different lanes you can run *catlanes* ahead of *STARlight*:

```bash
module load snakemake-7.12.1 fastqc-0.11.9 multiqc-1.9
snakemake -p 12 -s catlanes.smk --configfile /path/to/catlanes.json
```
<br>

## Important files:
In order to run the pipeline you need to modify the config file accordingly to your input files. Moreover, you need to provide a file containing the sample information and a table of comparisons to perform.

**sampleInfo file example:**
|name| condition|
|---|---|
| 1_S1|	CTR|
| 2_S2|	LNO|
| 3_S3|	HNO|
| 4_S4|	CTR|
| 5_S5|	LNO|
| 6_S6|	HNO|

**compsTab file example:**
|comp| cond1| cond2 |
|---|---|---|
| LNO_vs_CTR| LNO| CTR|
| HNO_vs_CTR| HNO| CTR|
| HNO_vs_LNO| HNO| LNO|

<br>

## How to select the species:
STARlight works with ENSEMBL. A reference genome and a GTF file for the interested species are needed and could be downloaded via FTP. The config file should be modified to point the right path to these files. The pipeline is using BioMart for the complete annotation so you have to specify the correct species in the config file. The species supported are listed in STARlight/biomaRt/listDatasets.biomaRt.txt. For example:

```bash=
human: hsapiens_gene_ensembl
mouse: mmusculus_gene_ensembl
```

# TO MODIFY!
