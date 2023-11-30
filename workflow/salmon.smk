#############################################
## TranscriptR paired-end RNA-seq Pipeline ##
#############################################

configfile:
    "config.json"

workdir:
   config['workdir']

R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']

SAMPLES, = glob_wildcards(config['data'] + "{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + 'results/'
LOGS = RESULTS + 'logs/'
BAMS = RESULTS + 'bams/'
BG = RESULTS + 'bedgraphs/'


rule all:
    input:
        ref = config['workdir'] + 'index/ref/ref.fa',
        gtf = config['workdir'] + 'index/gtf/anno.gtf',
        #faidx = config['workdir'] + 'index/ref/ref.fa.fai',
        idx = config['workdir'] + 'index/salmon/',
        fqc = expand(LOGS + 'fastqc/{sample}.fastqc.log', sample=SAMPLES),
        #bam = expand(BAMS + "{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        #bg = expand(BG + "{sample}.bedgraph", sample=SAMPLES),
        #fc = RESULTS + 'subread/featureCounts.txt',
        #tab = RESULTS + 'subread/featureCounts_sorted.txt',
        #multiqc = RESULTS + 'multiqc_report.html',
        #edgeR = RESULTS + "edgeR_results" + "/sessionInfo.txt",
        #DESeq2 = RESULTS + "deseq2_results" + "/sessionInfo.txt",
        #limma = RESULTS + "limma_results" + "/sessionInfo.txt",


rule reference:
    input:
        refR = config['refR'],
    output:
        ref = config['workdir'] + 'index/ref/ref.fa',
        rna = config['workdir'] + 'index/rna/rna.fa',
        gtf = config['workdir'] + 'index/gtf/anno.gtf',
        gff = config['workdir'] + 'index/gff/ref.gff3',
    params:
        species = config['species'],
        version = config['release'],
        path = config['workdir'],
    shell:
        """
        Rscript {input.refR} {params.species} {params.version} {params.path}
        """


rule salmon_index:
    input:
        fa = rules.reference.output.rna,
        gtf = rules.reference.output.gtf,
        #idx = config['workdir'] + 'index/',
    params:
        overhang = config['overhang'],
    output:
        #gp = config['workdir'] + 'index/genomeParameters.txt',
        idx = directory(config['workdir'] + 'index/salmon/'),
    threads:
        12 # set the maximum number of available cores
    shell:
        """
        salmon index -p {threads} -t {input.fa} -i {output.idx}
        """


rule fastqc:
    input:
        R1 = config['data'] + "{sample}" + R1SUFFIX,
        R2 = config['data'] + "{sample}" + R2SUFFIX,
    output:
        out = directory(RESULTS + "fastqc/{sample}")
    log:
        LOGS + 'fastqc/{sample}.fastqc.log'
    threads:
        12 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1
        """

#  ../../FastQC/fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1# docker requirement
# fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1 # without docker run


rule align_sort:
    input:
        R1 = config['data'] + "{sample}" + R1SUFFIX,
        R2 = config['data'] + "{sample}" + R2SUFFIX,
        idx = rules.index.output.idx
    output:
        bam = BAMS + '{sample}.Aligned.sortedByCoord.out.bam',
    params:
        out = BAMS + '{sample}.'
    log:
        LOGS + 'star/{sample}.STAR.log'
    threads:
        12 # set the maximum number of available cores
    shell:
#        '/STAR-2.7.10a/source/STAR --runThreadN {threads} ' # docker requirement
        'STAR --runThreadN {threads} '        
            '--genomeDir {input.idx} '
            '--readFilesIn <(zcat {input.R1}) <(zcat {input.R2}) '
            '--outSAMtype BAM SortedByCoordinate ' 
            #'--readFilesCommand zcat ' 
            '--quantMode GeneCounts '
            '--outFileNamePrefix {params.out} >> {log} 2>&1'
