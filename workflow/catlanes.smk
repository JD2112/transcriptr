##########################################################
## cat fastq files from different lanes and performs QC ##
##########################################################


configfile:
    "catlanes.json"

workdir:
   config['workdir']

R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']

R1_S1 = config['R1_L001']
R1_S2 = config['R1_L002']
R1_S3 = config['R1_L003']
R1_S4 = config['R1_L004']
R2_S1 = config['R2_L001']
R2_S2 = config['R2_L002']
R2_S3 = config['R2_L003']
R2_S4 = config['R2_L004']

SAMPLES, = glob_wildcards(config['data'] + "/{sample}" + R1_S1)
DATA = config['data']
RESULTS = config['workdir'] + '/catlanes'
FASTQ = RESULTS + '/fastq'
LOGS = RESULTS + '/logs'


rule all:
    input:
        #idx = directory('index'),
        fq1 = expand(FASTQ + '/{sample}' + R1SUFFIX, sample=SAMPLES),
        fq2 = expand(FASTQ + '/{sample}' + R2SUFFIX, sample=SAMPLES),
        fqc = expand(LOGS + '/{sample}.fastqc.log', sample=SAMPLES),
        multiqc = RESULTS + '/multiqc_report.html',

rule catlanes:
    input:
        R1_1 = DATA + "/{sample}" + R1_S1,
        R1_2 = DATA + "/{sample}" + R1_S2,
        R1_3 = DATA + "/{sample}" + R1_S3,
        R1_4 = DATA + "/{sample}" + R1_S4,
        R2_1 = DATA + "/{sample}" + R2_S1,
        R2_2 = DATA + "/{sample}" + R2_S2,
        R2_3 = DATA + "/{sample}" + R2_S3,
        R2_4 = DATA + "/{sample}" + R2_S4,
    output:
        R1 = FASTQ + "/{sample}" + R1SUFFIX,
        R2 = FASTQ + "/{sample}" + R2SUFFIX,
    log:
        LOGS + '/{sample}.cat.log'
    threads:
        24 # set the maximum number of available cores
    shell:
        """
        cat {input.R1_1} {input.R1_2} {input.R1_3} {input.R1_4} > {output.R1}
        cat {input.R2_1} {input.R2_2} {input.R2_3} {input.R2_4} > {output.R2}
        """


rule fastqc:
    input:
        R1 = rules.catlanes.output.R1,
        R2 = rules.catlanes.output.R2,
    output:
        out = directory(RESULTS + "/fastqc/{sample}")
    log:
        log = LOGS + '/{sample}.fastqc.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        fastqc {input.R1} {input.R2} -o {output.out} >> {log} 2>&1
        """


rule multiqc:
    input:
        fastqc = expand(LOGS + '/{sample}.fastqc.log', sample=SAMPLES),
    output:
        outf = RESULTS + '/multiqc_report.html',
    params:
        indir = RESULTS,
    log:
        log = LOGS + '/multiqc.log',
    threads: 
        24
    shell:
        """
        cd {params.indir}
        multiqc .
        """
