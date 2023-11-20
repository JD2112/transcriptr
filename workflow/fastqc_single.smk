############################
## Performs single-end QC ##
############################


configfile:
    "fastqc_single.json"


workdir:
   config['workdir']

R1SUFFIX = config['R1_suffix']
#R2SUFFIX = config['R2_suffix']

SAMPLES, = glob_wildcards(config['data'] + "{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + 'results_qc/'
LOGS = RESULTS + 'logs/'


rule all:
    input:
        fq1 = expand(RESULTS + 'fastqc/{sample}', sample=SAMPLES),
        #fq2 = expand(RESULTS + 'fastqc/{sample}', sample=SAMPLES),
        fqc = expand(LOGS + 'fastqc/{sample}.fastqc.log', sample=SAMPLES),
        multiqc = RESULTS + 'multiqc_report.html',


rule fastqc:
    input:
        R1 = config['data'] + "{sample}" + R1SUFFIX,
        #R2 = config['data'] + "{sample}" + R2SUFFIX,
    output:
        out = directory(RESULTS + "fastqc/{sample}")
    log:
        LOGS + 'fastqc/{sample}.fastqc.log'
    threads:
        12 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        ../../FastQC/fastqc {input.R1} -t {threads} -o {output.out} >> {log} 2>&1
        """


rule multiqc:
    input:
        fastqc = expand(LOGS + 'fastqc/{sample}.fastqc.log', sample=SAMPLES),
    output:
        outf = RESULTS + 'multiqc_report.html',
    params:
        indir = RESULTS,
    log:
        log = LOGS + 'multiqc/multiqc.log',
    threads: 
        12
    shell:
        """
        cd {params.indir}
        multiqc . >> {log} 2>&1
        """

