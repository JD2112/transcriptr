####################################################
## TranscriptR Salmon paired-end RNA-seq Pipeline ##
####################################################

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
        #log = LOGS + 'salmon/{sample}.salmon.log',
        res = expand(RESULTS + "salmon_quant/{sample}/", sample=SAMPLES),
        tab = RESULTS + 'salmon_quant/salmonCounts_sorted.txt',
        #bam = expand(BAMS + "{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        #bg = expand(BG + "{sample}.bedgraph", sample=SAMPLES),
        #fc = RESULTS + 'subread/featureCounts.txt',
        #tab = RESULTS + 'subread/featureCounts_sorted.txt',
        multiqc = RESULTS + 'multiqc_report.html',
        edgeR = RESULTS + "edgeR_results" + "/sessionInfo.txt",
        DESeq2 = RESULTS + "deseq2_results" + "/sessionInfo.txt",
        limma = RESULTS + "limma_results" + "/sessionInfo.txt",


rule reference:
    input:
        refR = config['refR'],
    output:
        ref = config['workdir'] + 'index/ref/ref.fa',
        cds = config['workdir'] + 'index/rna/cds.fa',
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


rule mergeRef:
    input:
        cds = rules.reference.output.cds,
        rna = rules.reference.output.rna,
    output:
        fa = config['workdir'] + 'index/rna/cdna.ncrna.fa'
    shell:
        """
        cat {input.cds} {input.rna} > {output.fa}
        """


rule salmon_t2g:
    input:
        gtf = rules.reference.output.gtf,
    output:
        t2g = config['workdir'] + 'index/gtf/t2g.tsv',
    params:
        script = config['t2g']
    threads:
        12 # set the maximum number of available cores
    shell:
        """
        Rscript {params.script} {input.gtf} {output.t2g}
        """


rule salmon_index:
    input:
        fa = rules.mergeRef.output.fa,
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
        /salmon-1.10.2/bin/salmon index -p {threads} -t {input.fa} -i {output.idx}
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
        ../../FastQC/fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1
        """

#  ../../FastQC/fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1# docker requirement
# fastqc {input.R1} {input.R2} -t {threads} -o {output.out} >> {log} 2>&1 # without docker run


rule salmon_quant:
    input:
        R1 = config['data'] + "{sample}" + R1SUFFIX,
        R2 = config['data'] + "{sample}" + R2SUFFIX,
        idx = rules.salmon_index.output.idx,
        #gtf = rules.reference.output.gtf,
        #gtf = rules.salmon_t2g.output.t2g,
    output:
        res = directory(RESULTS + "salmon_quant/{sample}/")
    params:
        lib = "A"
    log:
        log = LOGS + 'salmon/{sample}.salmon.log'
    threads:
        12 # set the maximum number of available cores
    shell:
        """
        /salmon-1.10.2/bin/salmon quant -p {threads} -i {input.idx} -l {params.lib} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o {output.res}
        """

#/salmon-1.10.2/bin/salmon quant -p {threads} -i {input.idx} -l {params.lib} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o {output.res} # docker requirement
#salmon quant -p {threads} -i {input.idx} -l {params.lib} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o {output.res} # without docker run

#Attribute choice has been moved into the rearrangeCounts rule
#salmon quant -p {threads} -i {input.idx} -l {params.lib} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o {output.res} -g {input.gtf}


rule rearrangeCounts:
    input:
        res = expand(RESULTS + "salmon_quant/{sample}/", sample=SAMPLES),
        sampleInfo = config['sampleInfo'],
        t2g = rules.salmon_t2g.output.t2g,
    output:
        tab = RESULTS + 'salmon_quant/salmonCounts_sorted.txt',
    params:
        rS = config['rearrangeSalmon'],
        attribute = config['attribute'],
        path = RESULTS + "salmon_quant/",
    threads: 
        12
    shell:
        """
        Rscript {params.rS} {params.path} {input.sampleInfo} {input.t2g} {params.attribute} {output.tab}
        """


rule multiqc:
    input:
        counts = rules.rearrangeCounts.output.tab,
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


rule edgeR:
    input:
        genes = rules.rearrangeCounts.output.tab,
        sampleInfo = config['sampleInfo'],
        compsTab = config['compsTab'],
        gtf = rules.reference.output.gtf,
    output:
        res = RESULTS + "edgeR_results" + "/sessionInfo.txt",
    params:
        species = config['species'],
        cpm = config['cpm'],
        nsamp = config['nsamples'],
        logFC = config['logFC'],
        FDR = config['FDR'],
        pval = config['pval'],
        path = RESULTS,
        edgeR = config['edgeR'],
        attribute = config['attribute'],
    threads: 
        12
    shell:
        """
        Rscript {params.edgeR} {input.genes} {input.sampleInfo} {input.compsTab} {params.species} {params.cpm} {params.nsamp} {params.logFC} {params.FDR} {params.pval} {params.path} {params.attribute} {input.gtf} 
        """


rule DESeq2:
    input:
        genes = rules.rearrangeCounts.output.tab,
        sampleInfo = config['sampleInfo'],
        compsTab = config['compsTab'],
        gtf = rules.reference.output.gtf,
    output:
        res = RESULTS + "deseq2_results" + "/sessionInfo.txt",
    params:
        species = config['species'],
        cpm = config['cpm'],
        nsamp = config['nsamples'],
        logFC = config['logFC'],
        FDR = config['FDR'],
        pval = config['pval'],
        path = RESULTS,
        deseq2 = config['DESeq2'],
        attribute = config['attribute'],
    threads: 
        12
    shell:
        """
        Rscript {params.deseq2} {input.genes} {input.sampleInfo} {input.compsTab} {params.species} {params.cpm} {params.nsamp} {params.logFC} {params.FDR} {params.pval} {params.path} {params.attribute} {input.gtf} 
        """


rule limma:
    input:
        genes = rules.rearrangeCounts.output.tab,
        sampleInfo = config['sampleInfo'],
        compsTab = config['compsTab'],
        gtf = rules.reference.output.gtf,
    output:
        res = RESULTS + "limma_results" + "/sessionInfo.txt",
    params:
        species = config['species'],
        cpm = config['cpm'],
        nsamp = config['nsamples'],
        logFC = config['logFC'],
        FDR = config['FDR'],
        pval = config['pval'],
        path = RESULTS,
        limma = config['limma'],
        attribute = config['attribute'],
    threads: 
        12
    shell:
        """
        Rscript {params.limma} {input.genes} {input.sampleInfo} {input.compsTab} {params.species} {params.cpm} {params.nsamp} {params.logFC} {params.FDR} {params.pval} {params.path} {params.attribute} {input.gtf} 
        """
