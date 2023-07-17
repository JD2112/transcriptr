################################
## STARlight RNA-seq Pipeline ##
################################

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
        #gp = config['workdir'] + 'index/genomeParameters.txt',
        idx = config['workdir'] + 'index/STAR/',
        fqc = expand(LOGS + 'fastqc/{sample}.fastqc.log', sample=SAMPLES),
        bam = expand(BAMS + "{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bg = expand(BG + "{sample}.bedgraph", sample=SAMPLES),
        fc = RESULTS + 'subread/featureCounts.txt',
        tab = RESULTS + 'subread/featureCounts_sorted.txt',
        multiqc = RESULTS + 'multiqc_report.html',
        #multiqc = LOGS + 'multiqc/multiqc.log',
        #res = RESULTS + "/edgeR_results" + "/raw_counts_cpm_after_norm.txt",
        res = RESULTS + "edgeR_results" + "/sessionInfo.txt",


rule reference:
    input:
        refR = config['refR'],
    output:
        ref = config['workdir'] + 'index/ref/ref.fa',
        gtf = config['workdir'] + 'index/gtf/anno.gtf'
    params:
        species = config['species'],
        version = config['release'],
        path = config['workdir'],
    shell:
        """
        Rscript {input.refR} {params.species} {params.version} {params.path}
        """


rule index:
    input:
        fa = rules.reference.output.ref,
        gtf = rules.reference.output.gtf,
        #idx = config['workdir'] + 'index/',
    params:
        overhang = config['overhang'],
    output:
        #gp = config['workdir'] + 'index/genomeParameters.txt',
        idx = directory(config['workdir'] + 'index/STAR/'),
    threads:
        12 # set the maximum number of available cores
    shell:
        '/STAR-2.7.10a/source/STAR --runThreadN {threads} ' # docker requirement
#        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output.idx} '
        '--genomeFastaFiles {input.fa} ' #'--genomeFastaFiles <(zcat {input.fa}) '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang {params.overhang}'


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
        ../../fastqc {input.R1} {input.R2} -o {output.out} >> {log} 2>&1 # docker requirement
       #fastqc {input.R1} {input.R2} -o {output.out} >> {log} 2>&1        
        """


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
        '/STAR-2.7.10a/source/STAR --runThreadN {threads} ' # docker requirement
#        'STAR --runThreadN {threads} '        
            '--genomeDir {input.idx} '
            '--readFilesIn <(zcat {input.R1}) <(zcat {input.R2}) '
            '--outSAMtype BAM SortedByCoordinate ' 
            #'--readFilesCommand zcat ' 
            '--quantMode GeneCounts '
            '--outFileNamePrefix {params.out} >> {log} 2>&1'


rule bamCoverage:
    input:
        bam = rules.align_sort.output.bam,
    output:
        bg = BG + '{sample}.bedgraph',
    log:
        LOGS + 'deeptools/{sample}.bamcoverage.log'
    threads:
        12
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage -p {threads} -b {input.bam} -o {output.bg} -of bedgraph --normalizeUsing RPKM
        """


rule featureCounts:
    input:
        bam = expand(BAMS + "{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf = rules.reference.output.gtf,
    output:
        counts = RESULTS + 'subread/featureCounts.txt',
    params:
        stranded = config['stranded'],
        attribute = config['attribute'],
    log:
        LOGS + 'subread/featureCounts.log',
    threads: 
        12
    shell:
        '/subread/bin/featureCounts -a {input.gtf} ' # docker requirement
#        'STAR --runThreadN {threads} '        
        '-g {params.attribute} '
        '-p -s {params.stranded} '
        '-o {output.counts} '
        '-T {threads} {input.bam} >> {log} 2>&1 \n'


rule multiqc:
    input:
        counts = rules.featureCounts.output.counts,
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


rule rearrangeCounts:
    input:
        counts = rules.featureCounts.output.counts,
        sampleInfo = config['sampleInfo']
    output:
        tab = RESULTS + 'subread/featureCounts_sorted.txt',
    params:
        rFC = config['rearrangeFC'],
    threads: 
        12
    shell:
        """
        Rscript {params} {input.counts} {input.sampleInfo} {output.tab}
        """


rule edgeR:
    input:
        #genes = rules.rearrangeCounts.output.tab1,
        #tranx = rules.rearrangeCounts.output.tab2,
        genes = rules.rearrangeCounts.output.tab,
        sampleInfo = config['sampleInfo'],
        compsTab = config['compsTab'],
        gtf = rules.reference.output.gtf,
    output:
        #res = RESULTS + "/edgeR_results" + "/raw_counts_cpm_after_norm.txt",
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