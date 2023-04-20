################################
## STARlight RNA-seq Pipeline ##
################################

configfile:
    "config.json"

workdir:
   config['workdir']

R1SUFFIX = config['R1_suffix']
R2SUFFIX = config['R2_suffix']

SAMPLES, = glob_wildcards(config['data'] + "/{sample}" + R1SUFFIX)
RESULTS = config['workdir'] + '/results'
LOGS = RESULTS + '/logs'
BAMS = RESULTS + '/bams'
BG = RESULTS + '/bedgraphs'


rule all:
    input:
        #idx = directory('index'),
        fqc = expand(LOGS + '/{sample}.fastqc.log', sample=SAMPLES),
        bam = expand(BAMS + "/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bg = expand(BG + "/{sample}.bedgraph", sample=SAMPLES),
        fc1 = RESULTS + '/featureCounts_genes.txt',
        #fc2 = RESULTS + '/featureCounts_transcripts.txt',
        m1 = RESULTS + '/featureCounts_genes_mod.txt',
        #m2 = RESULTS + '/featureCounts_transcripts_mod.txt',
        multiqc = RESULTS + '/multiqc_report.html',
        #res = RESULTS + "/edgeR_results" + "/raw_counts_cpm_after_norm.txt",
        res = RESULTS + "/edgeR_results" + "/sessionInfo.txt",


'''
rule index:
    input:
        fa = config['ref'], # provide your reference FASTA file
        gtf = config['gtf'] # provide your GTF file
    output:
        directory('index') # you can rename the index folder
    threads: 
        12 # set the maximum number of available cores
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {output} '
        '--genomeFastaFiles {input.fa} ' #'--genomeFastaFiles <(zcat {input.fa}) '
        '--sjdbGTFfile {input.gtf} '
        '--sjdbOverhang 100'
'''


rule fastqc:
    input:
        R1 = config['data'] + "/{sample}" + R1SUFFIX,
        R2 = config['data'] + "/{sample}" + R2SUFFIX,
    output:
        out = directory(RESULTS + "/fastqc/{sample}")
    log:
        LOGS + '/{sample}.fastqc.log'
    threads:
        48 # set the maximum number of available cores
    shell:
        """
        mkdir {output.out}
        fastqc {input.R1} {input.R2} -o {output.out} >> {log} 2>&1
        """


rule align_sort:
    input:
        R1 = config['data'] + "/{sample}" + R1SUFFIX,
        R2 = config['data'] + "/{sample}" + R2SUFFIX,
        idx = config['index']
        #idx = rules.index.output
        #idx = directory('index')
    output:
        #temp("{sample}.mapped.bam")
        bam = BAMS + '/{sample}.Aligned.sortedByCoord.out.bam',
        #sj = '{sample}.SJ.out.tab',
        #log = '{sample}.Log.final.out'
    params:
        out = BAMS + '/{sample}.'
    log:
        LOGS + '/{sample}.STAR.log'
    threads:
        12 # set the maximum number of available cores
    shell:
        'STAR --runThreadN {threads} '
            '--genomeDir {input.idx} '
            '--readFilesIn <(zcat {input.R1}) <(zcat {input.R2}) '
            '--outSAMtype BAM SortedByCoordinate ' 
            #'--readFilesCommand zcat ' 
            '--quantMode GeneCounts '
            '--outFileNamePrefix {params.out} >> {log} 2>&1'
        #'mv results/Aligned.sortedByCoord.out.bam {output.bam}'
        #'mv results/SJ.out.tab {output.sj}'
        #'mv results/Log.out {output.log}'


rule bamCoverage:
    input:
        bam = rules.align_sort.output.bam,
    output:
        bg = BG + '/{sample}.bedgraph',
    log:
        #log1 = LOGS + '/{sample}.index.log',
        #log2 = LOGS + '/{sample}.bamCoverage.log'
    threads:
        12
    shell:
        """
        samtools index {input.bam} -@ {threads}
        bamCoverage -p {threads} -b {input.bam} -o {output.bg} -of bedgraph --normalizeUsing RPKM
        """


rule featureCounts:
    input:
        bam = expand(BAMS + "/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        gtf = config['gtf'] # provide your GTF file
    output:
        res1 = RESULTS + '/featureCounts_genes.txt',
        #res2 = RESULTS + '/featureCounts_transcripts.txt'
    params:
        stranded = config['stranded'],
        attribute = config['attribute'],
    log:
        log1 = LOGS + '/featureCounts_genes.log',
        log2 = LOGS + '/featureCounts_transcripts.log'
    threads: 
        12
    shell:
        'featureCounts -a {input.gtf} '
        '-g {params.attribute} '
        '-p -s {params.stranded} '
        '-o {output.res1} '
        '-T {threads} {input.bam} >> {log.log1} 2>&1 \n'


rule multiqc:
    input:
        genes = rules.featureCounts.output.res1,
        fastqc = expand(LOGS + '/{sample}.fastqc.log', sample=SAMPLES),
    output:
        outf = RESULTS + '/multiqc_report.html',
    params:
        indir = RESULTS,
    log:
        log = LOGS + '/multiqc.log',
    threads: 
        12
    shell:
        """
        cd {params.indir}
        multiqc . >> {log} 2>&1
        """


rule rearrangeCounts:
    input:
        genes = rules.featureCounts.output.res1,
        #tranx = rules.featureCounts.output.res2,
        sampleInfo = config['sampleInfo']
    output:
        tab1 = RESULTS + '/featureCounts_genes_mod.txt',
        #tab2 = RESULTS + '/featureCounts_transcripts_mod.txt',
    params:
        rFC = config['rearrangeFC'],
    threads: 
        12
    shell:
        """
        Rscript {params} {input.genes} {input.sampleInfo} {output.tab1}
        """
        #Rscript {params} {input.tranx} {input.sampleInfo} {output.tab2}

'''
rule edgeR:
    input:
        genes = rules.rearrangeCounts.output.tab1,
        #tranx = rules.rearrangeCounts.output.tab2,
        sampleInfo = config['sampleInfo'],
        compsTab = config['compsTab'],
    output:
        #res = RESULTS + "/edgeR_results" + "/raw_counts_cpm_after_norm.txt",
        res = RESULTS + "/edgeR_results" + "/sessionInfo.txt",
    params:
        species = config['species'],
        cpm = config['cpm'],
        nsamp = config['nsamples'],
        logFC = config['logFC'],
        FDR = config['FDR'],
        pval = config['pval'],
        path = RESULTS,
        edgeR = config['edgeR'],
    threads: 
        12
    shell:
        """
        Rscript {params.edgeR} {input.genes} {input.sampleInfo} {input.compsTab} {params.species} {params.cpm} {params.nsamp} {params.logFC} {params.FDR} {params.pval} {params.path}
        """
'''
