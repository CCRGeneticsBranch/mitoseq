# Fastqc
rule fastqc:
    """
    Run FastQC on uBams
    """
    input:
        ubam=join(RESULTSDIR,"uBam","{sample}.unmapped.bam")
    output:
        html = join(RESULTSDIR,"qc","fastqc","{sample}.unmapped_fastqc.html")
    params:
        fastqc=join(RESULTSDIR,"qc","fastqc")
    envmodules:
        TOOLS["fastqc"]
    shell:
        """
        fastqc -o {params.fastqc} -f bam -t $SLURM_CPUS_PER_TASK {input.ubam}
        """

# Samtools bam stats
rule samtools_stats:
    """
    Run samtools flagstat and idxstats
    """
    input:
        bam=join(RESULTSDIR,"bam","{sample}.mapped.bam")
    output:
        flagstat=join(RESULTSDIR,"qc","samtools","{sample}.flagstat"),
        idxstats=join(RESULTSDIR,"qc","samtools","{sample}.idxstats")
    envmodules:
        TOOLS["samtools"]
    shell:
        """
        samtools flagstat {input.bam} > {output.flagstat}
        samtools idxstats {input.bam} > {output.idxstats}
        """

rule create_interval_list:
    """
    Create interval lists from bait and target files
    """
    input:
        bait = config[genome]['bait'],
        target = config[genome]['target']
    output:
        bait = join(RESULTSDIR,"bait.interval"),
        target = join(RESULTSDIR,"target.interval")
    params:
        ref = config[genome]['ref']
    envmodules:
        TOOLS["picard"]
    shell:
        """
        java -Xmx4g -jar $PICARDJARPATH/picard.jar BedToIntervalList \
            I={input.bait} \
            O={output.bait} \
            SD={params.ref}
        
        java -Xmx4g -jar $PICARDJARPATH/picard.jar BedToIntervalList \
            I={input.target} \
            O={output.target} \
            SD={params.ref}
        """

# Picard HsMetrics
rule picard_hsmetrics:
    """
    Run Picard CollectHsMetrics on aligned reads
    """
    input:
        bam=join(RESULTSDIR,"bam","{sample}.mapped.bam")
    output:
        metrics=join(RESULTSDIR,"qc","picard","{sample}.hs_metrics")
    params:
        ref = config[genome]['ref'],
        bait = join(RESULTSDIR,"bait.interval"),
        target = join(RESULTSDIR,"target.interval")
    envmodules:
        TOOLS["picard"]
    shell:
        """
        java -Xmx4g -jar $PICARDJARPATH/picard.jar CollectHsMetrics \
            -I {input.bam} \
            -R {params.ref} \
            -O {output.metrics} \
            -BAIT_INTERVALS {params.bait} \
            -TARGET_INTERVALS {params.target}
        """

# MultiQC
rule multiqc:
    """
    Run MultiQC
    """
    input:
        fastqc=expand(rules.fastqc.output.html,sample=SAMPLES),
        umis=expand(join(RESULTSDIR,"qc","fgbio","{samples}.grouped-family-sizes.txt"),samples=SAMPLES),
        flagstat=expand(rules.samtools_stats.output.flagstat,sample=SAMPLES),
        idxstats=expand(rules.samtools_stats.output.idxstats,sample=SAMPLES),
        hsmetrics=expand(rules.picard_hsmetrics.output.metrics,sample=SAMPLES),
    output:
        report = join(RESULTSDIR,"qc","multiqc_report.html")
    params:
        qc_dir=join(RESULTSDIR,"qc")
    envmodules:
        TOOLS["multiqc"]
    shell:
        """
        multiqc -f -o {params.qc_dir} {params.qc_dir} 
        """
