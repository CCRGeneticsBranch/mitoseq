# Call variants

rule mutect2_mt:
    """Runs GATK MuTect2 in mitochondrial mode (single sample only)"""
    input:
        bam = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam")
    output:
        vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.filtered.vcf.gz"),
        header = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.header"),
        reheader = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.reheader.vcf.gz")
    params:
        ref = config[genome]['ref']
    envmodules:
        TOOLS["gatk4"],
        TOOLS["samtools"]
    shell: 
        """
        gatk Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -L chrM \
            --mitochondria-mode

        gatk FilterMutectCalls \
            -R {params.ref} \
            -V {output.vcf} \
            -O {output.marked_vcf} \
            --mitochondria-mode

        gatk SelectVariants \
            -R {params.ref} \
            --variant {output.marked_vcf} \
            --exclude-filtered \
            --output {output.filtered_vcf}

        bcftools head {output.filtered_vcf} | sed 's/AS_FilterStatus,Number=A/AS_FilterStatus,Number=./' > {output.header}
        bcftools reheader -h {output.header} {output.filtered_vcf} > {output.reheader}
        tabix {output.reheader}
        """

rule mutect2_paired:
    """Runs GATK MuTect2 on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda wc: join(RESULTSDIR,"bam","%s.cons.filtered.bam") % PAIRS[wc.sample]
    output:
        vcf = join(RESULTSDIR,"variants","mutect2","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.filtered.vcf.gz"),
        header = join(RESULTSDIR,"variants","mutect2","{sample}.header"),
        reheader = join(RESULTSDIR,"variants","mutect2","{sample}.reheader.vcf.gz")
    params:
        ref = config[genome]['ref'],
        normalsample = lambda wc: "%s" % PAIRS[wc.sample]
    envmodules:
        TOOLS["gatk4"],
        TOOLS["samtools"]
    shell: 
        """
        gatk Mutect2 \
            -R {params.ref} \
            -L chrM \
            -I {input.tumor} \
            -I {input.normal} \
            -normal {params.normalsample} \
            -O {output.vcf}

        gatk FilterMutectCalls \
            -R {params.ref} \
            -V {output.vcf} \
            -O {output.marked_vcf}

        gatk SelectVariants \
            -R {params.ref} \
            --variant {output.marked_vcf} \
            --exclude-filtered \
            --output {output.filtered_vcf}

        bcftools head {output.filtered_vcf} | sed 's/AS_FilterStatus,Number=A/AS_FilterStatus,Number=./' > {output.header}
        bcftools reheader -h {output.header} {output.filtered_vcf} > {output.reheader}
        tabix {output.reheader}
        """

rule varscan_paired:
    """Runs VarScan on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda wc: join(RESULTSDIR,"bam","%s.cons.filtered.bam") % PAIRS[wc.sample]
    output:
        vcf = join(RESULTSDIR,"variants","varscan","{sample}.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","varscan","{sample}.filtered.vcf"),
        somatic_vcf = join(RESULTSDIR,"variants","varscan","{sample}.filtered.Somatic.vcf"),
        reheader = join(RESULTSDIR,"variants","varscan","{sample}.reheader.Somatic.vcf.gz")
    params:
        ref = config[genome]['ref'],
        tumorsample = "{sample}",
        normalsample = lambda wc: "%s" % PAIRS[wc.sample],
        sample_list = join(RESULTSDIR,"variants","varscan","{sample}.list"),
        tumor_pileup = join(RESULTSDIR,"variants","varscan","{sample}.pileup"),
        normal_pileup = lambda wc: join(RESULTSDIR,"variants","varscan","%s.pileup") % PAIRS[wc.sample]
    threads: 
        getthreads("varscan_paired")
    envmodules:
        TOOLS["samtools"],
        TOOLS["varscan"],
        TOOLS["gatk3"]
    shell: 
        """
        samtools mpileup -r chrM -d 0 -f {params.ref} {input.normal} > {params.normal_pileup}
        samtools mpileup -r chrM -d 0 -f {params.ref} {input.tumor} > {params.tumor_pileup}

        varscan somatic {params.normal_pileup} {params.tumor_pileup} {output.vcf} \
         --strand-filter 1 --output-vcf 1 --min-var-freq 0.001 

        java -Xmx3g -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $GATK_JAR -T CombineVariants \
            -R {params.ref} \
            --variant {output.vcf}.snp \
            --variant {output.vcf}.indel \
            --assumeIdenticalSamples \
            --filteredrecordsmergetype KEEP_UNCONDITIONAL \
            -o {output.vcf}
        
        varscan somaticFilter {output.vcf} \
            --min-strands2 2 --min-var-freq 0.001 \
            --min-reads2 2 \
            --output-file {output.filtered_vcf}

        varscan processSomatic {output.filtered_vcf} \
            --min-tumor-freq 0.001 

        echo {params.normalsample} > {params.sample_list}
        echo {params.tumorsample} >> {params.sample_list}

        bgzip < <(bcftools reheader {output.somatic_vcf} -s {params.sample_list}) > {output.reheader}
        tabix {output.reheader}
        """

rule combine_variants:
    input:
        chrM = expand(join(RESULTSDIR,"variants","mutect2_chrM","{samples}.reheader.vcf.gz"),samples=SAMPLES),
        mutect2 = expand(join(RESULTSDIR,"variants","mutect2","{samples}.reheader.vcf.gz"),samples=PAIRS.keys()),
        varscan = expand(join(RESULTSDIR,"variants","varscan","{samples}.reheader.Somatic.vcf.gz"),samples=PAIRS.keys())
    output:
        chrM = join(RESULTSDIR,"variants","all.chrM.filtered.vcf.gz"),
        mutect2 = join(RESULTSDIR,"variants","all.mutect2.filtered.vcf.gz"),
        varscan = join(RESULTSDIR,"variants","all.varscan.filtered.vcf.gz")
    envmodules:
        TOOLS["samtools"]
    shell:
        """
        bcftools merge --force-samples {input.chrM} -o {output.chrM}
        bcftools merge --force-samples {input.mutect2} -o {output.mutect2}
        bcftools merge --force-samples {input.varscan} -o {output.varscan}
        """
