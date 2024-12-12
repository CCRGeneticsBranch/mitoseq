# Call variants

rule mutect2_mt:
    """Runs GATK MuTect2 in mitochondrial mode (single sample only)"""
    input:
        bam = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam")
    output:
        vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.filtered.vcf")
    params:
        ref = config[genome]['ref']
    envmodules:
        TOOLS["gatk4"]
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
        """

rule mutect2_paired:
    """Runs GATK MuTect2 on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda wc: join(RESULTSDIR,"bam","%s.cons.filtered.bam") % PAIRS[wc.sample]
    output:
        vcf = join(RESULTSDIR,"variants","mutect2","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.filtered.vcf")
    params:
        ref = config[genome]['ref'],
        normalsample = lambda wc: "%s" % PAIRS[wc.sample]
    envmodules:
        TOOLS["gatk4"]
    shell: 
        """
        gatk Mutect2 \
            -R {params.ref} \
            -L chrM \
            -I {input.tumor} \
            -I {input.normal} \
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
         """

rule varscan_paired:
    """Runs VarScan on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda wc: join(RESULTSDIR,"bam","%s.cons.filtered.bam") % PAIRS[wc.sample]
    output:
        vcf = join(RESULTSDIR,"variants","varscan","{sample}.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","varscan","{sample}.filtered.vcf")
    params:
        ref = config[genome]['ref'],
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
        """