# Call variants

rule mutect2_mt:
    """Runs GATK MuTect2 in mitochondrial mode (single sample only)"""
    input:
        bam = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam")
    output:
        vcf = join(RESULTSDIR,"variants","mutect2","{sample}.chrM.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.chrM.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.chrM.filtered.vcf")
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
        normal = lambda w: [join(RESULTSDIR,"bam",PAIRS[w.samples]".cons.filtered.bam")]
    output:
        vcf = join(RESULTSDIR,"variants","mutect2","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.filtered.vcf")
    params:
        ref = config[genome]['ref'],
        normalsample = lambda w: [PAIRS[w.samples]]
    envmodules:
        TOOLS["gatk4"]
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
         """

rule varscan_paired:
    """Runs VarScan on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda w: [join(RESULTSDIR,"bam",PAIRS[w.samples]".cons.filtered.bam")]
    output:
        vcf = join(RESULTSDIR,"variants","varscan","{samples}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","varscan","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","varscan","{sample}.filtered.vcf")
    params:
        ref = config[genome]['ref'],
        tumorsample = '{samples}',
        normalsample = lambda w: [PAIRS[w.samples]]
    threads: 
        getthreads("varscan_paired")
    envmodules:
        TOOLS["samtools"],
        TOOLS["varscan"],
        TOOLS["gatk3"]
    shell: 
        """
        samtools mpileup -r chrM -d 0 -f {params.ref} {input.normal} > {params.normalsample}.pileup
        samtools mpileup -r chrM -d 0 -f {params.ref} {input.tumor} > {params.tumorsample}.pileup

        varscan somatic {params.normalsample}.pileup {params.tumorsample}.pileup {output.vcf} \
         --strand-filter 1 --output-vcf 1 --min-var-freq 0.001 

        java8 -Xmx12g -XX:ParallelGCThreads={threads} \\
            -jar $GATK_JAR -T CombineVariants \\
            -R {params.ref} \\
            --variant {output.vcf}.snp \\
            --variant {output.vcf}.indel \\
            --assumeIdenticalSamples \\
            --filteredrecordsmergetype KEEP_UNCONDITIONAL \\
            -o {output.vcf}
        
        varscan somaticFilter {output.vcf} \
            --min-strands2 2 --min-var-freq 0.001 \
            --output-file {output.marked_vcf}

        gatk SelectVariants \
            -R {params.ref} \
            --variant {output.marked_vcf} \
            --exclude-filtered \
            --output {output.filtered_vcf}
        """