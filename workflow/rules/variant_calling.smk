
# Call variants
localrules: combine_mafs

rule mutect2_mt:
    """Runs GATK MuTect2 in mitochondrial mode (single sample only)"""
    input:
        bam = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam")
    output:
        vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.filtered.vcf"),
        maf = join(RESULTSDIR,"variants","mutect2_chrM","{sample}.maf")
    params:
        ref = config[genome]['ref'],
        tumorsample = "{sample}",
        vep = config[genome]['vep']
    envmodules:
        TOOLS["gatk4"],
        TOOLS["samtools"],
        TOOLS["VEP"],
        TOOLS["vcf2maf"]
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

        vcf2maf.pl --input-vcf {output.filtered_vcf} --output-maf {output.maf} --tumor-id {params.tumorsample} \
        --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta {params.ref} \
        {params.vep} --vep-forks $SLURM_CPUS_ON_NODE
        """

rule mutect2_paired:
    """Runs GATK MuTect2 on paired tumor/normal (adapted from Xavier)"""
    input:
        tumor = join(RESULTSDIR,"bam","{sample}.cons.filtered.bam"),
        normal = lambda wc: join(RESULTSDIR,"bam","%s.cons.filtered.bam") % PAIRS[wc.sample]
    output:
        vcf = join(RESULTSDIR,"variants","mutect2","{sample}.vcf"),
        marked_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.marked.vcf"),
        filtered_vcf = join(RESULTSDIR,"variants","mutect2","{sample}.filtered.vcf"),
        marked_maf = join(RESULTSDIR,"variants","mutect2","{sample}.marked.maf"),
        maf = join(RESULTSDIR,"variants","mutect2","{sample}.maf")
    params:
        ref = config[genome]['ref'],
        tumorsample = "{sample}",
        normalsample = lambda wc: "%s" % PAIRS[wc.sample],
        vep = config[genome]['vep']
    envmodules:
        TOOLS["gatk4"],
        TOOLS["samtools"],
        TOOLS["VEP"],
        TOOLS["vcf2maf"]
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

        # Convert .marked.vcf to .maf (for all calls, not just filtered)
        vcf2maf.pl --input-vcf {output.marked_vcf} --output-maf {output.marked_maf} --tumor-id {params.tumorsample} --normal-id {params.normalsample} \
        --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta {params.ref} \
        {params.vep} --vep-forks $SLURM_CPUS_ON_NODE

        vcf2maf.pl --input-vcf {output.filtered_vcf} --output-maf {output.maf} --tumor-id {params.tumorsample} --normal-id {params.normalsample} \
        --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta {params.ref} \
        {params.vep} --vep-forks $SLURM_CPUS_ON_NODE
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
        reheader = join(RESULTSDIR,"variants","varscan","{sample}.reheader.Somatic.vcf"),
        maf = join(RESULTSDIR,"variants","varscan","{sample}.maf")
    params:
        ref = config[genome]['ref'],
        tumorsample = "{sample}",
        normalsample = lambda wc: "%s" % PAIRS[wc.sample],
        sample_list = join(RESULTSDIR,"variants","varscan","{sample}.list"),
        tumor_pileup = join(RESULTSDIR,"variants","varscan","{sample}.pileup"),
        normal_pileup = lambda wc: join(RESULTSDIR,"variants","varscan","%s.pileup") % PAIRS[wc.sample],
        vep = config[genome]['vep']
    threads: 
        getthreads("varscan_paired")
    envmodules:
        TOOLS["samtools"],
        TOOLS["varscan"],
        TOOLS["gatk3"],
        TOOLS["VEP"]
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

        bcftools reheader {output.somatic_vcf} -s {params.sample_list} > {output.reheader}

        ml vcf2maf
        vcf2maf.pl --input-vcf {output.reheader} --output-maf {output.maf} --tumor-id {params.tumorsample} --normal-id {params.normalsample} \
        --vep-path $VEP_HOME/bin --vep-data $VEP_CACHEDIR --cache-version $VEP_VERSION --ref-fasta {params.ref} \
        {params.vep} --vep-forks $SLURM_CPUS_ON_NODE
        """
      
rule combine_mafs:
    """Combine mafs across samples"""
    input:
        chrM = expand(join(RESULTSDIR,"variants","mutect2_chrM","{samples}.maf"),samples=SAMPLES),
        mutect2 = expand(join(RESULTSDIR,"variants","mutect2","{samples}.maf"),samples=PAIRS.keys()),
        mutect2_marked = expand(join(RESULTSDIR,"variants","mutect2","{samples}.marked.maf"),samples=PAIRS.keys()),
        varscan = expand(join(RESULTSDIR,"variants","varscan","{samples}.maf"),samples=PAIRS.keys())
    output:
        chrM = join(RESULTSDIR,"variants","all.chrM.maf"),
        mutect2 = join(RESULTSDIR,"variants","all.mutect2.maf"),
        mutect2_marked = join(RESULTSDIR,"variants","all.mutect2.marked.maf"),
        varscan = join(RESULTSDIR,"variants","all.varscan.maf")
    shell:
        """
        head -2 {input.chrM[0]} > {output.chrM}
        awk 'FNR>2 {{print}}' {input.chrM} >> {output.chrM}

        head -2 {input.mutect2[0]} > {output.mutect2}
        awk 'FNR>2 {{print}}' {input.mutect2} >> {output.mutect2}

        head -2 {input.mutect2_marked[0]} > {output.mutect2_marked}
        awk 'FNR>2 {{print}}' {input.mutect2_marked} >> {output.mutect2_marked}

        head -2 {input.varscan[0]} > {output.varscan}
        awk 'FNR>2 {{print}}' {input.varscan} >> {output.varscan}
        """
