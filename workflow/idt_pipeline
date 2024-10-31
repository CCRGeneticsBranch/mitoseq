#!/bin/bash

# Implements the xGen Prism DNA Library Prep Kit suggestions for 
# processing sequence data with UMIs (2020)

ml java
ml picard/2.26.9
ml bwa
ml vardict/1.5.8

# Set paths
WORK_DIR="/data/khanlab3/Erica/cao/targeted_umi/"
SAMPLE=$1
FGBIO_PATH="/data/pehrssonec/fgbio-0.7.0.jar"
REF="/data/khanlab/projects/ngs_pipeline_testing/References_4.0/New_GRCh37/Index/bwa_0.7.17/hg19.fa" #"/data/khanlab/projects/ngs_pipeline_testing/References_4.0/GRCh38/Index/bwa_0.7.17/hg38.fa"
TARGETS="xgen_aml_targets.bed"
VARDICT="/usr/local/apps/vardict/1.5.8/bin"

# Construct unmapped bam
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR FastqToSam \
  FASTQ=$WORK_DIR$SAMPLE\.R1.fastq.gz \
  FASTQ2=$WORK_DIR$SAMPLE\.R2.fastq.gz \
  O=$WORK_DIR$SAMPLE\.unmapped.bam \
  SM=$SAMPLE

# Extract UMIs and add to the RX tag
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH ExtractUmisFromBam \
    --input=$WORK_DIR$SAMPLE\.unmapped.bam --output=$WORK_DIR$SAMPLE\.unmapped.withUMI.bam \
    --read-structure=8M143T 8M143T --molecular-index-tags=ZA ZB --single-tag=RX

# Align reads and merge unmapped/mapped bams
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR SamToFastq \
    I=$WORK_DIR$SAMPLE\.unmapped.withUMI.bam F=/dev/stdout INTERLEAVE=true \
  | bwa mem -p -t 4 $REF /dev/stdin \
  | java -Xmx24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR MergeBamAlignment \
    UNMAPPED=$WORK_DIR$SAMPLE\.unmapped.withUMI.bam ALIGNED=/dev/stdin O=$WORK_DIR$SAMPLE\.mapped.bam R=$REF \
    SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
    ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

# Error correction by collapsing single read families (same start-stop, UMI)

## Correct errors in the UMI sequences
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH CorrectUmis \
  -i $WORK_DIR$SAMPLE\.mapped.bam -o $WORK_DIR$SAMPLE\.mapped.fixedumi.bam \
  --max-mismatches=3 \
  --min-distance=1 \
  -M $WORK_DIR$SAMPLE\.metrics.txt \
  -r $WORK_DIR$SAMPLE\.rejected.bam -t RX -u GAGACGAT TTCCAAGG CGCATGAT ACGGAACA CGGCTAAT GCTATCCT TGGACTCT \
    ATCCAGAG CTTAGGAC GTGCCATA TCGCTGTT TTCGTTGG AAGCACTG GTCGAAGA ACCACGAT GATTACCG GCACAACT \
    GCGTCATT GAAGGAAG ACTGAGGT TGAAGACG GTTACGCA AGCGTGTT GATCGAGT TTGCGAAG CTGTTGAC GATGTGTG \
    ACGTTCAG TTGCAGAC CAATGTGG ACGACTTG ACTAGGAG

## Group read families (assign source to each read, store in MI tag, output bam sorted by MI tag)
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH GroupReadsByUmi \
  --input=$WORK_DIR$SAMPLE\.mapped.fixedumi.bam --output=$WORK_DIR$SAMPLE\.grouped.bam \
  --strategy=paired --edits=0 --min-map-q=20

## Collapse single read families
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH CallMolecularConsensusReads \
  --input=$WORK_DIR$SAMPLE\.grouped.bam --output=$WORK_DIR$SAMPLE\.consensus.unmapped.bam \
  --error-rate-pre-umi=45 --error-rate-post-umi=30 \
  --min-consensus-base-quality=40 --min-input-base-quality=30 --min-reads=1

## Align collapsed families and merge mapped/unmapped bams
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR SamToFastq \
  I=$WORK_DIR$SAMPLE\.consensus.unmapped.bam F=/dev/stdout INTERLEAVE=true \
  | bwa mem -p -t 4 $REF /dev/stdin \
  | java -Xmx24g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $PICARD_JAR MergeBamAlignment \
    UNMAPPED=$WORK_DIR$SAMPLE\.consensus.unmapped.bam ALIGNED=/dev/stdin \
    O=$WORK_DIR$SAMPLE\.consensus.mapped.bam R=$REF \
    SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
    ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

## Filter collapsed families
java -Xmx8g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH FilterConsensusReads \
  --input=$WORK_DIR$SAMPLE\.consensus.mapped.bam \
  --output=$WORK_DIR$SAMPLE\.consensus.mapped.filtered.bam \
  --ref=$REF \
  --min-reads=3 \
  --max-read-error-rate=0.05 \
  --max-base-error-rate=0.1 \
  --min-base-quality=40 \
  --max-no-call-fraction=0.1

# Clip overlap between read pairs
java -Xmx4g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID -jar $FGBIO_PATH ClipBam \
  --input=$WORK_DIR$SAMPLE\.consensus.mapped.filtered.bam \
  --output=$WORK_DIR$SAMPLE\.consensus.mapped.filtered.clipped.bam \
  --ref=$REF --soft-clip=false --clip-overlapping-reads=true

# Call variants using VarDict
min_af=0.01

$VARDICT/VarDict \
  -G $REF \
  -N $SAMPLE \
  -f ${min_af} \
  -b $WORK_DIR$SAMPLE\.consensus.mapped.filtered.clipped.bam \
  -z -c 1 -S 2 -E 3 -g 4 -th 1 --nosv \
  $WORK_DIR$TARGETS \
| $VARDICT/teststrandbias.R \
| $VARDICT/var2vcf_valid.pl -N $SAMPLE -E -f ${min_af} \
| awk '{if ($1 ~ /^#/) print; else if ($4 !=5) print}' \
> $SAMPLE\.vcf
