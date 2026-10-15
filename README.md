# MitoSeq Pipeline

Snakemake pipeline for variant calling from high-coverage mitochondrial DNA capture sequencing with unique molecular identifiers (UMIs). Pipeline builds consensus reads by UMI based on [fgbio best practices](https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md); aligns to the full genome, an NuMT-masked genome, or chrM only; calls mitochondrial variants with GATK Mutect2 (mitochondria mode and with tumor/normal pairs) and VarScan2 (tumor/normal pairs); and reports capture QC. Human and mouse samples must be processed in separate runs.

## Pipeline steps

| Step | Rule | Tool | Description |
|---|---|---|---|
| 1 | `copy_uBam` | Picard `AddOrReplaceReadGroups` | Copy the unmapped BAM and set sample name |
| 2 | `align_bam` | BWA-MEM, fgbio `ZipperBams`, samtools | Align raw reads and carry UMI tags over from the unmapped BAM |
| 3 | `group_reads` | fgbio `GroupReadsByUmi` | Group reads by UMI/position (Adjacency strategy, 1 edit) |
| 4 | `call_consensus_reads` | fgbio `CallMolecularConsensusReads` | Build single-strand consensus reads |
| 5 | `align_consensus_bam` | BWA-MEM, fgbio `ZipperBams`, samtools | Align consensus reads |
| 6 | `filter_consensus_reads` | fgbio `FilterConsensusReads`, samtools | Filter consensus reads and bases, sort |
| 7 | `fastqc`, `samtools_stats`, `create_interval_list`, `picard_hsmetrics`, `multiqc` | FastQC, samtools, Picard, MultiQC | Read QC, alignment statistics, and capture (hybrid selection) metrics |
| 8 | `mutect2_mt` | GATK4 Mutect2 (`--mitochondria-mode`), vcf2maf/VEP | chrM variants for every sample |
| 9 | `mutect2_paired` | GATK4 Mutect2, vcf2maf/VEP | chrM somatic variants for tumor/normal pairs |
| 10 | `varscan_paired` | samtools `mpileup`, VarScan2, GATK3, vcf2maf/VEP | chrM somatic variants for tumor/normal pairs |
| 11 | `combine_mafs` | | Merge MAFs across samples; extra Mutect2 filter (VAF ≥ 1%, no `normal_artifact`/`weak_evidence`) |

## Requirements

- Runs on NIH [Biowulf](https://hpc.nih.gov/) (SLURM and environment modules).
- Must have access to dependencies listed in [config/config.yaml](config/config.yaml).

## Installation

```bash
git clone https://github.com/CCRGeneticsBranch/umi_correction.git mitoseq
```

## Quick start

The pipeline is designed to be run in a separate directory (/path/to/my_run).

```bash
# 1. Create a work folder with template config files
bash mitoseq/run_mitoseq --runmode=init --workdir=/path/to/my_run [--samplesheet=/path/to/samples.tsv]

# 2. Edit /path/to/my_run/config/config.yaml and config/samples.tsv (see Inputs)

# 3. Dry run (writes /path/to/my_run/dryrun.log)
bash mitoseq/run_mitoseq --runmode=dryrun --workdir=/path/to/my_run

# 4. Submit to SLURM
bash mitoseq/run_mitoseq --runmode=run --workdir=/path/to/my_run

# Unlock the work folder after an interrupted run
bash mitoseq/run_mitoseq --runmode=unlock --workdir=/path/to/my_run
```

`-h` shows all options; `-v` prints the pipeline and launcher versions. Per-rule SLURM resources are in `WORKDIR/config/cluster.yaml` (copied from [resources/cluster.yaml](resources/cluster.yaml) at `init`).

## Inputs

### Sample sheet (`config/samples.tsv`)

| Column | Description | Example |
|---|---|---|
| `Sample` | Sample name | `tumor1` |
| `uBam` | Absolute path to the unmapped BAM, with UMIs in the `RX` tag (e.g., from fgbio `FastqToBam` or demultiplexing to BAM) | `/path/to/tumor1.unmapped.bam` |
| `Normal` | For tumor samples, the `Sample` name of the matched normal; leave empty otherwise | `normal1` |

### Configuration (`config/config.yaml`)

| Key | Default | Description (source) | Edit? |
|---|---|---|---|
| `genome` | `hg38` | Unmasked genome (`hg38`, `mm10`), NuMT-masked genome (`hg38_masked`, `mm10_masked`), chrM only (`hg38_chrM`, `mm10_chrM`) | As needed |
| `hg38`, `mm10`: `ref`, `bwa_ref` | `/fdb/genomebrowser/fasta/`, `/fdb/bwa/indexes/` | UCSC reference FASTA and BWA indexes | No |
| `hg38_masked`, `mm10_masked`: `ref`, `bwa_ref` | `/data/khanlab3/gb_omics/resources/mitoblacklist/masked_genomes/` | NuMT-masked genome FASTA and BWA indexes (UCSC reference genome masked using [mitoblacklist](https://github.com/caleblareau/mitoblacklist) and bedtools makefasta) | No |
| `hg38_chrM`, `mm10_chrM`: `ref`, `bwa_ref` | `/data/khanlab3/gb_omics/resources/mitoblacklist/chrM_genome/` | chrM-only FASTA and BWA indexes (extracted from UCSC genomes) | No |
| `<genome>: bait`, `target` | `PIPELINE_HOME/resources/xgen-*.bed` | xGen panel probe and target BEDs for Picard HsMetrics (IDT xGen; exact product/catalog source and redistribution permission: **TODO**) | No |
| `<genome>: vep` | e.g. `--ncbi-build GRCh38 --species homo_sapiens` | vcf2maf/VEP species and build | No |
| `min_base_qual` | `30` | Minimum consensus base quality (`FilterConsensusReads`) | As needed |
| `max_reads_mutect2` | `0` | Mutect2 `--max-reads-per-alignment-start` (0 = no downsampling) | As needed |
| `workdir`, `manifest`, `tools`, `cluster` | `WORKDIR/...`, `PIPELINE_HOME/...` | Work folder and launcher/repository paths (filled by `init`) | No |

## Outputs

All paths are relative to the work folder.

| Path | Description |
|---|---|
| `results/bam/<sample>.cons.filtered.bam` | Filtered, aligned consensus reads |
| `results/qc/fastqc/`, `results/qc/samtools/`, `results/qc/picard/` | FastQC, flagstat/idxstats, and HsMetrics per sample |
| `results/qc/multiqc_report.html` | Combined QC report |
| `results/variants/mutect2_chrM/<sample>.filtered.vcf`, `.maf` | Mutect2 mitochondria-mode calls per sample |
| `results/variants/mutect2/<tumor>.filtered.vcf`, `.filtered.maf`, `.marked.maf` | Mutect2 tumor/normal calls (all and filtered) |
| `results/variants/varscan/<tumor>.filtered.Somatic.vcf`, `.maf` | VarScan2 somatic calls |
| `results/variants/all.*.maf` | Merged MAFs across samples (Mutect2 mitochondria-mode: `all.chrM.maf`, Mutect2 tumor/normal, standard filtering: `all.mutect2.filtered.maf`, Mutect2 tumor/normal, no filtering: `all.mutect2.marked.maf`, Mutect2 tumor/normal, custom filtering: `all.mutect2.filteredv2.maf`, Varscan: `all.varscan.maf`) |
| `results/pipeline_version.txt` | Pipeline version used for the run |
| `runinfo.yaml`, `snakemake.log`, `logs/` | Run information, Snakemake log and per-rule SLURM logs |

## Tools and references

Modules are listed in [resources/tools.yaml](resources/tools.yaml).

| Tool | Module | Used in | Reference |
|---|---|---|---|
| Snakemake | `snakemake/7.19.1` | workflow | Mölder et al. 2021, *F1000Research* [doi:10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2) |
| fgbio | `fgbio` | steps 2–6 | [fulcrumgenomics/fgbio](https://github.com/fulcrumgenomics/fgbio) |
| BWA-MEM | `bwa` | `align_bam`, `align_consensus_bam` | Li 2013, [arXiv:1303.3997](https://arxiv.org/abs/1303.3997) |
| SAMtools | `samtools` | alignment, QC, VarScan input | Danecek et al. 2021, *GigaScience* [doi:10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008) |
| Picard | `picard` | `copy_uBam`, `picard_hsmetrics` | [Broad Institute Picard](https://broadinstitute.github.io/picard/) |
| FastQC | `fastqc` | `fastqc` | [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| MultiQC | `multiqc` | `multiqc` | Ewels et al. 2016, *Bioinformatics* [doi:10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) |
| GATK4 Mutect2 | `GATK/4.4.0.0` | `mutect2_mt`, `mutect2_paired` | Benjamin et al. 2019, *bioRxiv* [doi:10.1101/861054](https://doi.org/10.1101/861054) |
| GATK3 | `GATK/3.8-1` | `varscan_paired` | McKenna et al. 2010, *Genome Res* [doi:10.1101/gr.107524.110](https://doi.org/10.1101/gr.107524.110) |
| VarScan2 | `VarScan` | `varscan_paired` | Koboldt et al. 2012, *Genome Res* [doi:10.1101/gr.129684.111](https://doi.org/10.1101/gr.129684.111) |
| vcf2maf | `vcf2maf` | MAF conversion | [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf) |
| Ensembl VEP | `VEP/102` | MAF annotation | McLaren et al. 2016, *Genome Biol* [doi:10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4) |
| CCBR Xavier |  | Reference for variant calling | https://zenodo.org/records/18421720 | 

## Versioning

The pipeline version is in [VERSION](VERSION) and changes are listed in [CHANGELOG.md](CHANGELOG.md). Each run records the pipeline version, git commit and launcher version in `runinfo.yaml` and `results/pipeline_version.txt`. The launch script is a copy of [biowulf-snakemake-launcher](https://github.com/CCRGeneticsBranch/biowulf-snakemake-launcher) (`bin/launcher.sh`).

## License

Original pipeline code and documentation are licensed under MIT; see [LICENSE](LICENSE). Reference data and bundled xGen panel files retain their own terms. The MIT license does not grant redistribution rights for external or bundled third-party resources.

## Contact

Erica Pehrsson, Genetics Branch, Center for Cancer Research, National Cancer Institute. Please report problems through GitHub issues.
