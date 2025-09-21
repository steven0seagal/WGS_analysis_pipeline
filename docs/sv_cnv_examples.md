# SV/CNV Pipeline Examples

This document provides practical examples of using the SV/CNV analysis pipeline for different scenarios.

## Table of Contents

- [Basic Single Sample Analysis](#basic-single-sample-analysis)
- [Multi-Sample Cohort Analysis](#multi-sample-cohort-analysis)
- [Clinical Analysis with High Specificity](#clinical-analysis-with-high-specificity)
- [Research Discovery with High Sensitivity](#research-discovery-with-high-sensitivity)
- [Parameter Optimization Examples](#parameter-optimization-examples)
- [Cluster Execution Examples](#cluster-execution-examples)

## Basic Single Sample Analysis

### Example 1: Quick Analysis with Convenience Script

```bash
# Simple analysis of a single sample
./run_sv_cnv.sh NA12878 /data/fastq /data/results /ref/GRCh38.fa

# With custom thread count
./run_sv_cnv.sh NA12878 /data/fastq /data/results /ref/GRCh38.fa 16

# Using bash mode instead of snakemake
./run_sv_cnv.sh NA12878 /data/fastq /data/results /ref/GRCh38.fa 8 bash
```

### Example 2: Manual Bash Script Execution

```bash
# Navigate to SV/CNV workflow directory
cd workflows/sv_cnv

# Create configuration file
cat > my_analysis_config.sh << EOF
SAMPLE="patient_001"
THREADS=12
REF_GENOME="/data/reference/GRCh38.fasta"
ADAPTERS="/data/adapters/TruSeq3-PE.fa"
RAW_READS_DIR="/data/fastq/patient_001"
OUTPUT_DIR="/data/analysis/patient_001"
BIN_SIZE=1000
REF_CHROM_DIR="/data/reference/GRCh38_chromosomes"
EOF

# Run the pipeline
./sv_cnv_pipeline.sh my_analysis_config.sh
```

### Example 3: Manual Snakemake Execution

```bash
cd workflows/sv_cnv

# Edit config.yaml
cat > config.yaml << EOF
samples_tsv: "samples.tsv"
raw_reads_dir: "/data/fastq"
reference_fasta: "/data/reference/GRCh38.fasta"
adapters_fasta: "/data/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/data/reference/GRCh38_chromosomes"
cnvnator_bin_size: 1000
survivor_min_callers: 2
sv_filter_expr: 'SVTYPE!="BND" && QUAL>30'
cnv_filter_expr: 'INFO/SVLEN>5000 || INFO/SVLEN<-5000'
EOF

# Create sample list
echo -e "sample\npatient_001" > samples.tsv

# Run pipeline
snakemake --use-conda -j 12
```

## Multi-Sample Cohort Analysis

### Example 4: Family Trio Analysis

```bash
# Create samples file for family trio
cat > family_trio.tsv << EOF
sample
proband
mother
father
EOF

# Update config for trio analysis
cat > config_trio.yaml << EOF
samples_tsv: "family_trio.tsv"
raw_reads_dir: "/data/trio_fastq"
reference_fasta: "/data/reference/GRCh38.fasta"
adapters_fasta: "/data/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/data/reference/GRCh38_chromosomes"

# Trio-specific parameters
cnvnator_bin_size: 1000
survivor_min_callers: 2
sv_filter_expr: 'SVTYPE!="BND" && QUAL>20'  # Slightly relaxed for trio
cnv_filter_expr: 'INFO/SVLEN>1000 || INFO/SVLEN<-1000'  # Smaller CNVs for trio
EOF

# Run on cluster with more resources
snakemake --configfile config_trio.yaml --use-conda \
    --cluster "sbatch -n {threads} --mem={resources.mem_mb}M -t 24:00:00" \
    --jobs 50
```

### Example 5: Population Cohort (100 samples)

```bash
# Generate samples file
for i in $(seq 1 100); do
    echo "sample_$(printf "%03d" $i)"
done > cohort_100.tsv
echo -e "sample\n$(cat cohort_100.tsv)" > cohort_100.tsv

# High-throughput configuration
cat > config_cohort.yaml << EOF
samples_tsv: "cohort_100.tsv"
raw_reads_dir: "/data/cohort_fastq"
reference_fasta: "/data/reference/GRCh38.fasta"
adapters_fasta: "/data/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/data/reference/GRCh38_chromosomes"

# Population-level parameters
cnvnator_bin_size: 1000
survivor_min_callers: 2
sv_filter_expr: 'SVTYPE!="BND" && QUAL>30 && INFO/PE>3'
cnv_filter_expr: 'INFO/SVLEN>5000 || INFO/SVLEN<-5000'

# Resource allocation
threads:
  alignment: 8
  sv_calling: 4
  cnv_calling: 4
memory:
  alignment: 8000
  sv_calling: 4000
  cnv_calling: 8000
EOF

# Submit to SLURM cluster
snakemake --configfile config_cohort.yaml --use-conda \
    --cluster-config cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -n {cluster.threads} --mem={cluster.mem}G -t {cluster.time}" \
    --jobs 200
```

## Clinical Analysis with High Specificity

### Example 6: Clinical Diagnostic Pipeline

```bash
# High-specificity configuration for clinical use
cat > config_clinical.yaml << EOF
samples_tsv: "clinical_samples.tsv"
raw_reads_dir: "/clinical/fastq"
reference_fasta: "/ref/GRCh38.fasta"
adapters_fasta: "/ref/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/ref/GRCh38_chromosomes"

# Conservative parameters for clinical accuracy
cnvnator_bin_size: 1000
survivor_max_dist: 500              # Tighter breakpoint matching
survivor_min_callers: 3             # Require 3 callers for high confidence
survivor_min_size: 100              # Larger minimum size

# Stringent filtering
sv_filter_expr: 'SVTYPE!="BND" && QUAL>50 && INFO/PE>10 && INFO/SR>2'
cnv_filter_expr: 'INFO/SVLEN>10000 || INFO/SVLEN<-10000'  # Only large CNVs

# Annotation parameters for clinical interpretation
vep_extra_args: "--symbol --terms SO --hgvs --variant_class --sift b --polyphen b --clinical_significance"
EOF

# Run with extra validation
snakemake --configfile config_clinical.yaml --use-conda -j 8 \
    --additional-rule validate_vcfs
```

### Example 7: Rare Disease Analysis

```bash
# Focus on rare, potentially pathogenic variants
cat > config_rare_disease.yaml << EOF
samples_tsv: "patient.tsv"
raw_reads_dir: "/data/patient_fastq"
reference_fasta: "/ref/GRCh38.fasta"
adapters_fasta: "/ref/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/ref/GRCh38_chromosomes"

# Parameters optimized for rare disease discovery
cnvnator_bin_size: 500               # Higher resolution
survivor_min_callers: 2
survivor_min_size: 50

# Filter for rare variants
sv_filter_expr: 'SVTYPE!="BND" && QUAL>30 && INFO/PE>5'
cnv_filter_expr: 'INFO/SVLEN>1000 || INFO/SVLEN<-1000'

# Comprehensive annotation for interpretation
vep_extra_args: "--symbol --terms SO --hgvs --variant_class --regulatory --custom gnomad_sv,/db/gnomad_sv.vcf.gz,vcf,exact,0,AF"
EOF
```

## Research Discovery with High Sensitivity

### Example 8: High-Sensitivity Discovery

```bash
# Configuration for maximum sensitivity
cat > config_discovery.yaml << EOF
samples_tsv: "research_samples.tsv"
raw_reads_dir: "/research/fastq"
reference_fasta: "/ref/GRCh38.fasta"
adapters_fasta: "/ref/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/ref/GRCh38_chromosomes"

# Sensitive parameters
cnvnator_bin_size: 1000
survivor_min_callers: 1             # Include variants from single caller
survivor_min_size: 50

# Permissive filtering (filter downstream)
sv_filter_expr: 'SVTYPE!="BND" && QUAL>10'
cnv_filter_expr: 'INFO/SVLEN>500 || INFO/SVLEN<-500'
EOF

# Run with all optional tools enabled
snakemake --configfile config_discovery.yaml --use-conda -j 16
```

## Parameter Optimization Examples

### Example 9: Coverage-Dependent Parameters

```bash
# For high coverage samples (60x)
cat > config_high_coverage.yaml << EOF
# Smaller bins for higher resolution
cnvnator_bin_size: 500
gatk_gcnv_bin_length: 500

# Can afford stricter filters with high coverage
sv_filter_expr: 'SVTYPE!="BND" && QUAL>40 && INFO/PE>10'
cnv_filter_expr: 'INFO/SVLEN>2000 || INFO/SVLEN<-2000'
EOF

# For low coverage samples (15x)
cat > config_low_coverage.yaml << EOF
# Larger bins for stability
cnvnator_bin_size: 2000
gatk_gcnv_bin_length: 2000

# More permissive filters for low coverage
sv_filter_expr: 'SVTYPE!="BND" && QUAL>20 && INFO/PE>3'
cnv_filter_expr: 'INFO/SVLEN>10000 || INFO/SVLEN<-10000'
EOF
```

### Example 10: Tumor Sample Analysis

```bash
# Special considerations for tumor samples
cat > config_tumor.yaml << EOF
samples_tsv: "tumor_samples.tsv"
raw_reads_dir: "/tumor/fastq"
reference_fasta: "/ref/GRCh38.fasta"
adapters_fasta: "/ref/adapters/TruSeq3-PE.fa"
reference_chrom_dir: "/ref/GRCh38_chromosomes"

# Tumor-specific parameters
cnvnator_bin_size: 1000
survivor_min_callers: 2
survivor_min_size: 50

# Account for tumor heterogeneity and aneuploidy
sv_filter_expr: 'SVTYPE!="BND" && QUAL>25'  # Slightly relaxed
cnv_filter_expr: 'INFO/SVLEN>1000 || INFO/SVLEN<-1000'

# Additional tumor-specific options
trimmomatic_extra: "ILLUMINACLIP:{adapters}:2:30:10 SLIDINGWINDOW:4:15 MINLEN:50"  # More aggressive trimming
EOF
```

## Cluster Execution Examples

### Example 11: SLURM Cluster Configuration

```bash
# Create cluster configuration file
cat > cluster.yaml << EOF
__default__:
  partition: "general"
  threads: 4
  mem: 8
  time: "4:00:00"

bwa_map:
  threads: 16
  mem: 16
  time: "8:00:00"

manta:
  threads: 8
  mem: 16
  time: "6:00:00"

cnvnator:
  threads: 4
  mem: 32
  time: "4:00:00"

annotsv:
  threads: 2
  mem: 8
  time: "2:00:00"
EOF

# Submit jobs to cluster
snakemake --use-conda \
    --cluster-config cluster.yaml \
    --cluster "sbatch -p {cluster.partition} -n {cluster.threads} --mem={cluster.mem}G -t {cluster.time}" \
    --jobs 100 \
    --rerun-incomplete
```

### Example 12: PBS/Torque Cluster

```bash
# PBS cluster submission
snakemake --use-conda \
    --cluster "qsub -q batch -l nodes=1:ppn={threads},mem={resources.mem_mb}mb,walltime=24:00:00" \
    --jobs 50 \
    --jobscript cluster_jobscript.sh
```

### Example 13: Cloud Execution (AWS)

```bash
# Using Snakemake with Kubernetes on AWS
snakemake --use-conda \
    --kubernetes \
    --default-remote-provider S3 \
    --default-remote-prefix s3://my-bucket/sv-cnv-analysis \
    --jobs 100
```

## Quality Control and Validation Examples

### Example 14: QC-Focused Run

```bash
# Run only QC stages to assess data quality first
snakemake --use-conda -j 8 qc/multiqc_report/multiqc_report.html

# Check the results before proceeding
firefox qc/multiqc_report/multiqc_report.html

# If QC passes, run full pipeline
snakemake --use-conda -j 8
```

### Example 15: Incremental Analysis

```bash
# Run stages incrementally for large datasets
# Stage 1: Preprocessing only
snakemake --use-conda -j 8 \
    $(for sample in $(cat samples.tsv | tail -n +2); do
        echo "analysis_ready_bam/${sample}.dedup.bam"
      done)

# Stage 2: SV calling only
snakemake --use-conda -j 8 \
    $(for sample in $(cat samples.tsv | tail -n +2); do
        echo "sv_calls/manta/${sample}/results/variants/diploidSV.vcf.gz"
        echo "sv_calls/delly/${sample}.merged.vcf.gz"
        echo "sv_calls/lumpy/${sample}.vcf"
      done)

# Stage 3: Complete pipeline
snakemake --use-conda -j 8
```

These examples demonstrate the flexibility of the SV/CNV pipeline for various research and clinical applications. Adjust parameters based on your specific requirements, sample characteristics, and computational resources.