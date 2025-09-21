# Structural and Copy Number Variant Analysis

This document provides detailed information about the SV/CNV analysis pipeline implemented in this repository.

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Getting Started](#getting-started)
- [Configuration](#configuration)
- [Understanding the Output](#understanding-the-output)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

## Overview

The SV/CNV pipeline implements a comprehensive analysis framework for detecting structural variants (SVs) and copy number variants (CNVs) from whole genome sequencing data. The pipeline follows a six-stage architecture based on established best practices:

### Key Features

- **Ensemble calling strategy**: Combines multiple algorithms to maximize sensitivity and specificity
- **Dual modality**: Available as both bash script and Snakemake workflow
- **Comprehensive QC**: Integrated quality control at multiple stages
- **Robust filtering**: Configurable filters to balance sensitivity vs. specificity
- **Rich annotation**: Functional annotation with AnnotSV and VEP

### Supported Variant Types

- **Structural Variants**: Deletions, duplications, insertions, inversions, translocations
- **Copy Number Variants**: Deletions and duplications affecting gene dosage

## Pipeline Architecture

### Six-Stage Framework

1. **Data Ingestion & Quality Assurance**
   - FastQC quality assessment
   - Adapter trimming with Trimmomatic
   - Post-trimming quality validation

2. **Alignment & Preparation**
   - BWA-MEM alignment to reference genome
   - BAM sorting and indexing
   - PCR duplicate marking with Picard

3. **Structural Variant Calling (Ensemble)**
   - **Manta**: High-precision breakpoint detection via local assembly
   - **Delly**: Integrated paired-end and split-read calling
   - **Lumpy**: Probabilistic framework with multiple evidence types

4. **Copy Number Variant Calling**
   - **CNVnator**: Read-depth binning with mean-shift algorithm
   - **GATK gCNV**: Bayesian HMM approach (optional)

5. **Integration & Refinement**
   - SURVIVOR merging of SV callsets
   - BCFtools filtering for quality and size thresholds
   - Confidence-based variant selection

6. **Functional Annotation**
   - **AnnotSV**: ACMG-compliant pathogenicity classification
   - **VEP**: Molecular consequence prediction

### Detection Signatures

The pipeline leverages four orthogonal detection signatures:

| Signature | Description | Used by | Detects |
|-----------|-------------|---------|---------|
| Paired-End (PE) | Abnormal insert size/orientation | Manta, Delly, Lumpy | All SV types |
| Split-Read (SR) | Reads spanning breakpoints | Manta, Delly, Lumpy | Precise breakpoints |
| Read-Depth (RD) | Coverage fluctuations | CNVnator, GATK gCNV | CNVs |
| De novo Assembly | Local contig assembly | Manta | Complex rearrangements |

## Getting Started

### Prerequisites

- Conda/Mamba package manager
- Reference genome (GRCh38 recommended)
- WGS FASTQ files (>20x coverage recommended)
- At least 32GB RAM and 8 CPU cores

### Quick Start

#### Option 1: Bash Script

```bash
# Navigate to the SV/CNV workflow directory
cd workflows/sv_cnv

# Copy and edit the configuration file
cp config_example.sh my_config.sh
# Edit my_config.sh with your file paths

# Run the pipeline
./sv_cnv_pipeline.sh my_config.sh
```

#### Option 2: Snakemake Workflow

```bash
# Navigate to the SV/CNV workflow directory
cd workflows/sv_cnv

# Edit the configuration file
vim config.yaml

# Create sample list
echo -e "sample\nNA12878" > samples.tsv

# Run with conda environments
snakemake --use-conda -j 8

# Or run on cluster
snakemake --cluster "sbatch -n {threads}" --jobs 100 --use-conda
```

## Configuration

### Bash Script Configuration

The bash script uses a shell configuration file with the following required parameters:

```bash
# Required parameters
SAMPLE="sample_name"           # Sample identifier
THREADS=8                      # Number of CPU threads
REF_GENOME="/path/to/ref.fa"   # Reference genome FASTA
ADAPTERS="/path/to/adapters.fa"# Adapter sequences
RAW_READS_DIR="/path/to/fastq" # Input FASTQ directory
OUTPUT_DIR="/path/to/output"   # Output directory

# Optional parameters
BIN_SIZE=1000                  # CNV analysis bin size
CNVNATOR_PATH="cnvnator"       # Path to CNVnator
ANNOTSV_ROOT="AnnotSV"         # Path to AnnotSV
```

### Snakemake Configuration

The Snakemake workflow uses a YAML configuration file:

```yaml
# Input data
samples_tsv: "samples.tsv"
raw_reads_dir: "/path/to/fastq"

# Reference files
reference_fasta: "/path/to/reference.fa"
adapters_fasta: "/path/to/adapters.fa"
reference_chrom_dir: "/path/to/chromosomes"

# Tool parameters
cnvnator_bin_size: 1000
survivor_min_callers: 2
sv_filter_expr: 'SVTYPE!="BND" && QUAL>30'
```

### Key Parameters

#### Coverage-dependent Parameters

| Coverage | CNVnator bin_size | GATK gCNV bin_length |
|----------|-------------------|---------------------|
| 20-30x   | 1000             | 1000               |
| 40-60x   | 500              | 500                |
| >80x     | 200              | 200                |

#### SURVIVOR Merge Parameters

- `max_dist`: Maximum distance between breakpoints (default: 1000bp)
- `min_callers`: Minimum supporting callers (1=sensitive, 3=specific, 2=balanced)
- `min_size`: Minimum SV size (default: 50bp)

## Understanding the Output

### Directory Structure

```
output/
├── qc/                          # Quality control reports
│   ├── raw/                     # Pre-trimming FastQC
│   ├── trimmed/                 # Post-trimming FastQC
│   └── multiqc_report/          # Aggregated QC report
├── analysis_ready_bam/          # Final processed BAM files
├── sv_calls/                    # Individual SV caller outputs
│   ├── manta/
│   ├── delly/
│   └── lumpy/
├── cnv_calls/                   # CNV caller outputs
│   ├── cnvnator/
│   └── gatk_gcnv/
├── merged_filtered_sv/          # Final variant callsets
│   ├── {sample}.sv.merged.filtered.vcf.gz
│   └── {sample}.cnv.filtered.vcf.gz
├── annotations/                 # Annotated variants
│   ├── {sample}.sv.annotsv.tsv
│   └── {sample}.sv.vep.vcf
└── summary_{sample}.txt         # Analysis summary
```

### Key Output Files

#### Primary Variants

- `{sample}.sv.merged.filtered.vcf.gz`: High-confidence structural variants
- `{sample}.cnv.filtered.vcf.gz`: Copy number variants
- `{sample}.sv.annotsv.tsv`: Functionally annotated SVs with pathogenicity scores

#### Quality Control

- `multiqc_report.html`: Comprehensive QC dashboard
- `{sample}.dedup_metrics.txt`: PCR duplicate statistics
- `summary_{sample}.txt`: Analysis overview and variant counts

### VCF Format Details

#### Important INFO Fields

- `SVTYPE`: Variant type (DEL, DUP, INV, INS, BND)
- `SVLEN`: Variant length (negative for deletions)
- `END`: End coordinate of the variant
- `PE`: Number of supporting paired-end reads
- `SR`: Number of supporting split reads
- `SUPP`: Number of supporting callers
- `SUPP_VEC`: Binary vector of supporting callers

#### AnnotSV Output

- `AnnotSV_type`: Full or split annotation
- `Gene_name`: Affected gene(s)
- `AnnotSV_ranking`: Pathogenicity score (1-5)
- `ACMG_class`: ACMG classification (1-5, benign to pathogenic)

## Best Practices

### Input Requirements

1. **Coverage**: Minimum 20x coverage for reliable SV detection
2. **Insert size**: Standard Illumina libraries (300-500bp insert) work best
3. **Read length**: 150bp paired-end reads recommended

### Parameter Tuning

1. **Conservative approach**: Use `min_callers=2` for initial analysis
2. **High sensitivity**: Use `min_callers=1` for discovery studies
3. **High specificity**: Use `min_callers=3` for clinical applications

### Quality Filters

1. **Size filters**: Remove very small variants (`SVLEN > 50bp`)
2. **Quality filters**: Filter low-quality calls (`QUAL > 30`)
3. **Support filters**: Require multiple supporting reads (`PE > 5`)

### Manual Validation

Always perform manual validation of candidate variants using IGV:

1. Load the analysis-ready BAM file
2. Navigate to variant coordinates
3. Look for expected signatures:
   - **Deletions**: Reduced coverage, larger insert sizes
   - **Duplications**: Increased coverage, smaller insert sizes
   - **Inversions**: Discordant read orientations
   - **Translocations**: Inter-chromosomal read pairs

## Troubleshooting

### Common Issues

#### Low Variant Counts

**Symptoms**: Fewer variants than expected
**Causes**:
- Low coverage data
- Overly stringent filters
- Reference genome mismatch

**Solutions**:
- Check coverage with `samtools depth`
- Relax filtering parameters
- Verify reference genome version

#### High False Positive Rate

**Symptoms**: Many low-quality variants
**Causes**:
- Poor quality data
- Permissive filters
- Repetitive regions

**Solutions**:
- Review MultiQC report
- Increase quality thresholds
- Apply blacklist filters

#### Tool-specific Issues

#### CNVnator Errors

```bash
# Check ROOT installation
root-config --version

# Verify chromosome directory structure
ls $REF_CHROM_DIR
```

#### Manta Configuration Issues

```bash
# Check reference indexing
ls reference.fa.fai

# Verify BAM indexing
samtools quickcheck analysis_ready_bam/*.bam
```

### Performance Optimization

#### Memory Requirements

- **Manta**: 8-16GB RAM
- **Delly**: 4-8GB RAM
- **CNVnator**: 16-32GB RAM
- **GATK gCNV**: 8-16GB RAM

#### Runtime Estimates

For 30x WGS human sample:

| Stage | Runtime (8 cores) |
|-------|------------------|
| QC & Trimming | 30 min |
| Alignment | 2-4 hours |
| SV Calling | 1-3 hours |
| CNV Calling | 30 min - 2 hours |
| Integration | 15 min |
| Annotation | 30 min |

### Getting Help

For technical support:

1. Check the [main troubleshooting guide](troubleshooting.md)
2. Review tool-specific documentation
3. Search existing GitHub issues
4. Create a new issue with detailed error logs

For biological interpretation:

1. Consult the AnnotSV manual for pathogenicity guidelines
2. Review population databases (gnomAD-SV, DGV)
3. Consider gene-disease associations (OMIM, ClinVar)
4. Seek clinical genetics consultation for medical cases