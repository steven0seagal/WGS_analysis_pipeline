# Setup Guide

## Environment Setup

1. Install Miniconda or Mamba if not already installed.

2. Create the WGS analysis environment:

```bash
conda env create -f environment.yml
conda activate wgs_analysis
```

## Reference Genome Preparation

Download the reference genome (e.g., GRCh38 from Ensembl or UCSC) and prepare it:

```bash
./scripts/prepare_reference.sh Homo_sapiens_assembly38.fasta
```

This will create the necessary index files (.fai, .dict, .bwt, etc.).

## Resource Files for GATK

Download the following resource files from the GATK resource bundle:

- dbSNP: dbsnp_146.hg38.vcf.gz
- HapMap: hapmap_3.3.hg38.vcf.gz
- 1000G Omni: 1000G_omni2.5.hg38.vcf.gz
- 1000G Phase 1: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
- Mills indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

Update the paths in the config.yaml files.

## DeepVariant Setup

Ensure Docker is installed. For GPU support, install NVIDIA Container Toolkit.

The DeepVariant pipeline will automatically pull the required Docker image.

## VEP Setup

For VEP annotation, download the cache:

```bash
# After activating environment
vep_install --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38 --CACHE_VERSION 105
```

Update the cache path in config.yaml.