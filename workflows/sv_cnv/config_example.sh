#!/bin/bash

# Configuration file for SV/CNV analysis pipeline
# Copy this file and modify the paths according to your setup

# === REQUIRED PARAMETERS ===

# Sample identifier
SAMPLE="NA12878"

# Number of threads to use
THREADS=8

# Path to reference genome FASTA file (must be indexed with BWA and samtools)
REF_GENOME="/path/to/reference/hg38.fasta"

# Path to adapter sequences for trimming
ADAPTERS="/path/to/adapters/TruSeq3-PE.fa"

# Directory containing raw FASTQ files
# Files should be named: ${SAMPLE}_R1.fastq.gz and ${SAMPLE}_R2.fastq.gz
RAW_READS_DIR="/path/to/raw_reads"

# Output directory for all analysis results
OUTPUT_DIR="/path/to/output"

# === OPTIONAL PARAMETERS ===

# Tool paths (leave as command name if installed via conda/in PATH)
TRIMMOMATIC_JAR="trimmomatic"  # or path to jar file: "/path/to/trimmomatic.jar"
PICARD_JAR="picard"            # or path to jar file: "/path/to/picard.jar"
CNVNATOR_PATH="cnvnator"       # or path to executable: "/path/to/cnvnator/src"
ANNOTSV_ROOT="AnnotSV"         # or path to installation: "/path/to/AnnotSV"

# CNVnator parameters
BIN_SIZE=1000                  # Bin size for CNV analysis (1000 for 30x coverage)
REF_CHROM_DIR="/path/to/reference/hg38_chrom_dirs"  # Directory with chromosome FASTA files for CNVnator

# GATK gCNV (optional)
GATK_GCNV_ENABLED="false"      # Set to "true" to enable GATK gCNV calling

# === EXAMPLE CONFIGURATION ===
# Uncomment and modify the following section for a typical setup:

# SAMPLE="sample_001"
# THREADS=16
# REF_GENOME="/data/reference/GRCh38.fasta"
# ADAPTERS="/data/adapters/TruSeq3-PE.fa"
# RAW_READS_DIR="/data/fastq"
# OUTPUT_DIR="/data/analysis/sample_001"
# REF_CHROM_DIR="/data/reference/GRCh38_chromosomes"