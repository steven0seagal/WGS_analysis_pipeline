#!/bin/bash
set -e

# Script to prepare reference genome for WGS analysis
# Usage: ./prepare_reference.sh <reference.fasta>

if [ $# -ne 1 ]; then
    echo "Usage: $0 <reference.fasta>"
    exit 1
fi

REF_GENOME=$1

echo "Preparing reference genome: $REF_GENOME"

# 1. Index for BWA
echo "Indexing for BWA..."
bwa index ${REF_GENOME}

# 2. Index for Samtools
echo "Indexing for Samtools..."
samtools faidx ${REF_GENOME}

# 3. Create Sequence Dictionary for GATK
echo "Creating sequence dictionary for GATK..."
gatk CreateSequenceDictionary -R ${REF_GENOME}

echo "Reference genome preparation completed."