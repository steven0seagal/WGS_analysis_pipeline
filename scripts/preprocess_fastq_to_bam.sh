#!/bin/bash
set -e

# Script to preprocess FASTQ files to analysis-ready BAM
# Usage: ./preprocess_fastq_to_bam.sh <sample_id> <fastq_r1> <fastq_r2> <ref_genome> <known_sites_vcf>

if [ $# -ne 5 ]; then
    echo "Usage: $0 <sample_id> <fastq_r1> <fastq_r2> <ref_genome> <known_sites_vcf>"
    exit 1
fi

SAMPLE_ID=$1
FASTQ_R1=$2
FASTQ_R2=$3
REF_GENOME=$4
KNOWN_SITES_VCF=$5

OUTPUT_DIR="results/bams"
mkdir -p ${OUTPUT_DIR}

# Step 1: Read Alignment
OUTPUT_SAM="${OUTPUT_DIR}/${SAMPLE_ID}.sam"
echo "Aligning reads for ${SAMPLE_ID}..."
bwa mem -t 8 \
    -R "@RG\tID:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1\tSM:${SAMPLE_ID}" \
    ${REF_GENOME} \
    ${FASTQ_R1} \
    ${FASTQ_R2} > ${OUTPUT_SAM}

# Step 2: Coordinate Sorting and Format Conversion
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
echo "Sorting and converting to BAM..."
samtools sort -@ 8 -o ${SORTED_BAM} ${OUTPUT_SAM}
rm ${OUTPUT_SAM}  # Remove SAM file

# Step 3: Duplicate Marking
MARKED_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.marked_duplicates.bam"
METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_ID}.metrics.txt"
echo "Marking duplicates..."
gatk MarkDuplicates \
    -I ${SORTED_BAM} \
    -O ${MARKED_BAM} \
    -M ${METRICS_FILE} \
    --CREATE_INDEX true

# Step 4: Base Quality Score Recalibration
RECAL_TABLE="${OUTPUT_DIR}/${SAMPLE_ID}.recal_data.table"
ANALYSIS_READY_BAM="${OUTPUT_DIR}/${SAMPLE_ID}.analysis_ready.bam"
echo "Performing BQSR..."
gatk BaseRecalibrator \
    -R ${REF_GENOME} \
    -I ${MARKED_BAM} \
    --known-sites ${KNOWN_SITES_VCF} \
    -O ${RECAL_TABLE}

gatk ApplyBQSR \
    -R ${REF_GENOME} \
    -I ${MARKED_BAM} \
    --bqsr-recal-file ${RECAL_TABLE} \
    -O ${ANALYSIS_READY_BAM}

echo "Preprocessing completed for ${SAMPLE_ID}. Analysis-ready BAM: ${ANALYSIS_READY_BAM}"