#!/bin/bash
set -e

# Run DeepVariant Pipeline
# Usage: ./run_deepvariant.sh <sample_id> <input_dir> <output_dir> <ref_genome> <bam_file>

if [ $# -ne 5 ]; then
    echo "Usage: $0 <sample_id> <input_dir> <output_dir> <ref_genome> <bam_file>"
    echo "Example: $0 Sample1 /path/to/inputs /path/to/outputs /path/to/ref.fasta /path/to/sample.bam"
    exit 1
fi

SAMPLE_ID=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
REF_GENOME=$4
BAM_FILE=$5

echo "Running DeepVariant pipeline for sample: $SAMPLE_ID"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Reference genome: $REF_GENOME"
echo "BAM file: $BAM_FILE"

# Run the DeepVariant pipeline script
./workflows/deepvariant/deepvariant_pipeline.sh "$SAMPLE_ID" "$INPUT_DIR" "$OUTPUT_DIR" "$REF_GENOME" "$BAM_FILE"

echo "DeepVariant pipeline completed successfully!"