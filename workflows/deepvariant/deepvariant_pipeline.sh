#!/bin/bash
set -e

# DeepVariant Pipeline using Docker
# Usage: ./deepvariant_pipeline.sh <sample_id> <input_dir> <output_dir> <ref_genome> <bam_file>

if [ $# -ne 5 ]; then
    echo "Usage: $0 <sample_id> <input_dir> <output_dir> <ref_genome> <bam_file>"
    exit 1
fi

SAMPLE_ID=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
REF_GENOME_HOST=$4
BAM_HOST=$5

BIN_VERSION="1.6.1"

# Define filenames for use INSIDE the container
REF_GENOME_CONTAINER="/input/$(basename ${REF_GENOME_HOST})"
BAM_CONTAINER="/input/$(basename ${BAM_HOST})"
OUTPUT_VCF_CONTAINER="/output/${SAMPLE_ID}.vcf.gz"
OUTPUT_GVCF_CONTAINER="/output/${SAMPLE_ID}.g.vcf.gz"
LOG_DIR_CONTAINER="/output/logs"

# Create the output directory on the host
mkdir -p ${OUTPUT_DIR}

echo "Running DeepVariant for ${SAMPLE_ID}..."
docker run --gpus all \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}":"/output" \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=${REF_GENOME_CONTAINER} \
  --reads=${BAM_CONTAINER} \
  --output_vcf=${OUTPUT_VCF_CONTAINER} \
  --output_gvcf=${OUTPUT_GVCF_CONTAINER} \
  --num_shards=$(nproc) \
  --logging_dir=${LOG_DIR_CONTAINER}

echo "DeepVariant pipeline completed successfully for ${SAMPLE_ID}."