#!/bin/bash
set -e

# Run GATK Pipeline using Snakemake
# Usage: ./run_gatk.sh [cores]

CORES=${1:-1}

echo "Running GATK pipeline using Snakemake with $CORES cores"

# Change to GATK workflow directory
cd workflows/gatk

# Run Snakemake
snakemake -j "$CORES" --use-conda

echo "GATK pipeline completed successfully!"

# Return to original directory
cd ../..