#!/bin/bash

# Convenience script for running SV/CNV analysis pipeline
# Usage: ./run_sv_cnv.sh <sample_name> <fastq_dir> <output_dir> <reference_fasta> [threads] [mode]

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Function to display usage
usage() {
    echo "Usage: $0 <sample_name> <fastq_dir> <output_dir> <reference_fasta> [threads] [mode]"
    echo ""
    echo "Arguments:"
    echo "  sample_name     Sample identifier (e.g., NA12878)"
    echo "  fastq_dir       Directory containing FASTQ files"
    echo "  output_dir      Output directory for results"
    echo "  reference_fasta Path to reference genome FASTA file"
    echo "  threads         Number of threads to use (default: 8)"
    echo "  mode           Pipeline mode: 'bash' or 'snakemake' (default: snakemake)"
    echo ""
    echo "Examples:"
    echo "  $0 NA12878 /data/fastq /data/output /ref/hg38.fa"
    echo "  $0 NA12878 /data/fastq /data/output /ref/hg38.fa 16 bash"
    echo ""
    echo "Requirements:"
    echo "  - FASTQ files should be named: {sample_name}_R1.fastq.gz and {sample_name}_R2.fastq.gz"
    echo "  - Reference genome must be indexed (BWA and samtools)"
    echo "  - Conda/Mamba must be available for environment management"
    exit 1
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        echo "Error: File not found: $1"
        exit 1
    fi
}

# Function to check if directory exists
check_dir() {
    if [ ! -d "$1" ]; then
        echo "Error: Directory not found: $1"
        exit 1
    fi
}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Parse command line arguments
if [ $# -lt 4 ]; then
    usage
fi

SAMPLE_NAME="$1"
FASTQ_DIR="$2"
OUTPUT_DIR="$3"
REFERENCE_FASTA="$4"
THREADS="${5:-8}"
MODE="${6:-snakemake}"

# Validate inputs
log "Validating inputs..."

# Check required directories and files
check_dir "$FASTQ_DIR"
check_file "$REFERENCE_FASTA"

# Check FASTQ files
R1_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_R1.fastq.gz"
R2_FILE="${FASTQ_DIR}/${SAMPLE_NAME}_R2.fastq.gz"
check_file "$R1_FILE"
check_file "$R2_FILE"

# Check reference indexing
check_file "${REFERENCE_FASTA}.fai"
BWA_INDEX="${REFERENCE_FASTA%.fa*}"
if [ ! -f "${BWA_INDEX}.bwt" ] && [ ! -f "${REFERENCE_FASTA}.bwt" ]; then
    echo "Error: BWA index not found. Please run 'bwa index $REFERENCE_FASTA'"
    exit 1
fi

# Validate mode
if [ "$MODE" != "bash" ] && [ "$MODE" != "snakemake" ]; then
    echo "Error: Mode must be 'bash' or 'snakemake'"
    exit 1
fi

# Check if conda/mamba is available
if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "Error: Conda or Mamba is required but not found in PATH"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

log "Starting SV/CNV analysis for sample: $SAMPLE_NAME"
log "Input FASTQ directory: $FASTQ_DIR"
log "Output directory: $OUTPUT_DIR"
log "Reference genome: $REFERENCE_FASTA"
log "Threads: $THREADS"
log "Mode: $MODE"

# Set up paths
SV_CNV_DIR="${SCRIPT_DIR}/workflows/sv_cnv"
cd "$OUTPUT_DIR"

if [ "$MODE" == "bash" ]; then
    log "Running bash script pipeline..."

    # Create temporary config file
    CONFIG_FILE="${OUTPUT_DIR}/sv_cnv_config_${SAMPLE_NAME}.sh"

    cat > "$CONFIG_FILE" << EOF
#!/bin/bash

# Auto-generated configuration for SV/CNV pipeline
# Sample: $SAMPLE_NAME
# Generated: $(date)

# Required parameters
SAMPLE="$SAMPLE_NAME"
THREADS=$THREADS
REF_GENOME="$REFERENCE_FASTA"
ADAPTERS="\${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"
RAW_READS_DIR="$FASTQ_DIR"
OUTPUT_DIR="$OUTPUT_DIR"

# Optional parameters
BIN_SIZE=1000
CNVNATOR_PATH="cnvnator"
ANNOTSV_ROOT="AnnotSV"
REF_CHROM_DIR="\$(dirname $REFERENCE_FASTA)/chromosomes"

# Tool paths (using conda versions)
TRIMMOMATIC_JAR="trimmomatic"
PICARD_JAR="picard"

# Enable GATK gCNV if desired
GATK_GCNV_ENABLED="false"
EOF

    log "Configuration file created: $CONFIG_FILE"

    # Run the bash pipeline
    bash "${SV_CNV_DIR}/sv_cnv_pipeline.sh" "$CONFIG_FILE"

    # Clean up config file
    rm "$CONFIG_FILE"

elif [ "$MODE" == "snakemake" ]; then
    log "Running Snakemake pipeline..."

    # Create temporary config and samples files
    CONFIG_FILE="${OUTPUT_DIR}/config.yaml"
    SAMPLES_FILE="${OUTPUT_DIR}/samples.tsv"

    # Determine adapter path
    ADAPTER_PATH="/usr/share/trimmomatic/adapters/TruSeq3-PE.fa"
    if [ ! -f "$ADAPTER_PATH" ]; then
        ADAPTER_PATH="\${CONDA_PREFIX}/share/trimmomatic/adapters/TruSeq3-PE.fa"
    fi

    cat > "$CONFIG_FILE" << EOF
# Auto-generated configuration for SV/CNV pipeline
# Sample: $SAMPLE_NAME
# Generated: $(date)

# Input data
samples_tsv: "samples.tsv"
raw_reads_dir: "$FASTQ_DIR"

# Reference files
reference_fasta: "$REFERENCE_FASTA"
adapters_fasta: "$ADAPTER_PATH"
reference_chrom_dir: "$(dirname $REFERENCE_FASTA)/chromosomes"

# Tool parameters
cnvnator_bin_size: 1000
survivor_max_dist: 1000
survivor_min_callers: 2
survivor_min_size: 50

# Filtering parameters
sv_filter_expr: 'SVTYPE!="BND" && QUAL>30'
cnv_filter_expr: 'INFO/SVLEN>5000 || INFO/SVLEN<-5000'

# Genome build
genome_build: "GRCh38"
EOF

    cat > "$SAMPLES_FILE" << EOF
sample
$SAMPLE_NAME
EOF

    log "Configuration files created in: $OUTPUT_DIR"

    # Check if snakemake is available
    if ! command -v snakemake &> /dev/null; then
        log "Snakemake not found. Installing in temporary environment..."
        if command -v mamba &> /dev/null; then
            mamba create -n snakemake_temp -c conda-forge -c bioconda snakemake -y
            source activate snakemake_temp
        else
            conda create -n snakemake_temp -c conda-forge -c bioconda snakemake -y
            source activate snakemake_temp
        fi
    fi

    # Run Snakemake pipeline
    log "Executing Snakemake workflow..."
    snakemake \
        --snakefile "${SV_CNV_DIR}/Snakefile" \
        --configfile "$CONFIG_FILE" \
        --use-conda \
        --cores "$THREADS" \
        --printshellcmds \
        --reason

    # Clean up temporary files
    rm "$CONFIG_FILE" "$SAMPLES_FILE"

    # Deactivate temporary environment if created
    if [ -n "$CONDA_DEFAULT_ENV" ] && [ "$CONDA_DEFAULT_ENV" == "snakemake_temp" ]; then
        conda deactivate
        conda env remove -n snakemake_temp -y
    fi
fi

log "SV/CNV analysis completed successfully!"
log "Results are available in: $OUTPUT_DIR"

# Display summary
if [ -f "${OUTPUT_DIR}/summary_${SAMPLE_NAME}.txt" ]; then
    echo ""
    echo "=== ANALYSIS SUMMARY ==="
    cat "${OUTPUT_DIR}/summary_${SAMPLE_NAME}.txt"
fi

echo ""
echo "Key output files:"
echo "- Analysis-ready BAM: analysis_ready_bam/${SAMPLE_NAME}.dedup.bam"
echo "- Structural variants: merged_filtered_sv/${SAMPLE_NAME}.sv.merged.filtered.vcf.gz"
echo "- Copy number variants: merged_filtered_sv/${SAMPLE_NAME}.cnv.filtered.vcf.gz"
echo "- Quality control report: qc/multiqc_report/multiqc_report.html"

if [ -f "annotations/${SAMPLE_NAME}.sv.annotsv.tsv" ]; then
    echo "- Annotated SVs: annotations/${SAMPLE_NAME}.sv.annotsv.tsv"
fi

echo ""
echo "Next steps:"
echo "1. Review the MultiQC report for data quality"
echo "2. Examine variant calls in IGV for manual validation"
echo "3. Filter variants based on your analysis criteria"
echo "4. Interpret results using the annotation files"