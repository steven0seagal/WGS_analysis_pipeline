#!/bin/bash

# Comprehensive pipeline for SV and CNV analysis from WGS data.
# Based on the research outlined in research2 document
# This script implements a 6-stage analysis framework:
# 1. Data Ingestion & Quality Assurance
# 2. Alignment & Preparation
# 3. Structural Variant (SV) Calling
# 4. Copy Number Variation (CNV) Calling
# 5. Integration & Refinement
# 6. Functional Annotation

# --- Configuration ---
set -e # Exit immediately if a command exits with a non-zero status
set -o pipefail # Return value of a pipeline is the status of the last command to exit with a non-zero status

# Check if config file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <config_file>"
    echo "Example: $0 config.sh"
    exit 1
fi

# Source the configuration file
source "$1"

# Validate required variables
required_vars=(
    "SAMPLE"
    "THREADS"
    "REF_GENOME"
    "ADAPTERS"
    "RAW_READS_DIR"
    "OUTPUT_DIR"
)

for var in "${required_vars[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Required variable $var is not set in config file"
        exit 1
    fi
done

# Set default values for optional variables
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-trimmomatic}"
PICARD_JAR="${PICARD_JAR:-picard}"
CNVNATOR_PATH="${CNVNATOR_PATH:-cnvnator}"
ANNOTSV_ROOT="${ANNOTSV_ROOT:-AnnotSV}"
BIN_SIZE="${BIN_SIZE:-1000}"

# Create output directory structure
mkdir -p "${OUTPUT_DIR}"/{qc/{raw,trimmed,multiqc_report},trimmed_reads,bam,analysis_ready_bam,sv_calls/{manta,delly,lumpy},cnv_calls/{cnvnator,gatk_gcnv},merged_filtered_sv,annotations}

# Change to output directory
cd "${OUTPUT_DIR}"

# --- Pipeline Functions ---

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: Command '$1' not found. Please ensure it is installed and in PATH."
        exit 1
    fi
}

# Stage 1: QC and Trimming
run_qc_and_trim() {
    log "--- STAGE 1: QC and Trimming ---"

    # Check required commands
    check_command fastqc
    check_command trimmomatic
    check_command multiqc

    # Initial QC
    log "Running FastQC on raw reads..."
    fastqc -o qc/raw --threads ${THREADS} ${RAW_READS_DIR}/${SAMPLE}_R1.fastq.gz
    fastqc -o qc/raw --threads ${THREADS} ${RAW_READS_DIR}/${SAMPLE}_R2.fastq.gz

    # Trimming with Trimmomatic
    log "Running Trimmomatic for read trimming..."
    if [[ "$TRIMMOMATIC_JAR" == "trimmomatic" ]]; then
        # Use trimmomatic command directly (conda version)
        trimmomatic PE -threads ${THREADS} \
            ${RAW_READS_DIR}/${SAMPLE}_R1.fastq.gz ${RAW_READS_DIR}/${SAMPLE}_R2.fastq.gz \
            trimmed_reads/${SAMPLE}_R1.trimmed.fastq.gz trimmed_reads/${SAMPLE}_R1.unpaired.fastq.gz \
            trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz trimmed_reads/${SAMPLE}_R2.unpaired.fastq.gz \
            ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
    else
        # Use java jar version
        java -jar ${TRIMMOMATIC_JAR} PE -threads ${THREADS} \
            ${RAW_READS_DIR}/${SAMPLE}_R1.fastq.gz ${RAW_READS_DIR}/${SAMPLE}_R2.fastq.gz \
            trimmed_reads/${SAMPLE}_R1.trimmed.fastq.gz trimmed_reads/${SAMPLE}_R1.unpaired.fastq.gz \
            trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz trimmed_reads/${SAMPLE}_R2.unpaired.fastq.gz \
            ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
    fi

    # Post-trimming QC
    log "Running FastQC on trimmed reads..."
    fastqc -o qc/trimmed --threads ${THREADS} trimmed_reads/${SAMPLE}_R1.trimmed.fastq.gz
    fastqc -o qc/trimmed --threads ${THREADS} trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz

    # Aggregate QC reports
    log "Generating MultiQC report..."
    multiqc qc -o qc/multiqc_report

    log "--- STAGE 1 Complete ---"
}

# Stage 2: Alignment and BAM Preparation
run_alignment() {
    log "--- STAGE 2: Alignment and BAM Preparation ---"

    # Check required commands
    check_command bwa
    check_command samtools

    # BWA-MEM alignment
    log "Running BWA-MEM alignment..."
    bwa mem -t ${THREADS} -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        ${REF_GENOME} \
        trimmed_reads/${SAMPLE}_R1.trimmed.fastq.gz \
        trimmed_reads/${SAMPLE}_R2.trimmed.fastq.gz | \
        samtools view -Sb - > bam/${SAMPLE}.bam

    # Sort and index
    log "Sorting and indexing BAM file..."
    samtools sort -@ ${THREADS} -o bam/${SAMPLE}.sorted.bam bam/${SAMPLE}.bam
    samtools index -@ ${THREADS} bam/${SAMPLE}.sorted.bam

    # Mark duplicates
    log "Marking duplicates with Picard..."
    if [[ "$PICARD_JAR" == "picard" ]]; then
        # Use picard command directly (conda version)
        picard MarkDuplicates \
            INPUT=bam/${SAMPLE}.sorted.bam \
            OUTPUT=analysis_ready_bam/${SAMPLE}.dedup.bam \
            METRICS_FILE=qc/${SAMPLE}.dedup_metrics.txt \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=STRICT
    else
        # Use java jar version
        java -jar ${PICARD_JAR} MarkDuplicates \
            INPUT=bam/${SAMPLE}.sorted.bam \
            OUTPUT=analysis_ready_bam/${SAMPLE}.dedup.bam \
            METRICS_FILE=qc/${SAMPLE}.dedup_metrics.txt \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=STRICT
    fi

    # Clean up intermediate files
    rm bam/${SAMPLE}.bam bam/${SAMPLE}.sorted.bam bam/${SAMPLE}.sorted.bam.bai

    log "--- STAGE 2 Complete ---"
}

# Stage 3: SV Calling
run_sv_calling() {
    log "--- STAGE 3: Ensemble SV Calling ---"

    # Manta
    log "Running Manta..."
    check_command configManta.py
    configManta.py --bam analysis_ready_bam/${SAMPLE}.dedup.bam \
                   --referenceFasta ${REF_GENOME} \
                   --runDir sv_calls/manta/${SAMPLE}
    python sv_calls/manta/${SAMPLE}/runWorkflow.py -m local -j ${THREADS}

    # Delly
    log "Running Delly..."
    check_command delly
    check_command bcftools
    delly call -t DEL -g ${REF_GENOME} -o sv_calls/delly/${SAMPLE}.del.bcf analysis_ready_bam/${SAMPLE}.dedup.bam
    delly call -t DUP -g ${REF_GENOME} -o sv_calls/delly/${SAMPLE}.dup.bcf analysis_ready_bam/${SAMPLE}.dedup.bam
    delly call -t INV -g ${REF_GENOME} -o sv_calls/delly/${SAMPLE}.inv.bcf analysis_ready_bam/${SAMPLE}.dedup.bam
    delly call -t TRA -g ${REF_GENOME} -o sv_calls/delly/${SAMPLE}.tra.bcf analysis_ready_bam/${SAMPLE}.dedup.bam

    # Merge Delly results
    bcftools concat -a sv_calls/delly/${SAMPLE}.*.bcf | \
        bcftools sort -Oz -o sv_calls/delly/${SAMPLE}.merged.vcf.gz
    tabix -p vcf sv_calls/delly/${SAMPLE}.merged.vcf.gz

    # Lumpy
    log "Running Lumpy..."
    check_command lumpyexpress
    check_command extractSplitReads_BwaMem

    # Extract discordant and split reads for Lumpy
    samtools view -b -F 1294 analysis_ready_bam/${SAMPLE}.dedup.bam > sv_calls/lumpy/${SAMPLE}.discordants.bam
    samtools view -h analysis_ready_bam/${SAMPLE}.dedup.bam | \
        extractSplitReads_BwaMem -i stdin | \
        samtools view -Sb - > sv_calls/lumpy/${SAMPLE}.splitters.bam

    # Sort the extracted reads
    samtools sort -@ ${THREADS} -o sv_calls/lumpy/${SAMPLE}.discordants.sorted.bam sv_calls/lumpy/${SAMPLE}.discordants.bam
    samtools sort -@ ${THREADS} -o sv_calls/lumpy/${SAMPLE}.splitters.sorted.bam sv_calls/lumpy/${SAMPLE}.splitters.bam

    # Run Lumpy
    lumpyexpress -B analysis_ready_bam/${SAMPLE}.dedup.bam \
                 -S sv_calls/lumpy/${SAMPLE}.splitters.sorted.bam \
                 -D sv_calls/lumpy/${SAMPLE}.discordants.sorted.bam \
                 -o sv_calls/lumpy/${SAMPLE}.vcf

    # Clean up intermediate files
    rm sv_calls/delly/${SAMPLE}.*.bcf
    rm sv_calls/lumpy/${SAMPLE}.discordants.bam sv_calls/lumpy/${SAMPLE}.splitters.bam

    log "--- STAGE 3 Complete ---"
}

# Stage 4: CNV Calling
run_cnv_calling() {
    log "--- STAGE 4: Read-Depth CNV Calling ---"

    # CNVnator
    log "Running CNVnator..."
    if [[ "$CNVNATOR_PATH" == "cnvnator" ]]; then
        # Use cnvnator command directly
        cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                 -tree analysis_ready_bam/${SAMPLE}.dedup.bam \
                 -chrom $(seq 1 22) X Y
        cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                 -his ${BIN_SIZE} \
                 -d ${REF_CHROM_DIR}
        cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                 -stat ${BIN_SIZE}
        cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                 -partition ${BIN_SIZE}
        cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                 -call ${BIN_SIZE} > cnv_calls/cnvnator/${SAMPLE}.cnvnator.out
    else
        # Use full path version
        ${CNVNATOR_PATH}/cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                                  -tree analysis_ready_bam/${SAMPLE}.dedup.bam \
                                  -chrom $(seq 1 22) X Y
        ${CNVNATOR_PATH}/cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                                  -his ${BIN_SIZE} \
                                  -d ${REF_CHROM_DIR}
        ${CNVNATOR_PATH}/cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                                  -stat ${BIN_SIZE}
        ${CNVNATOR_PATH}/cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                                  -partition ${BIN_SIZE}
        ${CNVNATOR_PATH}/cnvnator -root cnv_calls/cnvnator/${SAMPLE}.root \
                                  -call ${BIN_SIZE} > cnv_calls/cnvnator/${SAMPLE}.cnvnator.out
    fi

    # Convert CNVnator output to VCF
    if command -v cnvnator2VCF.pl &> /dev/null; then
        cnvnator2VCF.pl cnv_calls/cnvnator/${SAMPLE}.cnvnator.out > cnv_calls/cnvnator/${SAMPLE}.vcf
    else
        log "Warning: cnvnator2VCF.pl not found. CNVnator output kept in native format."
    fi

    # Optional: GATK gCNV (if GATK is available and configured)
    if command -v gatk &> /dev/null && [ -n "${GATK_GCNV_ENABLED}" ] && [ "${GATK_GCNV_ENABLED}" = "true" ]; then
        log "Running GATK gCNV..."

        # Preprocess intervals
        gatk PreprocessIntervals \
            -R ${REF_GENOME} \
            --bin-length ${BIN_SIZE} \
            --padding 0 \
            -O cnv_calls/gatk_gcnv/preprocessed.interval_list

        # Collect read counts
        gatk CollectReadCounts \
            -I analysis_ready_bam/${SAMPLE}.dedup.bam \
            -L cnv_calls/gatk_gcnv/preprocessed.interval_list \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O cnv_calls/gatk_gcnv/${SAMPLE}.counts.hdf5

        # Determine contig ploidy (simplified for single sample)
        gatk DetermineGermlineContigPloidy \
            -I cnv_calls/gatk_gcnv/${SAMPLE}.counts.hdf5 \
            --output cnv_calls/gatk_gcnv/${SAMPLE}-ploidy-calls \
            --output-prefix ${SAMPLE}

        # Call CNVs
        gatk GermlineCNVCaller \
            --run-mode CASE \
            -I cnv_calls/gatk_gcnv/${SAMPLE}.counts.hdf5 \
            --contig-ploidy-calls cnv_calls/gatk_gcnv/${SAMPLE}-ploidy-calls \
            --output cnv_calls/gatk_gcnv/${SAMPLE}-gcnv-calls \
            --output-prefix ${SAMPLE}

        # Post-process calls
        gatk PostprocessGermlineCNVCalls \
            --model-shard-path cnv_calls/gatk_gcnv/${SAMPLE}-gcnv-calls/shard-0 \
            --calls-shard-path cnv_calls/gatk_gcnv/${SAMPLE}-gcnv-calls/shard-0 \
            --allosomal-contig X --allosomal-contig Y \
            --contig-ploidy-calls cnv_calls/gatk_gcnv/${SAMPLE}-ploidy-calls \
            --sample-index 0 \
            --output-genotyped-intervals cnv_calls/gatk_gcnv/${SAMPLE}.intervals.vcf \
            --output-genotyped-segments cnv_calls/gatk_gcnv/${SAMPLE}.segments.vcf
    fi

    log "--- STAGE 4 Complete ---"
}

# Stage 5: Integration and Filtering
run_integration_and_filtering() {
    log "--- STAGE 5: Integration and Filtering ---"

    check_command SURVIVOR

    # Merge SV calls with SURVIVOR
    log "Merging SV calls with SURVIVOR..."
    ls sv_calls/manta/${SAMPLE}/results/variants/diploidSV.vcf.gz \
       sv_calls/delly/${SAMPLE}.merged.vcf.gz \
       sv_calls/lumpy/${SAMPLE}.vcf > merged_filtered_sv/sv_vcf_list.txt

    SURVIVOR merge merged_filtered_sv/sv_vcf_list.txt 1000 2 1 1 0 50 merged_filtered_sv/${SAMPLE}.sv.merged.vcf

    # Filter merged SV calls
    log "Filtering SV calls..."
    bcftools filter -i 'SVTYPE!="BND" && QUAL>30' \
        merged_filtered_sv/${SAMPLE}.sv.merged.vcf | \
        bgzip -c > merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz
    tabix -p vcf merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz

    # Filter CNV calls (if VCF format is available)
    if [ -f "cnv_calls/cnvnator/${SAMPLE}.vcf" ]; then
        log "Filtering CNV calls..."
        bcftools filter -i 'INFO/SVLEN>5000 || INFO/SVLEN<-5000' \
            cnv_calls/cnvnator/${SAMPLE}.vcf | \
            bgzip -c > merged_filtered_sv/${SAMPLE}.cnv.filtered.vcf.gz
        tabix -p vcf merged_filtered_sv/${SAMPLE}.cnv.filtered.vcf.gz
    fi

    log "--- STAGE 5 Complete ---"
}

# Stage 6: Annotation
run_annotation() {
    log "--- STAGE 6: Functional Annotation ---"

    # AnnotSV annotation
    if command -v AnnotSV &> /dev/null; then
        log "Running AnnotSV..."
        if [[ "$ANNOTSV_ROOT" == "AnnotSV" ]]; then
            AnnotSV -SVinputFile merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz \
                    -outputFile annotations/${SAMPLE}.sv.annotsv.tsv \
                    -genomeBuild GRCh38
        else
            ${ANNOTSV_ROOT}/bin/AnnotSV \
                -SVinputFile merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz \
                -outputFile annotations/${SAMPLE}.sv.annotsv.tsv \
                -genomeBuild GRCh38
        fi
    else
        log "Warning: AnnotSV not found. Skipping AnnotSV annotation."
    fi

    # VEP annotation
    if command -v vep &> /dev/null; then
        log "Running VEP..."
        vep --input_file merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz \
            --output_file annotations/${SAMPLE}.sv.vep.vcf \
            --vcf --cache --offline --fork ${THREADS} --fasta ${REF_GENOME} \
            --symbol --terms SO --hgvs
    else
        log "Warning: VEP not found. Skipping VEP annotation."
    fi

    log "--- STAGE 6 Complete ---"
}

# Function to generate summary report
generate_summary() {
    log "--- Generating Summary Report ---"

    SUMMARY_FILE="summary_${SAMPLE}.txt"

    cat > ${SUMMARY_FILE} << EOF
SV/CNV Analysis Summary for Sample: ${SAMPLE}
Date: $(date)
========================================

Input Files:
- R1: ${RAW_READS_DIR}/${SAMPLE}_R1.fastq.gz
- R2: ${RAW_READS_DIR}/${SAMPLE}_R2.fastq.gz
- Reference: ${REF_GENOME}

Output Files:
- Analysis-ready BAM: analysis_ready_bam/${SAMPLE}.dedup.bam
- Merged SV calls: merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz
- CNV calls: cnv_calls/cnvnator/${SAMPLE}.cnvnator.out
- MultiQC report: qc/multiqc_report/multiqc_report.html

Variant Calling Results:
EOF

    # Add counts if files exist
    if [ -f "merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz" ]; then
        SV_COUNT=$(bcftools view -H merged_filtered_sv/${SAMPLE}.sv.merged.filtered.vcf.gz | wc -l)
        echo "- Structural variants: ${SV_COUNT}" >> ${SUMMARY_FILE}
    fi

    if [ -f "cnv_calls/cnvnator/${SAMPLE}.cnvnator.out" ]; then
        CNV_COUNT=$(grep -v "^#" cnv_calls/cnvnator/${SAMPLE}.cnvnator.out | wc -l)
        echo "- Copy number variants: ${CNV_COUNT}" >> ${SUMMARY_FILE}
    fi

    echo "" >> ${SUMMARY_FILE}
    echo "Analysis completed successfully at $(date)" >> ${SUMMARY_FILE}

    log "Summary report saved to: ${SUMMARY_FILE}"
}

# --- Main Execution ---
main() {
    log "Starting SV/CNV analysis pipeline for sample: ${SAMPLE}"
    log "Output directory: ${OUTPUT_DIR}"

    # Check if input files exist
    if [ ! -f "${RAW_READS_DIR}/${SAMPLE}_R1.fastq.gz" ] || [ ! -f "${RAW_READS_DIR}/${SAMPLE}_R2.fastq.gz" ]; then
        echo "Error: Input FASTQ files not found in ${RAW_READS_DIR}"
        exit 1
    fi

    if [ ! -f "${REF_GENOME}" ]; then
        echo "Error: Reference genome not found: ${REF_GENOME}"
        exit 1
    fi

    # Run pipeline stages
    run_qc_and_trim
    run_alignment
    run_sv_calling
    run_cnv_calling
    run_integration_and_filtering
    run_annotation
    generate_summary

    log "--- PIPELINE FINISHED SUCCESSFULLY ---"
}

# Execute main function
main "$@"