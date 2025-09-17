#!/bin/bash
set -e

# GATK Germline Variant Discovery Pipeline
# Usage: ./gatk_pipeline.sh <config_file>

if [ $# -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    echo "Config file should define: GATK_PATH, REF_GENOME, DBSNP_VCF, HAPMAP_VCF, OMNI_VCF, PHASE1_VCF, MILLS_INDELS_VCF, BAM_DIR, GVCF_DIR, DB_DIR, VCF_DIR, RECAL_DIR, LOG_DIR, SAMPLES array"
    exit 1
fi

CONFIG_FILE=$1
source ${CONFIG_FILE}

# Create output directories
mkdir -p ${GVCF_DIR} ${DB_DIR} ${VCF_DIR} ${RECAL_DIR} ${LOG_DIR}

# --- 2. Per-Sample Variant Calling with HaplotypeCaller in GVCF mode ---
echo "Step 2: Running HaplotypeCaller for each sample..."
for SAMPLE_ID in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE_ID}..."
    gatk --java-options "-Xmx8g" HaplotypeCaller \
        -R ${REF_GENOME} \
        -I ${BAM_DIR}/${SAMPLE_ID}.analysis_ready.bam \
        -O ${GVCF_DIR}/${SAMPLE_ID}.g.vcf.gz \
        -ERC GVCF \
        -D ${DBSNP_VCF} > ${LOG_DIR}/${SAMPLE_ID}_haplotypecaller.log 2>&1
done
echo "HaplotypeCaller finished for all samples."

# --- 3. Consolidate GVCFs with GenomicsDBImport ---
echo "Step 3: Consolidating GVCFs with GenomicsDBImport..."
GVCF_ARGS=""
for SAMPLE_ID in "${SAMPLES[@]}"; do
    GVCF_ARGS+=" -V ${GVCF_DIR}/${SAMPLE_ID}.g.vcf.gz"
done

gatk --java-options "-Xmx8g -Xms8g" GenomicsDBImport \
    ${GVCF_ARGS} \
    --genomicsdb-workspace-path ${DB_DIR}/cohort_db_chr20 \
    -L chr20 \
    --reader-threads 4 > ${LOG_DIR}/genomicsdbimport.log 2>&1
echo "GenomicsDBImport finished."

# --- 4. Joint Genotyping with GenotypeGVCFs ---
echo "Step 4: Running GenotypeGVCFs for joint calling..."
gatk --java-options "-Xmx8g" GenotypeGVCFs \
    -R ${REF_GENOME} \
    -V gendb://${DB_DIR}/cohort_db_chr20 \
    -O ${VCF_DIR}/raw_variants_chr20.vcf.gz \
    -D ${DBSNP_VCF} > ${LOG_DIR}/genotypegvcfs.log 2>&1
echo "GenotypeGVCFs finished."

# --- 5. Variant Quality Score Recalibration (VQSR) ---
echo "Step 5a: VQSR for SNPs..."
gatk VariantRecalibrator \
    -R ${REF_GENOME} \
    -V ${VCF_DIR}/raw_variants_chr20.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP_VCF} \
    --resource:omni,known=false,training=true,truth=true,prior=12.0 ${OMNI_VCF} \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${PHASE1_VCF} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP_VCF} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -O ${RECAL_DIR}/cohort_snps.recal \
    --tranches-file ${RECAL_DIR}/cohort_snps.tranches \
    --rscript-file ${RECAL_DIR}/cohort_snps.plots.R > ${LOG_DIR}/vqsr_snp_recal.log 2>&1

gatk ApplyVQSR \
    -R ${REF_GENOME} \
    -V ${VCF_DIR}/raw_variants_chr20.vcf.gz \
    -O ${VCF_DIR}/recalibrated_snps_chr20.vcf.gz \
    --recal-file ${RECAL_DIR}/cohort_snps.recal \
    --tranches-file ${RECAL_DIR}/cohort_snps.tranches \
    --truth-sensitivity-filter-level 99.5 \
    -mode SNP > ${LOG_DIR}/vqsr_snp_apply.log 2>&1
echo "SNP recalibration finished."

echo "Step 5b: VQSR for Indels..."
gatk VariantRecalibrator \
    -R ${REF_GENOME} \
    -V ${VCF_DIR}/recalibrated_snps_chr20.vcf.gz \
    --resource:mills,known=false,training=true,truth=true,prior=12.0 ${MILLS_INDELS_VCF} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP_VCF} \
    -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
    -mode INDEL \
    -O ${RECAL_DIR}/cohort_indels.recal \
    --tranches-file ${RECAL_DIR}/cohort_indels.tranches \
    --rscript-file ${RECAL_DIR}/cohort_indels.plots.R > ${LOG_DIR}/vqsr_indel_recal.log 2>&1

gatk ApplyVQSR \
    -R ${REF_GENOME} \
    -V ${VCF_DIR}/recalibrated_snps_chr20.vcf.gz \
    -O ${VCF_DIR}/final_filtered_variants_chr20.vcf.gz \
    --recal-file ${RECAL_DIR}/cohort_indels.recal \
    --tranches-file ${RECAL_DIR}/cohort_indels.tranches \
    --truth-sensitivity-filter-level 99.0 \
    -mode INDEL > ${LOG_DIR}/vqsr_indel_apply.log 2>&1
echo "Indel recalibration finished."

echo "GATK pipeline completed successfully!"