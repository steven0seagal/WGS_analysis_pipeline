# Pipeline Descriptions

## Preprocessing Pipeline

The preprocessing pipeline converts raw FASTQ files into analysis-ready BAM files suitable for variant calling.

### Steps:
1. **Read Alignment**: BWA-MEM aligns paired-end reads to the reference genome
2. **Sorting and BAM Conversion**: Samtools sorts alignments and converts to BAM format
3. **Duplicate Marking**: GATK MarkDuplicates identifies and flags PCR duplicates
4. **Base Quality Score Recalibration (BQSR)**: GATK corrects systematic base quality errors

## GATK Germline Pipeline

Implements the GATK Best Practices for germline variant discovery.

### Key Features:
- **GVCF Workflow**: Scalable joint calling across multiple samples
- **HaplotypeCaller**: Local de novo assembly for accurate variant detection
- **VQSR**: Machine learning-based filtering using known variant resources
- **Joint Genotyping**: Combines evidence across samples

### Output:
- Filtered VCF with high-confidence variants
- Annotated variants with functional consequences

## DeepVariant Pipeline

Uses deep learning for variant calling, reframing the task as image classification.

### Key Features:
- **CNN-based Calling**: Convolutional neural network trained on diverse genomes
- **Pileup Images**: Transforms read alignments into multi-channel tensors
- **Model Selection**: Optimized models for different sequencing technologies
- **Containerized Execution**: Ensures reproducibility across environments

### Output:
- VCF and GVCF files with variant calls
- HTML report with quality metrics

## Annotation

### SnpEff
- Fast annotation with pre-built databases
- Predicts functional effects (missense, nonsense, etc.)
- Gene and transcript context

### Ensembl VEP
- Comprehensive annotation from Ensembl
- Extensible with plugins (CADD, SpliceAI, gnomAD)
- Population frequencies and clinical significance