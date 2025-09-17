# Examples

This document provides practical examples for running the WGS analysis pipelines.

## Complete Workflow Example

### 1. Setup

```bash
# Clone repository
git clone https://github.com/steven0seagal/WGS_analysis.git
cd WGS_analysis

# Create environment
conda env create -f environment.yml
conda activate wgs_analysis

# Download reference genome (GRCh38)
wget ftp://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa reference.fasta

# Prepare reference
./scripts/prepare_reference.sh reference.fasta
```

### 2. Download Test Data

```bash
# Create directories
mkdir -p data/fastq data/bam results

# Download sample FASTQ (small test dataset)
# Note: Replace with your actual data
wget https://example.com/sample_R1.fastq.gz -O data/fastq/sample_R1.fastq.gz
wget https://example.com/sample_R2.fastq.gz -O data/fastq/sample_R2.fastq.gz
```

### 3. Preprocessing

```bash
# Run preprocessing
./scripts/preprocess_fastq_to_bam.sh \
    Sample1 \
    data/fastq/sample_R1.fastq.gz \
    data/fastq/sample_R2.fastq.gz \
    reference.fasta \
    resources/dbsnp_146.hg38.vcf.gz
```

### 4. Run GATK Pipeline

```bash
# Edit config.yaml with your paths
vim workflows/gatk/config.yaml

# Run pipeline
./run_gatk.sh 4
```

### 5. Run DeepVariant Pipeline

```bash
# Edit config.yaml
vim workflows/deepvariant/config.yaml

# Run pipeline
./run_deepvariant.sh \
    Sample1 \
    data \
    results/deepvariant \
    reference.fasta \
    data/bam/Sample1.analysis_ready.bam
```

## Configuration Examples

### GATK Config (config.yaml)

```yaml
samples:
  - Sample1
  - Sample2

intervals:
  - chr20
  - chr21

paths:
  ref_genome: "/full/path/to/reference.fasta"
  bam_dir: "/full/path/to/data/bam"
  resources:
    dbsnp: "/full/path/to/resources/dbsnp_146.hg38.vcf.gz"
    hapmap: "/full/path/to/resources/hapmap_3.3.hg38.vcf.gz"
    omni: "/full/path/to/resources/1000G_omni2.5.hg38.vcf.gz"
    phase1: "/full/path/to/resources/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    mills: "/full/path/to/resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

params:
  java_opts: "-Xmx8g"
  vqsr:
    snp_sensitivity: "99.5"
    indel_sensitivity: "99.0"
```

### DeepVariant Config (config.yaml)

```yaml
samples:
  - Sample1
  - Sample2

paths:
  ref_genome: "/full/path/to/reference.fasta"
  bam_dir: "/full/path/to/data/bam"

params:
  deepvariant:
    version: "1.6.1"
    model_type: "WGS"
```

## Output Files

After successful pipeline execution, you should have:

### GATK Pipeline Outputs
```
results/
├── vcfs/
│   └── final_merged.vcf.gz          # Final filtered variants
├── annotated/
│   ├── final.snpeff.vcf            # SnpEff annotated
│   └── final.vep.vcf               # VEP annotated
└── logs/                           # Log files
```

### DeepVariant Pipeline Outputs
```
results/deepvariant/
├── Sample1/
│   ├── Sample1.vcf.gz              # Variants
│   ├── Sample1.g.vcf.gz           # GVCF
│   └── Sample1.vcf.gz.html        # Report
```

## Customizing for Your Data

### Multiple Samples

For cohort analysis with GATK:

1. Add all sample IDs to `samples` list in config.yaml
2. Ensure all BAM files are in the bam_dir
3. Run the pipeline - it will automatically handle joint genotyping

### Different Genome Builds

To use GRCh37 instead of GRCh38:

1. Download appropriate reference genome
2. Update resource file paths to GRCh37 versions
3. Change genome versions in annotation configs (e.g., GRCh37.87 for SnpEff)

### Exome Data

For WES instead of WGS:

1. For DeepVariant: Change model_type to "WES"
2. For GATK: Add exome capture intervals to config
3. Adjust VQSR parameters if needed

## Monitoring Pipeline Progress

### Snakemake

```bash
# Check status
snakemake -n

# View DAG
snakemake --dag | dot -Tpng > dag.png

# Unlock directory if needed
snakemake --unlock
```

### Logs

Check log files in `results/logs/` for detailed error messages and progress information.