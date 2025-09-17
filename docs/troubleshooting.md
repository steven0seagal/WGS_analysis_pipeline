# Troubleshooting Guide

This guide helps resolve common issues when running the WGS analysis pipelines.

## General Issues

### Conda Environment Issues

**Problem:** `conda env create -f environment.yml` fails

**Solutions:**
- Update conda: `conda update conda`
- Clear conda cache: `conda clean --all`
- Use mamba instead: `mamba env create -f environment.yml`

### Memory Issues

**Problem:** Pipeline runs out of memory

**Solutions:**
- Reduce threads in Snakemake: `snakemake -j 1`
- Increase Java heap size in config.yaml: `java_opts: "-Xmx16g"`
- Use a machine with more RAM

## GATK Pipeline Issues

### GenomicsDBImport Fails

**Problem:** `GenomicsDBImport` throws errors about intervals

**Solutions:**
- Ensure intervals in config.yaml are valid (e.g., chr1, chr2, not 1, 2)
- Check reference genome chromosome names match BAM headers
- Reduce number of samples or use smaller intervals

### VQSR Fails

**Problem:** VariantRecalibrator has insufficient data

**Solutions:**
- Ensure known variant resources are properly indexed
- Check that VCF files are bgzipped and tabix indexed
- For small datasets, consider skipping VQSR or using hard filtering

## DeepVariant Pipeline Issues

### Docker Issues

**Problem:** `docker run` fails with permission denied

**Solutions:**
- Add user to docker group: `sudo usermod -aG docker $USER`
- Run with sudo: `sudo docker run ...`
- For GPU support, ensure NVIDIA drivers are installed

**Problem:** GPU not detected

**Solutions:**
- Install NVIDIA Container Toolkit
- Check GPU availability: `nvidia-smi`
- Remove `--gpus all` flag for CPU-only execution

### Model Type Errors

**Problem:** Invalid model type specified

**Solutions:**
- Use valid model types: WGS, WES, PACBIO
- Ensure model matches your data type

## Preprocessing Issues

### BWA Alignment Fails

**Problem:** `bwa mem` produces empty SAM

**Solutions:**
- Check FASTQ quality and format
- Verify reference genome is properly indexed
- Ensure read group information is correct

### MarkDuplicates Fails

**Problem:** Picard MarkDuplicates throws errors

**Solutions:**
- Ensure BAM is coordinate-sorted
- Check for corrupted input files
- Verify sufficient disk space

## Annotation Issues

### SnpEff Database Missing

**Problem:** SnpEff cannot find database

**Solutions:**
- Download database: `snpEff download GRCh38.105`
- Check genome version in config matches downloaded database
- Update paths in config.yaml

### VEP Cache Issues

**Problem:** VEP cannot find cache

**Solutions:**
- Download cache: `vep_install --AUTO c --SPECIES homo_sapiens --ASSEMBLY GRCh38`
- Set correct cache directory in config.yaml
- Use offline mode for faster processing

## Performance Optimization

### Slow Pipeline Execution

**Solutions:**
- Increase cores: `snakemake -j 8`
- Use SSD storage for I/O intensive steps
- Pre-download all required resources
- Use local VEP cache instead of remote

### Large Output Files

**Solutions:**
- Compress intermediate files where possible
- Clean up temporary files regularly
- Use appropriate storage with sufficient space

## Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](https://github.com/steven0seagal/WGS_analysis/issues) for similar problems
2. Create a new issue with:
   - Error messages and logs
   - Your environment details
   - Steps to reproduce
   - Sample data if possible

3. Include relevant tool versions and system information