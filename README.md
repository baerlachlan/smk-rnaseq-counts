# Snakemake workflow for generating gene-level counts from RNA-seq data

This Snakemake workflow implements the pre-processing steps to achieve gene-level read counts from raw RNA-seq data.

## Features

- Quality reports (`FastQC`)
- Downloading/indexing of genome and annotations
- Trimming (`fastp`)
- Merging of sequencing units (e.g. samples split across multiple lanes)
- Single- or paired-end read compatibility
- Alignment (`STAR` 2-pass mode)
- Counting (`Subread` featureCounts)

## Usage

Standardised usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=baerlachlan/smk-rnaseq-de-counts).

However, if the intention is to run this workflow in an HPC environment where internet access is not be available on compute nodes, downloading the workflow from [Releases](https://github.com/baerlachlan/smk-rnaseq-de-counts/releases) is recommended.


## Testing

Test data and configurations are available for quick testing of this workflow.

Simply copy the desired config directory as follows:

```bash
## Remove example config
rm -r config/
## Paired-end
cp -r .test/config_pe/ config/
## Single-end
cp -r .test/config_se/ config/
```
