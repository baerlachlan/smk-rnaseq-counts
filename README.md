# Snakemake workflow for generating gene-level counts from RNA-seq data

This Snakemake workflow implements the pre-processing steps to achieve gene-level read counts from raw RNA-seq data.

## Features

- Downloading/indexing of genome and annotations
- Trimming (fastp)
- Merging of sequencing units
- Alignment (STAR 2-pass mode)
- Counting (Subread featureCounts)
- Quality reports (FastQC)

## Usage

Standardised usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=baerlachlan/smk-rnaseq-de-counts).

However, if the intention is to run this workflow in an HPC environment where internet access is not be available on compute nodes, downloading the workflow from [Releases](https://github.com/baerlachlan/smk-rnaseq-de-counts/releases) is recommended.
