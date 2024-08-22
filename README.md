# Snakemake workflow for generating gene-level counts from RNA-seq data

This Snakemake workflow implements the pre-processing steps to achieve gene-level read counts from raw RNA-seq data.

## Features

- Quality reports (`FastQC`)
- Downloading/indexing of genome and annotations
- Trimming (`fastp`)
- Merging of sequencing units (e.g. samples split across multiple lanes)
- Single- or paired-end read compatibility
- Alignment (`STAR` 2-pass mode)
- Gene-level counting (`Subread` `featureCounts`)

## Standardised usage

Standardised usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=baerlachlan/smk-rnaseq-star-featurecounts).

However, Snakemake standardised usage requires internet access which is commonly unavailable in an HPC environment.
If the intention is to run the workflow in an offline environment, please see [Recommended usage](#recommended-usage).

## Recommended usage

For compatibility across environments, the source code of this workflow is available via [Releases](https://github.com/baerlachlan/smk-rnaseq-star-featurecounts/releases).

1. Download and extract the workflow's [latest release](https://github.com/baerlachlan/smk-rnaseq-gatk-variants/releases/latest)
1. Remove unrequired files
    ```bash
    rm .gitignore .snakemake-workflow-catalog.yml
    ```
1. Follow the instructions in [`config/README.md`](config/README.md) to modify `config/samples.tsv` and `config/units.tsv`
1. Follow the comments in `config/config.yaml` to configure the workflow parameters
1. Use the example profile in `workflow/profiles/example/config.v8+.yaml` as a guide to fine-tune workflow-specific resource configuration
    - NOTE: the example profile has been designed for compatibility with my [SLURM profile](https://github.com/baerlachlan/smk-cluster-generic-slurm)
    - To automatically enable a workflow-specific profile, move it to `workflow/profiles/default`, for example:
        ```bash
        cp -r workflow/profiles/example workflow/profiles/default
        ```
1. Execute the workflow
    ```bash
    snakemake
    ```

## Testing

Example data and configurations are available for testing of this workflow.

Backup the example configuration, and copy the test configuration in its place: 

```bash
## Backup example config
mv config/ config_example/
## Test paired-end
cp -r .test/config_pe/ config/
## Or test single-end
cp -r .test/config_se/ config/
```

The example data is small, so the test workflow profile can be used upon execution.
To keep intermediate files, specify the `--notemp` flag.

```bash
snakemake --notemp --workflow-profile workflow/profiles/test
```
