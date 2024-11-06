# Snakemake workflow for estimating read counts from RNA-seq data

This Snakemake workflow implements the pre-processing steps to estimate gene- and transcript-level read counts from raw RNA-seq data.

## Contents

- [Workflow summary](#workflow-summary)
- [Standardised usage](#standardised-usage)
- [Recommended usage](#recommended-usage)
- [Testing](#testing)

## Workflow summary

### 1. Raw data

Raw RNA-seq data is expected in `FASTQ` format.
No specific location is required for the raw data.
This is left flexible for the user, however the chosen location must be specified in [`config/units.tsv`](config/units.tsv).

### 2. Trim

Trimming is performed with `fastp`

Module: [`trim.smk`](workflow/rules/trim.smk)

- Input: Raw data at chosen location (`FASTQ`)
- Output: `results/trim` (`FASTQ`)

### 3. Merge (optional)

Merging involves concatenating files of the same sample.
This step will only be performed when multiple sequencing units are specified for the same sample identifier in [`config/units.tsv`](config/units.tsv).
A common example for multiple sequencing units is when a sample is split across multiple lanes.

Module: [`merge.smk`](workflow/rules/merge.smk)

- Input: `results/trim` (`FASTQ`)
- Output: `results/merge` (`FASTQ`)

### 4. Align

Aligment to the genome is performed with `STAR` in 2-pass mode.

Module: [`align.smk`](workflow/rules/align.smk)

- Input: `results/trim` and/or `results/merge` (`FASTQ`)
- Output: `results/align` (`BAM`)

### 5. Gene-level counts

Summarisation of read counts to the gene-level is performed with `featureCounts`.

Module: [`featureCounts.smk`](workflow/rules/featureCounts.smk)

- Input: `results/align` (`BAM`)
- Output: `results/featureCounts` (`TSV`)
  - NOTE: by default, the workflow will run featurecounts on all samples 3 times. It will produce 3 folders as subdirectories of `results/featureCounts`: `unstranded`, `stranded` and `reverse`. This is for the purpose of inferring strandedness from the count summary stats. If the strandedness of the library is already known, this can be specified in [`config/config.yaml`](config/config.yaml), and only a single subdirectory will be produced.

### 6. Transcript-level counts

Summarisation of read counts to the transcript-level is performed with `Salmon`.

Module: [`salmon.smk`](workflow/rules/salmon.smk)

- Input: `results/trim` and/or `results/merge` (`FASTQ`)
- Output: `results/salmon` (`One directory per sample`)

### Other features

#### Quality control

Quality reports are produced with `FastQC` for raw, trimmed and aligned data.

Module: [`fastqc.smk`](workflow/rules/fastqc.smk)

#### Reference files & indexing

Genome, transcriptome and annotation files are downloaded from Ensembl.
The user may provide their own reference files by copying them into the `resources/`, with filenames `genome.fa`, `transcriptome.fa` and `annotation.gtf` respectively.
The workflow produces reference indices for both `STAR` and `Salmon` as required.

Module: [`refs.smk`](workflow/rules/refs.smk)

#### Single and paired end data compatibility

Both single and paired end data is compatible with this workflow.
This is specified by how one configures the [`config/units.tsv`](config/units.tsv) file (see [`config/README.md`](config/README.md)).

## Standardised usage

Standardised usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=baerlachlan/smk-rnaseq-star-featurecounts).

However, Snakemake standardised usage requires internet access which is commonly unavailable in an HPC environment.
If the intention is to run the workflow in an offline environment, please see [Recommended usage](#recommended-usage).

## Recommended usage

For compatibility across environments, the source code of this workflow is available via [Releases](https://github.com/baerlachlan/smk-rnaseq-star-featurecounts/releases).

1. Download and extract the workflow's [latest release](https://github.com/baerlachlan/smk-rnaseq-star-featurecounts/releases/latest)
2. Follow the instructions in [`config/README.md`](config/README.md) to modify [`config/samples.tsv`](config/samples.tsv) and [`config/units.tsv`](config/units.tsv)
3. Follow the comments in [`config/config.yaml`](config/config.yaml) to configure the workflow parameters
4. Use the [example profile](workflow/profiles/default/config.v8+.yaml) as a guide to fine-tune workflow-specific resource configuration
    - NOTE: the example profile has been designed for compatibility with my [SLURM profile](https://github.com/baerlachlan/smk-cluster-generic-slurm)
5. Execute the workflow
    ```bash
    snakemake
    ```

## Testing

Example data and configurations are available in the `.test` directory for testing this workflow.
The example data is small, so the test workflow profile can be used upon execution.
To keep and examine intermediary files, specify the `--notemp` flag.

```bash
## Test paired-end
snakemake --configfile .test/config_pe/config.yaml --workflow-profile workflow/profiles/test --notemp
## or single-end
snakemake --configfile .test/config_pe/config.yaml --workflow-profile workflow/profiles/test --notemp
```
