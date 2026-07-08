# Snakemake Workflow For RNA-Seq Counts

This workflow processes raw RNA-seq FASTQ files into gene-level and transcript-level count outputs for downstream analysis in R or similar environments.

The standard workflow performs raw/processed read QC, read trimming, optional merging of multiple sequencing units, STAR genome alignment, featureCounts gene-level quantification, and Salmon transcript-level quantification.
Optional modules support UMI-based deduplication, splice junction extraction, RSeQC summaries, rRNA alignment checks, and genome coverage summaries.

## Contents

- [Workflow Summary](#workflow-summary)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Testing](#testing)
- [Outputs](#outputs)
- [Optional Modules](#optional-modules)
- [Reference Files](#reference-files)
- [HPC And Profiles](#hpc-and-profiles)
- [Important Assumptions](#important-assumptions)

## Workflow Summary

The main processing steps are:

1. Raw FASTQ quality control with FastQC.
2. Read trimming with fastp.
3. Optional merging of multiple sequencing units from the same sample.
4. Genome alignment with STAR.
5. BAM sorting and indexing with samtools.
6. Optional UMI-based deduplication with UMI-tools.
7. Gene-level quantification with featureCounts.
8. Transcript-level quantification with Salmon.
9. Optional downstream QC and diagnostic modules.

The workflow supports paired-end and single-end RNA-seq data, but a single run must not mix paired-end and single-end samples.

## Quick Start

Activate an environment containing Snakemake before running the workflow.

```bash
conda activate snakemake
```

For a new project, edit these files first:

```text
config/config.yaml
config/samples.tsv
config/units.tsv
```

Then perform a dry-run:

```bash
snakemake -n
```

The default `config/samples.tsv` and `config/units.tsv` contain placeholder sample names and FASTQ filenames.
These are intended as templates and will not run against real data unless matching files exist under the configured `data_dir`.

## Configuration

The workflow is configured with three main files:

| File | Purpose |
| --- | --- |
| `config/config.yaml` | Main workflow settings, reference settings, tool options, and optional module activation. |
| `config/samples.tsv` | Sample names to process. |
| `config/units.tsv` | FASTQ files and sequencing units associated with each sample. |

See [`config/README.md`](config/README.md) for detailed examples and field descriptions.

## Testing

Small paired-end and single-end test configurations are provided in `.test/`.
These can be run without editing the default config files.

Dry-run the paired-end test workflow:

```bash
snakemake -n \
    --configfile .test/config_pe/config.yaml \
    --workflow-profile workflow/profiles/test
```

Dry-run the single-end test workflow:

```bash
snakemake -n \
    --configfile .test/config_se/config.yaml \
    --workflow-profile workflow/profiles/test
```

To run the test workflow and keep temporary/intermediate outputs for inspection, omit `-n` and add `--notemp`:

```bash
snakemake \
    --configfile .test/config_pe/config.yaml \
    --workflow-profile workflow/profiles/test \
    --notemp
```

## Outputs

The exact outputs depend on which modules are activated in `config/config.yaml`.

| Module | Main outputs |
| --- | --- |
| FastQC | `results/raw_data/FastQC/`, `results/trim/FastQC/`, `results/align/FastQC/` |
| fastp trimming | `results/trim/fastq/`, `results/trim/log/` |
| Merge | `results/merge/fastq/` |
| STAR genome alignment | `results/align/bam/{sample}.bam`, `results/align/bam/{sample}.bam.bai`, `results/align/log/` |
| UMI deduplication | `results/deduplicate/bam/{sample}.bam`, `results/deduplicate/bam/{sample}.bam.bai`, `results/deduplicate/log/` |
| featureCounts | `results/featureCounts/{unstranded,stranded,reverse}/all.featureCounts` |
| Salmon | `results/salmon/{sample}/quant.sf` |
| Junctions | `results/junctions/{sample}.adj.bed` |
| RSeQC read distribution | `results/rseqc/read_distribution/{sample}.read_distribution.txt` |
| RSeQC inner distance | `results/rseqc/inner_distance/{sample}.inner_distance.txt` |
| rRNA alignment | `results/rrna/bam/{sample}.bam`, `results/rrna/bam/{sample}.bam.bai` |
| Coverage | `results/coverage/{sample}.coverage.summary` |

By default, featureCounts runs for all three strandedness settings: unstranded, stranded, and reverse-stranded.
This is useful for inferring library strandedness from the assignment summaries.
If strandedness is already known, set only the appropriate value in `config/config.yaml`.

## Optional Modules

Optional modules are controlled in `config/config.yaml`.

| Config section | Purpose |
| --- | --- |
| `deduplicate.activate` | Deduplicate aligned reads using UMI-tools. |
| `junctions.activate` | Extract splice junctions using regtools. |
| `read_distribution.activate` | Run RSeQC read distribution. |
| `inner_distance.activate` | Run RSeQC inner distance analysis. |
| `rrna.activate` | Align reads to an rDNA/rRNA reference and index the resulting BAM files. |
| `coverage.activate` | Produce coverage summaries across genome, exon, intron, and intergenic regions. |

When deduplication is activated, downstream BAM-consuming modules use deduplicated BAM files where supported.

## Reference Files

Reference genome, transcriptome, and annotation files are downloaded from Ensembl using the species, release, and build specified in `config/config.yaml`.
The workflow builds STAR and Salmon indices as required.

Custom or spike-in sequences can be merged into the Ensembl reference using the `ref.merge_with` settings in `config/config.yaml`.
This requires both a FASTA file and a matching GTF file.

Derived annotation BED files are created for optional RSeQC and coverage modules.
Intergenic regions are defined as the complement of annotated gene intervals, and intronic regions are defined as gene regions outside merged exons.

## HPC And Profiles

Workflow profiles provide rule-specific thread and resource settings.
Snakemake will automatically use the workflow profile in `workflow/profiles/default/` when running from this repository.

The default profile is:

```text
workflow/profiles/default/config.v8+.yaml
```

This workflow is run on the University of Adelaide HPC system, Phoenix.
For that type of SLURM-based setup, an external Snakemake profile is available at <https://github.com/baerlachlan/smk-slurm-profile>.
Users should treat that profile as a starting point and modify it for their own scheduler, account, partition, memory limits, temporary storage, and runtime requirements.

The test profile is:

```text
workflow/profiles/test/config.v8+.yaml
```

It is intended for the small example data under `.test/`.

## Important Assumptions

- All samples in a single workflow run must be either paired-end or single-end.
- `fq1` and `fq2` values in `config/units.tsv` are resolved relative to `data_dir` in `config/config.yaml`.
- Multiple rows for the same sample in `config/units.tsv` are treated as multiple sequencing units and are merged after trimming.
- `featureCounts.extra` defaults to paired-end settings and should be reviewed for single-end projects.
- `trim.extra`, `featureCounts.extra`, and `deduplicate.extra` are passed directly to their respective tools and should be reviewed for each dataset.
- Sample names should be unique, simple identifiers without whitespace.
- Optional modules may require additional reference-derived files and will add jobs to the DAG when activated.

## Standardised Usage

This workflow can also be described through the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/).
For HPC systems without internet access, downloading a release or cloning the repository and running with local profiles is usually more practical.
