# Configuration

This directory contains the user-editable configuration files for the workflow.

Most projects require editing:

```text
config/config.yaml
config/samples.tsv
config/units.tsv
```

## `config.yaml`

`config.yaml` controls workflow inputs, references, tool settings, and optional modules.

| Section | Purpose |
| --- | --- |
| `samples` | Path to the sample table. Defaults to `config/samples.tsv`. |
| `units` | Path to the sequencing unit table. Defaults to `config/units.tsv`. |
| `data_dir` | Directory containing FASTQ files listed in `units.tsv`. |
| `ref` | Ensembl reference species, release, build, and optional custom sequence merging. |
| `read_length` | Raw read length, used for STAR splice junction index settings. |
| `fastqc` | FastQC activation and options. |
| `trim` | fastp trimming options and UMI trimming defaults. |
| `align` | STAR genome alignment options. |
| `featureCounts` | Gene-level counting options. |
| `salmon` | Salmon indexing and quantification options. |
| `junctions` | Optional regtools junction extraction. |
| `deduplicate` | Optional UMI-tools deduplication. |
| `read_distribution` | Optional RSeQC read distribution. |
| `inner_distance` | Optional RSeQC inner distance. |
| `rrna` | Optional rRNA/rDNA alignment. |
| `coverage` | Optional genome/exon/intron/intergenic coverage summaries. |

The comments in `config.yaml` describe the main options in place.
Always review paired-end/single-end-specific tool options before running a new dataset.

## `samples.tsv`

The sample table lists the biological or analysis-level samples to process.

Required columns:

| Column | Description |
| --- | --- |
| `sample` | Unique sample identifier. |

Example:

```tsv
sample
WT_1
WT_2
KO_1
KO_2
```

Sample names should be unique and should not contain whitespace.
Biological replicates should be listed as separate samples.
Technical sequencing units for the same sample should be represented in `units.tsv`, not as separate samples.

## `units.tsv`

The unit table links samples to FASTQ files.
Each row corresponds to one sequencing unit, such as one lane or one sequencing run for a sample.

Required columns:

| Column | Description |
| --- | --- |
| `sample` | Sample identifier matching a row in `samples.tsv`. |
| `unit` | Sequencing unit identifier, unique within each sample. |
| `fq1` | FASTQ file for read 1, or the only FASTQ file for single-end data. |
| `fq2` | FASTQ file for read 2. Leave blank for single-end data. |

Optional UMI override columns:

| Column | Description |
| --- | --- |
| `umi_trim` | Override whether fastp UMI processing is used for this unit. |
| `umi_loc` | Override UMI location for this unit. |
| `umi_len` | Override UMI length for this unit. |
| `umi_skip` | Override number of bases removed after the UMI for this unit. |

Per-unit `umi_len` and `umi_skip` values should be written as integers.

FASTQ paths in `fq1` and `fq2` are resolved relative to `data_dir` in `config.yaml`.

For example, with:

```yaml
data_dir: "../raw_data"
```

this unit entry:

```tsv
sample	unit	fq1	fq2
WT_1	1	WT_1_R1.fastq.gz	WT_1_R2.fastq.gz
```

refers to:

```text
../raw_data/WT_1_R1.fastq.gz
../raw_data/WT_1_R2.fastq.gz
```

### Paired-End Example

One sequencing unit per sample:

```tsv
sample	unit	fq1	fq2
WT_1	1	WT_1_R1.fastq.gz	WT_1_R2.fastq.gz
WT_2	1	WT_2_R1.fastq.gz	WT_2_R2.fastq.gz
KO_1	1	KO_1_R1.fastq.gz	KO_1_R2.fastq.gz
KO_2	1	KO_2_R1.fastq.gz	KO_2_R2.fastq.gz
```

Multiple sequencing units for the same sample:

```tsv
sample	unit	fq1	fq2
WT_1	L001	WT_1_L001_R1.fastq.gz	WT_1_L001_R2.fastq.gz
WT_1	L002	WT_1_L002_R1.fastq.gz	WT_1_L002_R2.fastq.gz
WT_2	L001	WT_2_L001_R1.fastq.gz	WT_2_L001_R2.fastq.gz
WT_2	L002	WT_2_L002_R1.fastq.gz	WT_2_L002_R2.fastq.gz
```

Multiple units for the same sample are trimmed separately and then merged before alignment and Salmon quantification.

### Single-End Example

For single-end data, keep the `fq2` column but leave it blank:

```tsv
sample	unit	fq1	fq2
WT_1	1	WT_1.fastq.gz	
WT_2	1	WT_2.fastq.gz	
KO_1	1	KO_1.fastq.gz	
KO_2	1	KO_2.fastq.gz	
```

All samples in one run must be either paired-end or single-end.
Do not mix paired-end and single-end samples in the same config.

## References

Reference files are configured under `ref` in `config.yaml`.

Example:

```yaml
ref:
  species: homo_sapiens
  release: 112
  build: GRCh38
```

The workflow downloads the matching Ensembl genome FASTA, transcriptome FASTA, and GTF annotation, then builds STAR and Salmon indices as needed.

To merge custom or spike-in sequences into the reference, activate `ref.merge_with` and provide matching FASTA and GTF files:

```yaml
ref:
  merge_with:
    activate: True
    fasta: "resources/genome_to_merge.fa"
    gtf: "resources/annotation_to_merge.gtf"
```

The FASTA and GTF should describe the same additional sequences.

## Trimming And UMIs

Global fastp options are configured under `trim.extra`.

UMI trimming defaults are configured under `trim.umi`:

```yaml
trim:
  umi:
    activate: False
    umi_loc: "per_read"
    umi_len: "5"
    umi_skip: "2"
```

These defaults can be overridden for individual rows in `units.tsv` with the optional UMI columns described above.

If UMI deduplication is required after alignment, also activate the `deduplicate` module.

## featureCounts Strandedness

featureCounts strandedness is configured as a list:

```yaml
featureCounts:
  strandedness: [0, 1, 2]
```

Values are:

| Value | Meaning | Output directory |
| --- | --- | --- |
| `0` | Unstranded | `results/featureCounts/{align,deduplicate}/unstranded/` |
| `1` | Stranded | `results/featureCounts/{align,deduplicate}/stranded/` |
| `2` | Reverse-stranded | `results/featureCounts/{align,deduplicate}/reverse/` |

Using all three values can help infer library strandedness from featureCounts summary statistics.
If strandedness is known, provide only the relevant value.
Aligned featureCounts outputs are always generated when featureCounts is active.
If `deduplicate.activate` is `True`, featureCounts also generates a separate count set from the deduplicated BAMs.

For example:

```yaml
featureCounts:
  strandedness: [2]
```

## Optional Modules

Optional modules are disabled by default unless otherwise specified in `config.yaml`.

### UMI Deduplication

```yaml
deduplicate:
  activate: True
```

Deduplicates aligned BAM files using UMI-tools.
When active, featureCounts produces both aligned and deduplicated count sets.
Other supported downstream modules use deduplicated BAM files.

### Splice Junctions

```yaml
junctions:
  activate: True
```

Extracts splice junctions from aligned BAM files using regtools and writes adjusted BED files under `results/junctions/`.

### RSeQC

```yaml
read_distribution:
  activate: True

inner_distance:
  activate: True
```

Runs selected RSeQC summaries using annotation-derived BED files.

### rRNA Alignment

```yaml
rrna:
  activate: True
```

Aligns reads to the bundled rDNA/rRNA reference with STAR and writes sorted, indexed BAM files under `results/rrna/bam/`.

### Coverage

```yaml
coverage:
  activate: True
```

Computes coverage summaries across genome, exon, intron, and intergenic regions.
Intergenic regions are defined as the complement of annotated gene intervals.
Intronic regions are defined as gene regions outside merged exons.

## Common Checks Before Running

- Confirm every `sample` in `units.tsv` appears in `samples.tsv`.
- Confirm each `(sample, unit)` pair is unique.
- Confirm all FASTQ paths resolve correctly relative to `data_dir`.
- Confirm the whole dataset is either paired-end or single-end.
- Review `trim.extra`, `featureCounts.extra`, and `deduplicate.extra` for paired-end/single-end compatibility.
- Run a dry-run with `snakemake -n` before submitting a full workflow.
