# Configuration

## Workflow config

The workflow requires configuration by modification of `config/config.yaml`.
Follow the explanations provided as comments in the file.

## Sample & unit config

The configuration of samples and units is specified as tab-separated value (`.tsv`) files.
Each `.tsv` requires specific columns (see below), but extra columns may be present (however, will not be used).

### `samples.tsv`

The default path for the sample sheet is `config/samples.tsv`.
This may be changed via configuration in `config/config.yaml`.

`samples.tsv` requires only one column named `sample`, which contains the desired names of the samples.
Sample names must be unique, corresponding to a physical sample.
Biological and technical replicates should be specified as separate samples.

### `units.tsv`

The default path for the unit sheet is `config/units.tsv`.
This may be changed via configuration in `config/config.yaml`.

`units.tsv` requires four columns, named `sample`, `unit`, `fq1` and `fq2`.
Each row of the units sheet corresponds to a single sequencing unit.
Therefore, for each sample specified in `samples.tsv`, one or more sequencing units should be present.
`unit` values must be unique within each sample.
A common example of an experiment with multiple sequencing units is a sample split across several runs/lanes.

For each unit, the respective path to `FASTQ` files must be specified in the `fq1` and `fq2` columns.
Both columns must exist, however, the `fq2` column may be left empty in the case of single-end sequencing experiments.
This is how one specifies whether single- or paired-end rules are run by the workflow.
