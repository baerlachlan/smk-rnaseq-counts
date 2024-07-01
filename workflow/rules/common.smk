import pandas as pd
from snakemake.utils import min_version, validate

min_version("8.14.0")


configfile: "config/config.yaml"


samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample": str}).set_index(
    "sample", drop=False
)
validate(samples, "../schemas/samples.schema.yml")

units = pd.read_csv(
    config["units"], sep="\t", dtype={"sample": str, "unit": str}
).set_index(["sample", "unit"], drop=False)
validate(units, "../schemas/units.schema.yml")

####
## Helper functions
####


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), f"all units for sample {sample} must be single or paired end"
    return all_paired


if all(is_paired_end(i) for i in samples["sample"]):
    pair_tags = ["R1", "R2"]
else:
    pair_tags = ["R0"]

####
## Wildcard constraints
####


wildcard_constraints:
    SAMPLE="|".join(samples["sample"]),
    UNIT="|".join(units["unit"]),
    PAIRTAG="|".join(pair_tags),


####
## Input functions
####


def fastqc_raw_inputs(wildcards):
    unit = units.loc[wildcards.SAMPLE, wildcards.UNIT]
    if is_paired_end(wildcards.SAMPLE):
        if wildcards.PAIRTAG == pair_tags[0]:
            return f"{unit.fq1}"
        elif wildcards.PAIRTAG == pair_tags[1]:
            return f"{unit.fq2}"
    else:
        return f"{unit.fq1}"


def fastqc_trim_inputs(wildcards):
    if is_paired_end(wildcards.SAMPLE):
        if wildcards.PAIRTAG == pair_tags[0]:
            return "results/trim/fastq/{SAMPLE}_{UNIT}_R1.fastq.gz"
        elif wildcards.PAIRTAG == pair_tags[1]:
            return "results/trim/fastq/{SAMPLE}_{UNIT}_R2.fastq.gz"
    else:
        return "results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz"


def trim_inputs(wildcards):
    unit = units.loc[wildcards.SAMPLE, wildcards.UNIT]
    if is_paired_end(wildcards.SAMPLE):
        return {"sample": [f"{unit.fq1}", f"{unit.fq2}"]}
    else:
        return {"sample": [f"{unit.fq1}"]}


def merge_inputs(wildcards):
    sample_units = units.loc[wildcards.SAMPLE]
    return {
        "fq": expand(
            "results/trim/fastq/{{SAMPLE}}_{UNIT}_{{PAIRTAG}}.fastq.gz",
            UNIT=sample_units["unit"],
        )
    }


def align_inputs(wildcards):
    if is_paired_end(wildcards.SAMPLE):
        return {
            "fq1": "results/merge/fastq/{SAMPLE}_R1.fastq.gz",
            "fq2": "results/merge/fastq/{SAMPLE}_R2.fastq.gz",
        }
    else:
        return {"fq1": "results/merge/fastq/{SAMPLE}_R0.fastq.gz"}


def trim_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/trim/fastq/{SAMPLE}_{UNIT}_{PAIRTAG}.fastq.gz",
                SAMPLE=sample_units["sample"],
                UNIT=sample_units["unit"],
                PAIRTAG=pair_tags,
            )
        )
    return inputs


def merge_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/merge/fastq/{SAMPLE}_{PAIRTAG}.fastq.gz",
                SAMPLE=sample_units["sample"],
                PAIRTAG=pair_tags,
            )
        )
    return inputs


def align_md5_inputs():
    inputs = []
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        inputs.extend(
            expand(
                "results/align/bam/{SAMPLE}.bam",
                SAMPLE=sample_units["sample"],
            )
        )
    return inputs


####
## Workflow output files (Rule all inputs)
####


def workflow_outputs():
    """
    Returns all file endpoints for the workflow
    """

    outputs = []

    ## FastQC outputs
    for sample in samples["sample"]:
        sample_units = units.loc[sample]
        ## Raw
        outputs.extend(
            expand(
                "results/raw_data/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.{EXT}",
                SAMPLE=sample,
                UNIT=sample_units["unit"],
                PAIRTAG=pair_tags,
                EXT=["html", "zip"],
            )
        )
        ## Trim
        outputs.extend(
            expand(
                "results/trim/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.{EXT}",
                SAMPLE=sample,
                UNIT=sample_units["unit"],
                PAIRTAG=pair_tags,
                EXT=["html", "zip"],
            )
        )
        ## Align
        outputs.extend(
            expand(
                "results/align/FastQC/{SAMPLE}_fastqc.{EXT}",
                SAMPLE=sample,
                EXT=["html", "zip"],
            )
        )

    ## md5sums
    outputs.extend(
        [
            "results/trim/fastq/md5.txt",
            "results/merge/fastq/md5.txt",
            "results/align/bam/md5.txt",
        ]
    )

    ## Processed counts
    outputs.append("results/count/all.featureCounts")

    return outputs
