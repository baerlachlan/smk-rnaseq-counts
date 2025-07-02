import sys
import pandas as pd
import os
import math
from snakemake.utils import min_version, validate


min_version("8.14.0")


####
## Config
####


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
## Helper functions and code
####


species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]


genome_fa = f"resources/{species.capitalize()}.{build}.dna.primary_assembly.fa"
transcriptome_fa = f"resources/{species.capitalize()}.{build}.cdna.all.fa"
annotation_gtf = f"resources/{species.capitalize()}.{build}.{str(release)}.gtf"
star_index_dir = "resources/star_index/"
gentrome_fa = f"resources/{species.capitalize()}.{build}.gentrome.fa"
decoys_txt = f"resources/{species.capitalize()}.{build}.decoys.txt"
salmon_index_dir = "resources/salmon_index/"


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), f"ERROR: All units for sample {sample} must be single or paired end"
    return all_paired


paired_end_samples = [is_paired_end(i) for i in samples["sample"]]
single_end_samples = [not i for i in paired_end_samples]


assert all(paired_end_samples) or all(
    single_end_samples
), "ERROR: The dataset must be entirely paired or single end, not a combination."


if all(paired_end_samples):
    pair_tags = ["R1", "R2"]
elif all(single_end_samples):
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
            return os.path.join(config["data_dir"], f"{unit.fq1}")
        elif wildcards.PAIRTAG == pair_tags[1]:
            return os.path.join(config["data_dir"], f"{unit.fq2}")
    else:
        return os.path.join(config["data_dir"], f"{unit.fq1}")


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
        return {"sample": [os.path.join(config["data_dir"], f"{unit.fq1}"), os.path.join(config["data_dir"], f"{unit.fq2}")]}
    else:
        return {"sample": [os.path.join(config["data_dir"], f"{unit.fq1}")]}


def merge_inputs(wildcards):
    sample_units = units.loc[wildcards.SAMPLE]
    return {
        "fq": expand(
            "results/trim/fastq/{{SAMPLE}}_{UNIT}_{{PAIRTAG}}.fastq.gz",
            UNIT=sample_units["unit"],
        )
    }


def align_inputs(wildcards):
    sample_units = units.loc[wildcards.SAMPLE]
    if is_paired_end(wildcards.SAMPLE):
        if len(sample_units) == 1:
            return {
                "fq1": expand(
                    "results/trim/fastq/{{SAMPLE}}_{UNIT}_R1.fastq.gz",
                    UNIT=sample_units["unit"],
                ),
                "fq2": expand(
                    "results/trim/fastq/{{SAMPLE}}_{UNIT}_R2.fastq.gz",
                    UNIT=sample_units["unit"],
                ),
            }
        else:
            return {
                "fq1": "results/merge/fastq/{SAMPLE}_R1.fastq.gz",
                "fq2": "results/merge/fastq/{SAMPLE}_R2.fastq.gz",
            }
    else:
        if len(sample_units) == 1:
            return {
                "fq1": expand(
                    "results/trim/fastq/{{SAMPLE}}_{UNIT}_R0.fastq.gz",
                    UNIT=sample_units["unit"],
                )
            }
        else:
            return {"fq1": "results/merge/fastq/{SAMPLE}_R0.fastq.gz"}


def featureCounts_inputs(wildcards):
    if config["deduplicate"]["activate"]:
        return {
            "samples": expand(
                "results/deduplicate/bam/{SAMPLE}.bam",
                SAMPLE=samples["sample"]
            ),
            "bai": expand(
                "results/deduplicate/bam/{SAMPLE}.bam.bai",
                SAMPLE=samples["sample"]
            )
        }
    else:
        return {
            "samples": expand(
                "results/align/bam/{SAMPLE}.bam",
                SAMPLE=samples["sample"]
            ),
            "bai": expand(
                "results/align/bam/{SAMPLE}.bam.bai",
                SAMPLE=samples["sample"]
            )
        }



def salmon_inputs(wildcards):
    sample_units = units.loc[wildcards.SAMPLE]
    if is_paired_end(wildcards.SAMPLE):
        if len(sample_units) == 1:
            return {
                "r1": expand(
                    "results/trim/fastq/{SAMPLE}_{UNIT}_R1.fastq.gz",
                    SAMPLE=sample_units["sample"],
                    UNIT=sample_units["unit"],
                ),
                "r2": expand(
                    "results/trim/fastq/{SAMPLE}_{UNIT}_R2.fastq.gz",
                    SAMPLE=sample_units["sample"],
                    UNIT=sample_units["unit"],
                ),
            }
        else:
            return {
                "r1": "results/merge/fastq/{SAMPLE}_R1.fastq.gz",
                "r2": "results/merge/fastq/{SAMPLE}_R2.fastq.gz",
            }
    else:
        if len(sample_units) == 1:
            return {
                "r": expand(
                    "results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz",
                    SAMPLE=sample_units["sample"],
                    UNIT=sample_units["unit"],
                ),
            }
        else:
            return {
                "r": "results/merge/fastq/{SAMPLE}_R0.fastq.gz",
            }


####
## Workflow output files (Rule all inputs)
####


def workflow_outputs():
    """
    Returns all file endpoints for the workflow
    """

    outputs = []

    ## FastQC outputs
    if config["fastqc"]["activate"]:
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

    ## Aligned reads
    outputs.extend(expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]))
    outputs.extend(expand("results/align/bam/{SAMPLE}.bam.bai", SAMPLE=samples["sample"]))

    ## Gene-level counts (featureCounts)
    if config["featureCounts"]["activate"]:
        strandedness_labels = ["unstranded", "stranded", "reverse"]
        for i in config["featureCounts"]["strandedness"]:
            outputs.append(
                f"results/featureCounts/{strandedness_labels[i]}/all.featureCounts"
            )

    ## Transcript-level counts (Salmon)
    if config["salmon"]["activate"]:
        outputs.extend(expand("results/salmon/{SAMPLE}/quant.sf", SAMPLE=samples["sample"]))

    return outputs
