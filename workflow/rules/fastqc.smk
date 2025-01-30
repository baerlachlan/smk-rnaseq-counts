rule fastqc_raw:
    input:
        fastqc_raw_inputs,
    output:
        html="results/raw_data/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.html",
        zip="results/raw_data/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.zip",
    params:
        extra="--noextract",
    wrapper:
        "v5.5.2/bio/fastqc"


rule fastqc_trim:
    input:
        fastqc_trim_inputs,
    output:
        html="results/trim/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.html",
        zip="results/trim/FastQC/{SAMPLE}_{UNIT}_{PAIRTAG}_fastqc.zip",
    params:
        extra="--noextract",
    wrapper:
        "v5.5.2/bio/fastqc"


rule fastqc_align:
    input:
        "results/align/bam/{SAMPLE}.bam",
    output:
        html="results/align/FastQC/{SAMPLE}_fastqc.html",
        zip="results/align/FastQC/{SAMPLE}_fastqc.zip",
    params:
        extra="--noextract",
    wrapper:
        "v5.5.2/bio/fastqc"
