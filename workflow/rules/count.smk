## Unstranded
rule featureCounts_s0:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/unstranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=0,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v4.0.0/bio/subread/featurecounts"


## Stranded
rule featureCounts_s1:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/stranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=1,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v4.0.0/bio/subread/featurecounts"


## Reverse-stranded
rule featureCounts_s2:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/reverse/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=2,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v4.0.0/bio/subread/featurecounts"
