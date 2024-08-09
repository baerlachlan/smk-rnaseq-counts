## Unstranded
rule count_s0:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/unstranded/all",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    params:
        strand=0,
        extra=config["count"]["extra"],
    wrapper:
        "v3.7.0/bio/subread/featurecounts"

## Stranded
rule count_s1:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/stranded/all",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    params:
        strand=1,
        extra=config["count"]["extra"],
    wrapper:
        "v3.7.0/bio/subread/featurecounts"

## Reverse-stranded
rule count_s2:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/reverse/all",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    params:
        strand=2,
        extra=config["count"]["extra"],
    wrapper:
        "v3.7.0/bio/subread/featurecounts"
