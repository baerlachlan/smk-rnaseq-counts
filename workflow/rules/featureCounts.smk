## Unstranded
rule featureCounts_s0:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        bai=expand("results/align/bam/{SAMPLE}.bam.bai", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/featureCounts/unstranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=0,
        extra=config["featureCounts"]["extra"],
    resources:
        slurm_extra=lambda wc, input: f"'--gres=tmpfs:{math.ceil((input.size_mb / 1024) * 2)}G'"
    wrapper:
        "v5.5.2/bio/subread/featurecounts"


## Stranded
rule featureCounts_s1:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        bai=expand("results/align/bam/{SAMPLE}.bam.bai", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/featureCounts/stranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=1,
        extra=config["featureCounts"]["extra"],
    resources:
        slurm_extra=lambda wc, input: f"'--gres=tmpfs:{math.ceil((input.size_mb / 1024) * 2)}G'"
    wrapper:
        "v5.5.2/bio/subread/featurecounts"


## Reverse-stranded
rule featureCounts_s2:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        bai=expand("results/align/bam/{SAMPLE}.bam.bai", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/featureCounts/reverse/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=2,
        extra=config["featureCounts"]["extra"],
    resources:
        slurm_extra=lambda wc, input: f"'--gres=tmpfs:{math.ceil((input.size_mb / 1024) * 2)}G'"
    wrapper:
        "v5.5.2/bio/subread/featurecounts"
