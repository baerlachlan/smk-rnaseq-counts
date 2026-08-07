## Unstranded
rule featureCounts_s0:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/{BAM_SOURCE}/unstranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    log:
        "logs/featureCounts_s0/{BAM_SOURCE}.log",
    params:
        strand=0,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


## Stranded
rule featureCounts_s1:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/{BAM_SOURCE}/stranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    log:
        "logs/featureCounts_s1/{BAM_SOURCE}.log",
    params:
        strand=1,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v7.2.0/bio/subread/featurecounts"


## Reverse-stranded
rule featureCounts_s2:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/{BAM_SOURCE}/reverse/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    log:
        "logs/featureCounts_s2/{BAM_SOURCE}.log",
    params:
        strand=2,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v7.2.0/bio/subread/featurecounts"
