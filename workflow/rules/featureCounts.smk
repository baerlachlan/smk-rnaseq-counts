## Unstranded
rule featureCounts_s0:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/unstranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=0,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v5.5.2/bio/subread/featurecounts"


## Stranded
rule featureCounts_s1:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/stranded/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=1,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v5.5.2/bio/subread/featurecounts"


## Reverse-stranded
rule featureCounts_s2:
    input:
        unpack(featureCounts_inputs),
        annotation=annotation_gtf,
    output:
        multiext(
            "results/featureCounts/reverse/all",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    params:
        strand=2,
        extra=config["featureCounts"]["extra"],
    wrapper:
        "v5.5.2/bio/subread/featurecounts"
