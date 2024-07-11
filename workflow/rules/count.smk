rule count:
    input:
        samples=expand("results/align/bam/{SAMPLE}.bam", SAMPLE=samples["sample"]),
        annotation="resources/annotation.gtf",
    output:
        multiext(
            "results/count/all",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    params:
        strand=config["count"]["strandedness"],
        extra=config["count"]["extra"],
    wrapper:
        "v3.7.0/bio/subread/featurecounts"
