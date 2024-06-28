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
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra="-Q 1 --minOverlap 35 --fracOverlap 0.9 -J",
    wrapper:
        "v3.7.0/bio/subread/featurecounts"
