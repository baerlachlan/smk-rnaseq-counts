rule genome_coverage:
    input:
        bam_inputs(),
    output:
        "results/coverage/{SAMPLE}.bedGraph"
    params:
        "-bg"
    wrapper:
        "v7.2.0/bio/bedtools/genomecov"