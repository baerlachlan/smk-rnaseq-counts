rule align:
    input:
        unpack(align_inputs),
        idx=star_index_dir,
    output:
        aln=temp("results/align/bam/{SAMPLE}.unsorted.bam"),
        log="results/align/log/{SAMPLE}.log",
        log_final="results/align/log/{SAMPLE}.log.final.out",
    params:
        extra=config["align"]["extra"],
    wrapper:
        "v7.2.0/bio/star/align"


rule align_sort:
    input:
        "results/align/bam/{SAMPLE}.unsorted.bam",
    output:
        "results/align/bam/{SAMPLE}.bam" if config["align"]["keep_bam"] else temp("results/align/bam/{SAMPLE}.bam")
    wrapper:
        "v7.2.0/bio/samtools/sort"


rule align_index:
    input:
        "results/align/bam/{SAMPLE}.bam",
    output:
        "results/align/bam/{SAMPLE}.bam.bai",
    wrapper:
        "v7.2.0/bio/samtools/index"
