rule deduplicate:
    input:
        bam="results/align/bam/{SAMPLE}.bam" if config["align"]["keep_bam"] else temp("results/align/bam/{SAMPLE}.bam"),
        bai="results/align/bam/{SAMPLE}.bam.bai",
    output:
        bam="results/deduplicate/bam/{SAMPLE}.bam",
        log="results/deduplicate/log/{SAMPLE}.log",
    params:
        extra=config["deduplicate"]["extra"],
    conda:
        "../envs/umitools.yml"
    shell:
        """
        umi_tools dedup --stdin={input.bam} --stdout={output.bam} \
        --log={output.log} {params.extra}
        """

rule deduplicate_index:
    input:
        "results/deduplicate/bam/{SAMPLE}.bam",
    output:
        "results/deduplicate/bam/{SAMPLE}.bam.bai",
    wrapper:
        "v5.5.2/bio/samtools/index"