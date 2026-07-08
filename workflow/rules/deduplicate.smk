rule deduplicate:
    input:
        bam="results/align/bam/{SAMPLE}.bam",
        bai="results/align/bam/{SAMPLE}.bam.bai",
    output:
        bam="results/deduplicate/bam/{SAMPLE}.bam" if config["deduplicate"]["keep_bam"] else temp("results/deduplicate/bam/{SAMPLE}.bam"),
        tool_log="results/deduplicate/log/{SAMPLE}.log",
    log:
        "logs/deduplicate/{SAMPLE}.log",
    params:
        extra=config["deduplicate"]["extra"],
    conda:
        "../envs/umitools.yml"
    shell:
        """
        umi_tools dedup --stdin={input.bam} --stdout={output.bam} \
        --log={output.tool_log} {params.extra} 2> {log}
        """


rule deduplicate_index:
    input:
        "results/deduplicate/bam/{SAMPLE}.bam",
    output:
        "results/deduplicate/bam/{SAMPLE}.bam.bai",
    log:
        "logs/deduplicate_index/{SAMPLE}.log",
    wrapper:
        "v7.2.0/bio/samtools/index"
