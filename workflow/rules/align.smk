rule align:
    input:
        unpack(align_inputs),
        idx="resources/genome",
    output:
        aln=temp("results/align/bam/{SAMPLE}.bam"),
        log="results/align/log/{SAMPLE}.log",
        log_final="results/align/log/{SAMPLE}.log.final.out",
    params:
        extra=f"--sjdbOverhang {int(config["read_length"]) - 1} --outSAMtype BAM SortedByCoordinate --twopassMode Basic",
    wrapper:
        "v3.7.0/bio/star/align"


rule align_md5:
    input:
        align_md5_inputs(),
    output:
        "results/align/bam/md5.txt",
    shell:
        """
        md5sum {input} > {output}
        """
