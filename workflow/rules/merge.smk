rule merge:
    input:
        unpack(merge_inputs),
    output:
        fq=temp("results/merge/fastq/{SAMPLE}_{PAIRTAG}.fastq.gz"),
    log:
        "logs/merge/{SAMPLE}_{PAIRTAG}.log",
    shell:
        """
        cat {input.fq} > {output.fq} 2> {log}
        """
