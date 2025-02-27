rule merge:
    input:
        unpack(merge_inputs),
    output:
        fq=temp("results/merge/fastq/{SAMPLE}_{PAIRTAG}.fastq.gz"),
    shell:
        """
        cat {input.fq} > {output.fq}
        """
