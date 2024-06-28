rule merge:
    input:
        unpack(merge_inputs),
    output:
        fq=temp("results/merge/fastq/{SAMPLE}_{PAIRTAG}.fastq.gz"),
    shell:
        """
        cat {input.fq} > {output.fq}
        """


rule merge_md5:
    input:
        merge_md5_inputs(),
    output:
        "results/merge/fastq/md5.txt",
    shell:
        """
        md5sum {input} > {output}
        """
