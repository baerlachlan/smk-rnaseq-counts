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
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
