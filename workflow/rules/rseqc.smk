rule read_distribution:
    input:
        aln=bam_inputs(),
        refgene=annotation_bed,
    output:
        "results/rseqc/read_distribution/{SAMPLE}.read_distribution.txt"
    wrapper:
        "v7.2.0/bio/rseqc/read_distribution"


rule inner_distance:
    input:
        aln=bam_inputs(),
        refgene=annotation_bed,
    output:
        reads_inner_distance="results/rseqc/inner_distance/{SAMPLE}.inner_distance.txt",
        freq="results/rseqc/inner_distance/{SAMPLE}.inner_distance_freq.txt",
        pdf="results/rseqc/inner_distance/{SAMPLE}.inner_distance_plot.pdf",
        plot_r="results/rseqc/inner_distance/{SAMPLE}.inner_distance_plot.r",
    conda:
        "../envs/rseqc.yml"
    shell:
        """
        inner_distance.py -i {input.aln} -o results/rseqc/inner_distance/{wildcards.SAMPLE} -r {input.refgene}
        """
