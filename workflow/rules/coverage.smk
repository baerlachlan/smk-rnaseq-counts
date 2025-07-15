rule coverage_bedGraph:
    input:
        bam_inputs(),
    output:
        "results/coverage/{SAMPLE}.bedGraph" if config["coverage"]["keep_bedGraphs"] else temp("results/coverage/{SAMPLE}.bedGraph")
    params:
        "-bg"
    wrapper:
        "v7.2.0/bio/bedtools/genomecov"


rule coverage_summary:
    input:
        bedGraph="results/coverage/{SAMPLE}.bedGraph",
        exons=annotation_exon,
        introns=annotation_intron,
        intergenic=annotation_intergenic,
    output:
        "results/coverage/{SAMPLE}.coverage.summary"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        ## Header
        echo -e "region\ttotal_covered_bases" > {output}
        ## Whole genome
        cov_genome=$(awk '{{sum += ($3 - $2) * $4}} END {{print sum}}' {input.bedGraph})
        echo -e "genome\t$cov_genome" >> {output}
        ## Exonic regions
        cov_exonic=$(bedtools intersect -a {input.bedGraph} -b {input.exons} -wa | \
            awk '{{sum += ($3 - $2) * $4}} END {{print sum}}')
        echo -e "exons\t$cov_exonic" >> {output}
        ## Intronic regions
        cov_intronic=$(bedtools intersect -a {input.bedGraph} -b {input.introns} -wa | \
            awk '{{sum += ($3 - $2) * $4}} END {{print sum}}')
        echo -e "introns\t$cov_intronic" >> {output}
        ## Intergenic regions
        cov_intergenic=$(bedtools intersect -a {input.bedGraph} -b {input.intergenic} -wa | \
            awk '{{sum += ($3 - $2) * $4}} END {{print sum}}')
        echo -e "intergenic\t$cov_intergenic" >> {output}
        ## Chromosomes
        awk '{{cov[$1] += ($3 - $2) * $4}} END {{for (chr in cov) print chr, cov[chr]}}' {input.bedGraph} | \
            awk '{{print $1 "\t" $2}}' >> {output}
        """

# rule coverage_bigWig:
#     input:
#         bedgraph="results/coverage/{SAMPLE}.bedGraph",
#         chromsizes=genome_chrom_sizes,
#     output:
#         "results/coverage/{SAMPLE}.bw"
#     params:
#         ""
#     wrapper:
#         "v7.2.0/bio/ucsc/bedGraphToBigWig"
