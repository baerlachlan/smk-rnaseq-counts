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
        ####
        ## Header
        ####
        echo -e "chromosome\tregion\ttotal_covered_bases\tmean_per_covered_base" > {output}.tmp
        ####
        ## Whole genome
        ####
        ## All chromosomes
        awk '{{cov += ($3 - $2) * $4; len += ($3 - $2)}} END {{print cov, cov / (len + 1)}}' {input.bedGraph} | \
            awk '{{print "All" "\t" "genome" "\t" $1 "\t" $2}}' >> {output}.tmp
        ## Per chromosome
        awk '{{cov[$1] += ($3 - $2) * $4; len[$1] += ($3 - $2)}} END {{for (chr in cov) print chr, cov[chr], cov[chr] / len[chr]}}' {input.bedGraph} | \
            awk '{{print $1 "\t" "genome" "\t" $2 "\t" $3}}' >> {output}.tmp
        ####
        ## Exonic regions
        ####
        ## All chromosomes
        bedtools intersect -a {input.bedGraph} -b {input.exons} | \
            awk '{{cov += ($3 - $2) * $4; len += ($3 - $2)}} END {{print cov, cov / (len + 1)}}' | \
            awk '{{print "All" "\t" "exons" "\t" $1 "\t" $2}}' >> {output}.tmp
        ## Per chromosome
        bedtools intersect -a {input.bedGraph} -b {input.exons} | \
            awk '{{cov[$1] += ($3 - $2) * $4; len[$1] += ($3 - $2)}} END {{for (chr in cov) print chr, cov[chr], cov[chr] / len[chr]}}' | \
            awk '{{print $1 "\t" "exons" "\t" $2 "\t" $3}}' >> {output}.tmp
        ####
        ## Intronic regions
        ####
        ## All chromosomes
        bedtools intersect -a {input.bedGraph} -b {input.introns} | \
            awk '{{cov += ($3 - $2) * $4; len += ($3 - $2)}} END {{print cov, cov / (len + 1)}}' | \
            awk '{{print "All" "\t" "introns" "\t" $1 "\t" $2}}' >> {output}.tmp
        ## Per chromosome
        bedtools intersect -a {input.bedGraph} -b {input.introns} | \
            awk '{{cov[$1] += ($3 - $2) * $4; len[$1] += ($3 - $2)}} END {{for (chr in cov) print chr, cov[chr], cov[chr] / len[chr]}}' | \
            awk '{{print $1 "\t" "introns" "\t" $2 "\t" $3}}' >> {output}.tmp
        ####
        ## Intergenic regions
        ####
        ## All chromosomes
        bedtools intersect -a {input.bedGraph} -b {input.intergenic} | \
            awk '{{cov += ($3 - $2) * $4; len += ($3 - $2)}} END {{print cov, cov / (len + 1)}}' | \
            awk '{{print "All" "\t" "intergenic" "\t" $1 "\t" $2}}' >> {output}.tmp
        ## Per chromosome
        bedtools intersect -a {input.bedGraph} -b {input.intergenic} | \
            awk '{{cov[$1] += ($3 - $2) * $4; len[$1] += ($3 - $2)}} END {{for (chr in cov) print chr, cov[chr], cov[chr] / len[chr]}}' | \
            awk '{{print $1 "\t" "intergenic" "\t" $2 "\t" $3}}' >> {output}.tmp
        ####
        ## Final sort
        ####
        awk 'NR==1{{print $0; next}} {{print $0 | "sort -k1,1n"}}' {output}.tmp > {output}
        rm {output}.tmp
        """
