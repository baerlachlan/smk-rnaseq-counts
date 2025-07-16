rule genome_get:
    output:
        genome_fa,
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.2.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome_fa,
    output:
        genome_fai,
    params:
        extra="",
    wrapper:
        "v7.2.0/bio/samtools/faidx"


rule genome_chrom_sizes:
    input:
        genome_fai,
    output:
        genome_chrom_sizes,
    shell:
        """
        cut -f1,2 {input} | sort -k1,1 > {output}
        """


rule transcriptome_get:
    output:
        transcriptome_fa,
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.2.0/bio/reference/ensembl-sequence"


rule annotation_get:
    output:
        annotation_gtf,
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v7.2.0/bio/reference/ensembl-annotation"


rule annotation_sort:
    input:
        annotation_gtf
    output:
        annotation_sorted,
    shell:
        """
        cat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' > {output}
        """


rule annotation_genePred:
    input:
        annotation_gtf,
    output:
        temp(annotation_genePred),
    params:
        extra="-genePredExt",
    wrapper:
        "v7.2.0/bio/ucsc/gtfToGenePred"


rule annotation_bed:
    input:
        annotation_genePred
    output:
        temp(annotation_bed)
    params:
        extra="",
    wrapper:
        "v7.2.0/bio/ucsc/genePredToBed"


rule annotation_intergenic:
    input:
        gtf=annotation_sorted,
        chromsizes=genome_chrom_sizes,
    output:
        temp(annotation_intergenic),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} {{print $1, $4-1, $5}}' {input.gtf} | \
            bedtools merge -i - | \
            bedtools complement -i - -g {input.chromsizes} > {output}
        """


rule annotation_exon:
    input:
        annotation_sorted,
    output:
        temp(annotation_exon),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk '$3 == "exon"' {input} | \
            awk 'BEGIN{{OFS="\t"}} {{print $1, $4-1, $5}}' | \
            bedtools merge -i - > {output}
        """


rule annotation_intron:
    input:
        annotation_exon=annotation_exon,
        annotation_intergenic=annotation_intergenic,
        chromsizes=genome_chrom_sizes,
    output:
        temp(annotation_intron),
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input.annotation_exon} {input.annotation_intergenic} | \
            sort -k1,1 -k2,2n | \
            bedtools complement -i - -g {input.chromsizes} | \
            bedtools merge -i - > {output}
        """


rule star_index:
    input:
        fasta=genome_fa,
        gtf=annotation_gtf,
    output:
        directory(star_index_dir),
    params:
        sjdbOverhang=int(config["read_length"]) - 1,
        extra="",
    wrapper:
        "v7.2.0/bio/star/index"


rule salmon_decoy:
    input:
        transcriptome=transcriptome_fa,
        genome=genome_fa,
    output:
        gentrome=temp(gentrome_fa),
        decoys=temp(decoys_txt),
    wrapper:
        "v7.2.0/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences=gentrome_fa,
        decoys=decoys_txt,
    output:
        multiext(
            salmon_index_dir,
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json",
        ),
        directory(salmon_index_dir),  # Added for dependency
    params:
        extra=config["salmon"]["index"]["extra"],
    wrapper:
        "v7.2.0/bio/salmon/index"
