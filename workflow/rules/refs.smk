rule genome_get:
    output:
        temp("resources/genome.fa") if config["ref"]["merge_with"]["activate"] else genome_fa,
    log:
        "logs/genome_get/genome_get.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.2.0/bio/reference/ensembl-sequence"


rule genome_merge:
    input:
        genome="resources/genome.fa",
        merge=config["ref"]["merge_with"]["fasta"],
    output:
        genome_fa,
    log:
        "logs/genome_merge/genome_merge.log",
    shell:
        """
        cat {input.genome} {input.merge} > {output} 2> {log}
        """


rule genome_faidx:
    input:
        genome_fa,
    output:
        genome_fai,
    log:
        "logs/genome_faidx/genome_faidx.log",
    params:
        extra="",
    wrapper:
        "v7.2.0/bio/samtools/faidx"


rule genome_chrom_sizes:
    input:
        genome_fai,
    output:
        genome_chrom_sizes,
    log:
        "logs/genome_chrom_sizes/genome_chrom_sizes.log",
    shell:
        """
        cut -f1,2 {input} | sort -k1,1 > {output} 2> {log}
        """


rule transcriptome_get:
    output:
        temp("resources/transcriptome.fa") if config["ref"]["merge_with"]["activate"] else transcriptome_fa,
    log:
        "logs/transcriptome_get/transcriptome_get.log",
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v7.2.0/bio/reference/ensembl-sequence"


rule transcriptome_fasta:
    input:
        fasta=config["ref"]["merge_with"]["fasta"],
        annotation=config["ref"]["merge_with"]["gtf"],
    output:
        transcript_fasta=temp("resources/transcriptome_to_merge.fa"),
    log:
        "logs/transcriptome_fasta/transcriptome_fasta.log",
    params:
        fasta_flag="-w",
        extra="",
    wrapper:
        "v7.2.0/bio/gffread"


rule transcriptome_merge:
    input:
        transcriptome="resources/transcriptome.fa",
        merge="resources/transcriptome_to_merge.fa",
    output:
        transcriptome_fa,
    log:
        "logs/transcriptome_merge/transcriptome_merge.log",
    shell:
        """
        cat {input.transcriptome} {input.merge} > {output} 2> {log}
        """


rule annotation_get:
    output:
        temp("resources/annotation.gtf") if config["ref"]["merge_with"]["activate"] else annotation_gtf,
    log:
        "logs/annotation_get/annotation_get.log",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v7.2.0/bio/reference/ensembl-annotation"


rule annotation_merge:
    input:
        annotation="resources/annotation.gtf",
        merge=config["ref"]["merge_with"]["gtf"],
    output:
        annotation_gtf,
    log:
        "logs/annotation_merge/annotation_merge.log",
    shell:
        """
        cat {input.annotation} {input.merge} > {output} 2> {log}
        """


rule annotation_sort:
    input:
        annotation_gtf
    output:
        annotation_sorted,
    log:
        "logs/annotation_sort/annotation_sort.log",
    shell:
        """
        cat {input} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' > {output} 2> {log}
        """


rule annotation_genePred:
    input:
        annotation_gtf,
    output:
        temp(annotation_genePred),
    log:
        "logs/annotation_genePred/annotation_genePred.log",
    params:
        extra="-genePredExt",
    wrapper:
        "v7.2.0/bio/ucsc/gtfToGenePred"


rule annotation_bed:
    input:
        annotation_genePred
    output:
        temp(annotation_bed)
    log:
        "logs/annotation_bed/annotation_bed.log",
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
    log:
        "logs/annotation_intergenic/annotation_intergenic.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} $1 !~ /^#/ && $3 == "gene" {{print $1, $4-1, $5}}' {input.gtf} | \
            bedtools sort -i - | \
            bedtools merge -i - | \
            bedtools complement -i - -g {input.chromsizes} > {output} 2> {log}
        """


rule annotation_exon:
    input:
        annotation_sorted,
    output:
        temp(annotation_exon),
    log:
        "logs/annotation_exon/annotation_exon.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} $1 !~ /^#/ && $3 == "exon" {{print $1, $4-1, $5}}' {input} | \
            bedtools sort -i - | \
            bedtools merge -i - > {output} 2> {log}
        """


rule annotation_intron:
    input:
        annotation_exon=annotation_exon,
        annotation_intergenic=annotation_intergenic,
        chromsizes=genome_chrom_sizes,
    output:
        temp(annotation_intron),
    log:
        "logs/annotation_intron/annotation_intron.log",
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        cat {input.annotation_exon} {input.annotation_intergenic} | \
            sort -k1,1 -k2,2n | \
            bedtools complement -i - -g {input.chromsizes} | \
            bedtools merge -i - > {output} 2> {log}
        """


rule star_index:
    input:
        fasta=genome_fa,
        gtf=annotation_gtf,
    output:
        directory(star_index_dir),
    log:
        "logs/star_index/star_index.log",
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
    log:
        "logs/salmon_decoy/salmon_decoy.log",
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
    log:
        "logs/salmon_index/salmon_index.log",
    params:
        extra=config["salmon"]["index"]["extra"],
    wrapper:
        "v7.2.0/bio/salmon/index"
