rule genome_get:
    output:
        genome_fa,
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v5.5.2/bio/reference/ensembl-sequence"


rule transcriptome_get:
    output:
        transcriptome_fa,
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v5.5.2/bio/reference/ensembl-sequence"


rule annotation_get:
    output:
        annotation_gtf,
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v5.5.2/bio/reference/ensembl-annotation"


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
        "v5.5.2/bio/star/index"


rule salmon_decoy:
    input:
        transcriptome=transcriptome_fa,
        genome=genome_fa,
    output:
        gentrome=temp(gentrome_fa),
        decoys=temp(decoys_txt),
    wrapper:
        "v5.5.2/bio/salmon/decoys"


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
        "v5.5.2/bio/salmon/index"
