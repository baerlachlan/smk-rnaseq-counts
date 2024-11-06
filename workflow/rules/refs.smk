rule genome_get:
    output:
        "resources/genome.fa",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v4.0.0/bio/reference/ensembl-sequence"


rule transcriptome_get:
    output:
        "resources/transcriptome.fa",
    params:
        species=config["ref"]["species"],
        datatype="cdna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v4.0.0/bio/reference/ensembl-sequence"


rule annotation_get:
    output:
        "resources/annotation.gtf",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v4.0.0/bio/reference/ensembl-annotation"


rule star_index:
    input:
        fasta="resources/genome.fa",
        gtf="resources/annotation.gtf",
    output:
        temp(directory("resources/star")),
    params:
        sjdbOverhang=int(config["read_length"]) - 1,
        extra="",
    wrapper:
        "v4.0.0/bio/star/index"


rule salmon_decoy:
    input:
        transcriptome="resources/transcriptome.fa",
        genome="resources/genome.fa",
    output:
        gentrome=temp("resources/gentrome.fa"),
        decoys=temp("resources/decoys.txt"),
    wrapper:
        "v4.0.0/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences="resources/gentrome.fa",
        decoys="resources/decoys.txt",
    output:
        multiext(
            "resources/salmon_index/",
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
        directory("resources/salmon_index/"),  # Added for dependency
    params:
        extra=config["salmon"]["index"]["extra"],
    wrapper:
        "v4.0.0/bio/salmon/index"
