rule genome_get:
    output:
        temp("resources/genome.fa"),
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    wrapper:
        "v3.7.0/bio/reference/ensembl-sequence"


rule genome_index:
    input:
        "resources/genome.fa",
    output:
        temp("resources/genome.fa.fai"),
    params:
        extra="",
    wrapper:
        "v3.7.0/bio/samtools/faidx"


rule annotation_get:
    output:
        temp("resources/annotation.gtf"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    wrapper:
        "v3.7.0/bio/reference/ensembl-annotation"


rule star_index:
    input:
        fasta="resources/genome.fa",
        gtf="resources/annotation.gtf",
    output:
        temp(directory("resources/genome")),
    params:
        sjdbOverhang=int(config["read_length"]) - 1,
        extra="",
    wrapper:
        "v3.7.0/bio/star/index"
