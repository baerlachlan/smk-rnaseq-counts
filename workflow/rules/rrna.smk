rule rrna_get:
    output:
        "resources/rdna.fa"
    log:
        "logs/rrna_get/rrna_get.log",
    shell:
        """
        wget -O {output} \
            "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=555853&extrafeat=null&conwithfeat=on&hide-cdd=on&ncbi_phid=CE89B83286F157810000000000B500AC" 2> {log}
        """


rule rrna_index:
    input:
        fasta="resources/rdna.fa",
    output:
        directory("resources/rrna_index/"),
    log:
        "logs/rrna_index/rrna_index.log",
    params:
        extra="--genomeSAindexNbases 6",
    wrapper:
        "v7.2.0/bio/star/index"


rule rrna_align:
    input:
        unpack(align_inputs),
        idx="resources/rrna_index/",
    output:
        aln=temp("results/rrna/bam/{SAMPLE}.unsorted.bam"),
        log="results/rrna/log/{SAMPLE}.log",
        log_final="results/rrna/log/{SAMPLE}.log.final.out",
    log:
        "logs/rrna_align/{SAMPLE}.log",
    params:
        extra=config["align"]["extra"],
    wrapper:
        "v7.2.0/bio/star/align"


rule rrna_align_sort:
    input:
        "results/rrna/bam/{SAMPLE}.unsorted.bam",
    output:
        "results/rrna/bam/{SAMPLE}.bam" if config["align"]["keep_bam"] else temp("results/rrna/bam/{SAMPLE}.bam")
    log:
        "logs/rrna_align_sort/{SAMPLE}.log",
    wrapper:
        "v7.2.0/bio/samtools/sort"


rule rrna_align_index:
    input:
        "results/rrna/bam/{SAMPLE}.bam",
    output:
        "results/rrna/bam/{SAMPLE}.bam.bai",
    log:
        "logs/rrna_align_index/{SAMPLE}.log",
    wrapper:
        "v7.2.0/bio/samtools/index"
