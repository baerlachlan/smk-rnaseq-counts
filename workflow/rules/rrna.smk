rule rrna_get:
    output:
        "resources/rdna.fa"
    shell:
        """
        wget -O {output} \
            "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=555853&extrafeat=null&conwithfeat=on&hide-cdd=on&ncbi_phid=CE89B83286F157810000000000B500AC"
        """


rule rrna_index:
    input:
        fasta="resources/rdna.fa",
    output:
        directory("resources/rrna_index/"),
    params:
        extra="",
    wrapper:
        "v7.2.0/bio/star/index"


rule rrna_align:
    input:
        unpack(align_inputs),
        idx="resources/rrna_index/",
    output:
        aln="results/rrna/bam/{SAMPLE}.bam" if config["align"]["keep_bam"] else temp("results/rrna/bam/{SAMPLE}.bam"),
        log="results/rrna/log/{SAMPLE}.log",
        log_final="results/rrna/log/{SAMPLE}.log.final.out",
    params:
        extra=f"{config["align"]["extra"]}",
    wrapper:
        "v5.5.2/bio/star/align"