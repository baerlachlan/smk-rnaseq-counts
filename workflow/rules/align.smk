rule align:
    input:
        unpack(align_inputs),
        idx="resources/star",
    output:
        aln="results/align/bam/{SAMPLE}.bam",
        log="results/align/log/{SAMPLE}.log",
        log_final="results/align/log/{SAMPLE}.log.final.out",
    params:
        extra=f"--sjdbOverhang {int(config["read_length"]) - 1} --outSAMtype BAM SortedByCoordinate --twopassMode Basic",
    wrapper:
        "v4.0.0/bio/star/align"


rule align_index:
	input:
		"results/align/bam/{SAMPLE}.bam",
	output:
		"results/align/bam/{SAMPLE}.bam.bai",
	threads: 8
	resources:
		runtime=60,
		mem_mb=32768,
	wrapper:
		"v4.0.0/bio/samtools/index"


rule align_md5:
    input:
        align_md5_inputs(),
    output:
        "results/align/bam/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
