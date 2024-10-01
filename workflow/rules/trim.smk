rule trim_se:
    input:
        unpack(trim_inputs),
    output:
        trimmed=temp(["results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz"]),
        html="results/trim/log/{SAMPLE}_{UNIT}.html",
        json="results/trim/log/{SAMPLE}_{UNIT}.json",
    params:
        extra="--qualified_quality_phred 20 --length_required 35 --trim_poly_g",
    wrapper:
        "v4.0.0/bio/fastp"


rule trim_pe:
    input:
        unpack(trim_inputs),
    output:
        trimmed=temp(
            [
                "results/trim/fastq/{SAMPLE}_{UNIT}_R1.fastq.gz",
                "results/trim/fastq/{SAMPLE}_{UNIT}_R2.fastq.gz",
            ]
        ),
        html="results/trim/log/{SAMPLE}_{UNIT}.html",
        json="results/trim/log/{SAMPLE}_{UNIT}.json",
    params:
        extra="--detect_adapter_for_pe --qualified_quality_phred 20 --length_required 35 --trim_poly_g",
    wrapper:
        "v4.0.0/bio/fastp"


rule trim_md5:
    input:
        trim_md5_inputs(),
    output:
        "results/trim/fastq/md5.txt",
    conda:
        "../envs/parallel.yml"
    shell:
        """
        echo "{input}" | tr " " "\n" | parallel -j {threads} md5sum  > {output}
        """
