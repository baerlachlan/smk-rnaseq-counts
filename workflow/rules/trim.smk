rule trim_se:
    input:
        unpack(trim_inputs),
    output:
        trimmed=temp(["results/trim/fastq/{SAMPLE}_{UNIT}_R0.fastq.gz"]),
        html="results/trim/log/{SAMPLE}_{UNIT}.html",
        json="results/trim/log/{SAMPLE}_{UNIT}.json",
    log:
        "logs/trim_se/{SAMPLE}_{UNIT}.log",
    params:
        extra=fastp_args,
    wrapper:
        "v7.2.0/bio/fastp"


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
    log:
        "logs/trim_pe/{SAMPLE}_{UNIT}.log",
    params:
        extra=fastp_args,
    wrapper:
        "v7.2.0/bio/fastp"
