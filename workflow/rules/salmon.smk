rule salmon_quant:
    input:
        unpack(salmon_inputs),
        index="resources/salmon_index/"
    output:
        quant="results/salmon/{SAMPLE}/quant.sf",
        lib="results/salmon/{SAMPLE}/lib_format_counts.json",
    params:
        libtype="A",
        extra=config["salmon"]["quant"]["extra"],
    wrapper:
        "v4.0.0/bio/salmon/quant"
