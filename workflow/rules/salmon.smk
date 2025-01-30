rule salmon_quant:
    input:
        unpack(salmon_inputs),
        index="resources/salmon_index/",
    output:
        quant="results/salmon/{SAMPLE}/quant.sf",
        lib="results/salmon/{SAMPLE}/lib_format_counts.json",
    params:
        libtype=config["salmon"]["quant"]["libtype"],
        extra=config["salmon"]["quant"]["extra"],
    wrapper:
        "v5.5.2/bio/salmon/quant"
