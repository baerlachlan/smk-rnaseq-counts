rule salmon_quant:
    input:
        unpack(salmon_inputs),
        index=salmon_index_dir,
    output:
        quant="results/salmon/{SAMPLE}/quant.sf",
        lib="results/salmon/{SAMPLE}/lib_format_counts.json",
    params:
        libtype=config["salmon"]["quant"]["libtype"],
        extra=config["salmon"]["quant"]["extra"],
    wrapper:
        "v7.2.0/bio/salmon/quant"
