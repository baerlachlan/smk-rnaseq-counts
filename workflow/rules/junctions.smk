rule junctions_count:
    input:
        unpack(junctions_inputs),
    output:
        "results/junctions/{SAMPLE}.bed",
    params:
        extra=config["junctions"]["regtools_extra"],
    conda:
        "../envs/regtools.yml"
    shell:
        """
        regtools junctions extract {params.extra} {input.bam} -o {output}
        """


## Regtools reports the chromStart and chromEnd (column 2 & 3),
## which corresponds to the start and end positions of the junction anchor,
## i.e. the maximum exonic read overhang.
## To get the actual junction coordinates, the blocksizes in column 11 need
## to be added/subtracted
rule junctions_adjust:
    input:
        "results/junctions/{SAMPLE}.bed"
    output:
        "results/junctions/{SAMPLE}.adj.bed"
    shell:
        """
        awk '{{split($11, sizes, ","); print $1, $2+sizes[1], $3-sizes[2], $4, $5, $6}}' OFS="\t" {input} > {output}
        """