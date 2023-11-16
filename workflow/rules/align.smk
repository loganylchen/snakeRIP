rule align:
    input:
        unpack(get_clean_data),
    output:
        aln="results/star/{sample}/{sample}.star.bam",
    log:
        "logs/star/{sample}.log",
    benchmark:
        "benchmarks/{sample}.star_align.benchmark.txt"
    params:
        idx=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: config["threads"]["star"]
    wrapper:
        "v1.21.4/bio/star/align"
