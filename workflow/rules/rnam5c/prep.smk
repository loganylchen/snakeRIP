rule un_compressed:
    input:
        unpack(get_clean_fq),
    output:
        tmp_rev=temp("{project}/m5C/{sample}/{sample}_rev.fastq"),
        tmp_forward=temp("{project}/m5C/{sample}/{sample}_forward.fastq"),
    log:
        "logs/{project}/prep_m5C_uncompressed/{sample}.log",
    benchmark:
        "benchmarks/prep_uncompressed/m5C/{project}_{sample}.benchmark.txt"
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    params:
        prefix="{project}/m5C/{sample}/{sample}"
    threads: config["threads"]["m5C_align"]
    shell:
        "gzip -dc {input.fq1} > {output.tmp_rev} && "
        "gzip -dc {input.fq2} > {output.tmp_forward}  2>{log}"