rule align_m5C:
    input:
        unpack(get_clean_fq),
        index="resources/rnam5c_hisat2_genome_index",
    output:
        mapping_result="{project}/m5C/{sample}/{sample}.bam",
        multimapping_result="{project}/m5C/{sample}/{sample}.multimappers.bam",
        config=temp("{project}/m5C/{sample}/{sample}.config"),
        tmp_rev=temp("{project}/m5C/{sample}/{sample}_rev.fastq"),
        tmp_forward=temp("{project}/m5C/{sample}/{sample}_R2.fastq"),
    log:
        "logs/{project}/align_m5C/{sample}.log",
    benchmark:
        "benchmarks/align/m5C/{project}_{sample}.benchmark.txt"
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    params:
        prefix="{project}/m5C/{sample}/{sample}"
    threads: config["threads"]["m5C_align"]
    shell:
        "echo '-p {threads}' > {output.config} &&"
        "gzip -dc {input.fq1} > {output.tmp_rev} && "
        "gzip -dc {input.fq2} > {output.tmp_forward} && "
        "python /opt/conda/RNA-m5C/2_m5C_step-by-step_hisat2/BS_hisat2.py "
        "-F {output.tmp_forward} "
        "-R {output.tmp_rev} "
        "--del-convert "
        "--del-sam "
        "--hisat2-path /opt/conda/bin "
        "-I {input.index} "
        "-o {params.prefix} "
        "--continue-prefix {params.prefix}_tmp "
        "--index-prefix HISAT2 " 
        "--hisat2-param {output.config} 2>{log}"
