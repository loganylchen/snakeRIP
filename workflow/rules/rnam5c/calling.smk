rule m5c_calling:
    input:
        genome_reference="resources/genome.fasta",
        merged_bam="{project}/m5C/{sample}/{sample}.merged.bam",
        noredundance_base="resources/genome.nd_base",
    output:
        m5c_pileup_tmp="{project}/m5C/{sample}/{sample}.pileup.tmp",
        m5c_calling_result="{project}/m5C/{sample}/{sample}.m5c.3.txt",
    log:
        "logs/{project}/m5c_calling/{sample}.log",
    benchmark:
        "benchmarks/m5c_calling/m5C/{project}_{sample}.benchmark.txt",
    threads: config["threads"]["m5C_align"],
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    params:
        output_file_prefix="{project}/m5C/{sample}/{sample}.m5c"
    script:
        "../../scripts/m5c_calling.sh"
