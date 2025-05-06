rule align_m5C_genome:
    input:
        tmp_rev="{project}/m5C/{sample}/{sample}_rev.fastq",
        tmp_forward="{project}/m5C/{sample}/{sample}_forward.fastq",
        index="resources/rnam5c_hisat2_genome_index",
    output:
        mapping_result="{project}/m5C/{sample}/{sample}_genome.bam",
        multimapping_result="{project}/m5C/{sample}/{sample}_genome.multimappers.bam",
        config=temp("{project}/m5C/{sample}/{sample}.hisat2.config"),
    log:
        "logs/{project}/align_m5C/{sample}.log",
    benchmark:
        "benchmarks/align/m5C/{project}_{sample}.benchmark.txt"
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    params:
        prefix="{project}/m5C/{sample}/{sample}_genome"
    threads: config["threads"]["m5C_align"]
    shell:
        "echo '-p {threads}' > {output.config} &&"
        "python /opt/conda/RNA-m5C/2_m5C_step-by-step_hisat2/BS_hisat2.py "
        "-F {input.tmp_forward} "
        "-R {input.tmp_rev} "
        "--del-convert "
        "--del-sam "
        "--hisat2-path /opt/conda/bin "
        "-I {input.index} "
        "-o {params.prefix} "
        "--continue-prefix {params.prefix}_tmp "
        "--index-prefix HISAT2 " 
        "--hisat2-param {output.config} 2>{log}"



rule align_m5C_transcriptome:
    input:
        tmp_rev="{project}/m5C/{sample}/{sample}_rev.fastq",
        tmp_forward="{project}/m5C/{sample}/{sample}_forward.fastq",
        c2t_index=multiext("resources/transcriptome.c2t",'.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2','.rev.2.bt2' ),
        gene_list="resources/genelist.txt",
        anno="resources/genome.anno",
        genomesize="resources/genomesize.txt",
    output:
        mapping_result="{project}/m5C/{sample}/{sample}_transcriptome.bam",
        coverted_mapping_result="{project}/m5C/{sample}/{sample}_transcriptome.converted.bam",
        multimapping_result="{project}/m5C/{sample}/{sample}_transcriptome.multimappers.bam",
        config=temp("{project}/m5C/{sample}/{sample}.bowtie2.config"),
    log:
        "logs/{project}/align_m5C_transcriptome/{sample}.log",
    benchmark:
        "benchmarks/align_transcriptome/m5C/{project}_{sample}.benchmark.txt"
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    params:
        index_prefix='resources/transcriptome.c2t',
        prefix="{project}/m5C/{sample}/{sample}_transcriptome"
    threads: config["threads"]["m5C_align"]
    shell:
        "echo '-p {threads}' > {output.config} &&"
        "python /opt/conda/RNA-m5C/3_m5C_step-by-step_bowtie2/BS_bowtie2.py "
        "-F {input.tmp_forward} "
        "-R {input.tmp_rev} "
        "--del-convert "
        "--del-sam "
        "--bowtie2-path /opt/conda/bin "
        "-I {input.index} "
        "-g {input.gene_list}"
        "-o {params.prefix} "
        "--continue-prefix {params.prefix}_tmp "
        "--bowtie2-param {output.config} 2>{log} && "
        "python /opt/conda/RNA-m5C/3_m5C_step-by-step_bowtie2/Bam_transcriptome_to_genome_v1.0.py "
        "-i {output.mapping_result} "
        "-o {output.coverted_mapping_result} "
        "-a {input.anno} "
        "--dict {input.genomesize} 2>>{log} "

rule merge_bams:
    input:
        genome_bam="{project}/m5C/{sample}/{sample}_genome.bam",
        transcriptome_bam="{project}/m5C/{sample}/{sample}_transcriptome.converted.bam",
    output:
        merged_bam="{project}/m5C/{sample}/{sample}.merged.bam",
    log:
        "logs/{project}/merge_bams/{sample}.log",
    benchmark:
        "benchmarks/merge_bams/m5C/{project}_{sample}.benchmark.txt"
    container:
        "docker://btrspg/rnam5c:4c6656b36e5f88116a5a2df8c23897891cc887f5"
    threads:
        config["threads"]["m5C_align"]
    shell:
        "python /opt/conda/RNA-m5C/4_m5C_step-by-step_pileup/concat_bam.py "
        "-t {threads} "
        "-i {input.genome_bam} {input.transcriptome_bam} "
        "-o {output.merged_bam} "
        "--sort --index  2>{log}"

