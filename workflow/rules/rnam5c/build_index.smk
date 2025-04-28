rule build_genome_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        directory('resources/rnam5c_hisat2_genome_index')
    resources:
        mem_mb=250 * 1024
    log:
        "logs/rnam5c/build_genome_index.log",
    container:
        "docker://btrspg/rnam5c:913a09dee6d2d414e9ca0d63c84755c12f943a82"
    benchmark:
        "benchmarks/rnam5c_build_genome_index.benchmark.txt"
    threads: config['threads']['rnam5c_build_genome_index']
    shell:
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/BS_hisat2_index.py "
        "-i {input.fasta} "
        "--gtf {input.gtf} "
        "--hisat2-path /opt/conda/bin "
        "-o {output} 2>{log} 1>&2 "

rule build_transcriptome_index:
    input:
        fasta="resources/transcriptome.fasta",
    output:
        c2t_fa="resources/transcriptome.c2t.fa",
        c2t_index=multiext("resources/transcriptome.c2t",'.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2','.rev.2.bt2' )
    params:
        prefix="resources/transcriptome.c2t"
    resources:
        mem_mb=30 * 1024
    log:
        "logs/rnam5c/build_transcriptome_index.log",
    container:
        "docker://btrspg/rnam5c:913a09dee6d2d414e9ca0d63c84755c12f943a82"
    benchmark:
        "benchmarks/rnam5c_build_transcriptome_index.benchmark.txt"
    shell:
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/fasta_c2t.py "
        "-i {input.fasta} "
        "> {output.c2t_fa} 2>{log} && "
        "bowtie2-build "
        "{output.c2t_fa} "
        "{params.prefix} 2>>{log} 1>> {log}"

rule get_metadata:
    input:
        genome_fasta="resources/genome.fasta",
        transcriptome_fasta="resources/transcriptome.c2t.fa",
        gtf="resources/genome.gtf",
    output:
        anno="resources/genome.anno",
        base="resources/genome.base",
        nd_base="resources/genome.nd_base",
        genelist="resources/genelist.txt",
        genomesize="resources/genomesize.txt",
    resources:
        mem_mb=250 * 1024
    log:
        "logs/rnam5c/get_metadata.log",
    container:
        "docker://btrspg/rnam5c:913a09dee6d2d414e9ca0d63c84755c12f943a82"
    benchmark:
        "benchmarks/rnam5c_get_metadata.benchmark.txt"
    shell:
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/gtf2anno.py "
        "-i {input.gtf} > {output.anno} && "
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/gtf2genelist.py "
        "-i {input.gtf} "
        "-f {input.transcriptome_fasta} > {output.genelist} && "
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/ref_sizes.py "
        "-i {input.genome_fasta} -o {output.genomesize} &&"
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/anno_to_base.py "
        "-i {output.anno} -o {output.base} && "
        "python /opt/conda/RNA-m5C/0_m5C_step-by-step_metadata/anno_to_base_remove_redundance_v1.0.py "
        "-i {output.base} "
        "-o {output.nd_base} "
        "-g {output.genelist} 2>{log}"


