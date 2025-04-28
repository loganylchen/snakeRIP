rule get_genome:
    output:
        "resources/raw_genome.fasta",
    log:
        "logs/ref/get-genome.log",
    container:
        "docker://btrspg/lftp:latest"
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    threads: config['threads']['lftp']
    script:
        "../scripts/get_ensembl_sequence.sh"


rule extract_transcripts:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf"
    output:
        transcripts="resources/transcriptome.fasta"
    log:
        "logs/ref/get_transcripts.log"
    container:
        "docker://btrspg/gffread:0.12.7"
    benchmark:
        "benchmarks/extract_transcripts.benchmark.txt"
    shell:
        "gffread -w {output.transcripts} -g {input.fasta} {input.gtf} 2> {log}"

rule get_annotation:
    output:
        "resources/raw_genome.gtf",
    params:
        species=config["reference"]["species"],
        fmt="gtf",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        flavor="",
    container:
        "docker://btrspg/lftp:latest"
    threads: 1
    log:
        "logs/ref/get_annotation.log",
    benchmark:
        "benchmarks/get_annotation.benchmark.txt"
    script:
        "../scripts/get_ensembl_annotation.sh"


rule filtering_genome_and_annotation:
    input:
        fasta="resources/raw_genome.fasta",
        gtf="resources/raw_genome.gtf",
    output:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    log:
        "logs/ref/filtering_references.log",
    params:
        contigs=config["reference"]["contigs"],
    container:
        "docker://btrspg/biopython:1.85"

    benchmark:
        "benchmarks/filtering_references.benchmark.txt"
    script:
        "../scripts/reference_filtering.py"





