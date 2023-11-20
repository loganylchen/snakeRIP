rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/ref/get-genome.log",
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    cache: True
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["reference"]["species"],
        fmt="gtf",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        flavor="",
    cache: True
    log:
        "logs/ref/get_annotation.log",
    benchmark:
        "benchmarks/get_annotation.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/ref/genome-faidx.log",
    cache: True
    wrapper:
        "v1.21.4/bio/samtools/faidx"

rule create_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/picard/create_dict.log",
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=10240,
    wrapper:
        "v2.3.1/bio/picard/createsequencedictionary"



rule star_index:
    input:
        fasta='resources/genome.fasta',
        gtf='resources/genome.gtf'
    output:
        directory("resources/star_genome"),
        "resources/star_genome/SAindex",
        "resources/star_genome/SA"
    threads: config['threads']['star']
    params:
        sjdb_overhang=100,
        extra="",
    log:
        "logs/star_index_genome.log",
    benchmark:
        "benchmarks/star_index.benchmark.txt"
    wrapper:
        "v1.21.4/bio/star/index"