rule fastp_fastqc_pe:
    input:
        unpack(get_fq)
    output:
        unpack(get_clean_fq),
        qc_html="results/qc/{sample}/{sample}.fastp.html",
        qc_json="results/qc/{sample}/{sample}.fastp.json",
    log:
        "logs/clean_data/{sample}.log"
    benchmark:
        "benchmarks/{sample}.fastp_qc.benchmark.txt"
    params:
        extra=config['params']["fastp"],
        fq_n=get_fq_n,
    container:
        "docker://btrspg/fastp:0.24.0"
    threads:
        config['threads']["fastp"]
    script:
        "../scripts/fastp_clean_data.bash"

