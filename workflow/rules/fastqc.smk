rule qc_fastp:
    input:
        unpack(get_fq),
    output:
        reads = ["{project}/clean_data/{sample}/{sample}_1.fastq.gz"] if get_fq_n == 1 else ["{project}/clean_data/{sample}/{sample}_1.fastq.gz", "{project}/clean_data/{sample}/{sample}_2.fastq.gz"],
        qc_html="{project}/qc/{sample}/{sample}.fastp.html",
        qc_json="{project}/qc/{sample}/{sample}.fastp.json",
    log:
        "logs/{project}/clean_data/{sample}.log"
    benchmark:
        "benchmarks/qc/fastp/{project}_{sample}.benchmark.txt"
    params:
        extra=config['params']["fastp"],
        fq_n=get_fq_n,
    container:
        "docker://btrspg/fastp:0.24.0"
    threads:
        config['threads']["fastp"]
    script:
        "../scripts/fastp_clean_data.sh"

