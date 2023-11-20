rule qualimap_rnaseq:
    input:
        bam="results/star/{sample}/{sample}.star.bam",
    output:
        report = "results/qc/{sample}/{sample}_rnaseq.pdf"
    params:
        output_dir = "results/qc/{sample}/",
        reference_gtf="resources/genome.gtf",
        output_file="{sample}_rnaseq.pdf"
    benchmark:
        "benchmarks/{sample}.qualimap_rnaseq.benchmark.txt"
    log:
        "logs/qc/{sample}_qualimap_rnaseq.log"
    conda:
        "../envs/qualimap.yaml"
    resources:
        mem='40G',
        javaopt="-Djava.awt.headless=true -Djava.io.tmpdir=results/tmp/"
    shell:
        "export JAVA_OPTS='{resources.javaopt}' && "
        "qualimap rnaseq -bam {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {params.output_dir} "
        "-outfile {params.output_file} "
        "-outformat PDF:HTML --java-mem-size={resources.mem}> {log}"

rule qualimap_bamqc:
    input:
        bam="results/star/{sample}/{sample}.star.bam",
    output:
        report = "results/qc/{sample}/{sample}_bamqc.pdf"
    params:
        output_dir = "results/qc/{sample}/",
        reference_gtf="resources/genome.gtf",
        output_file="{sample}_bamqc.pdf"
    benchmark:
        "benchmarks/{sample}.qualimap_bamqc.benchmark.txt"
    log:
        "logs/qc/{sample}_qualimap_bamqc.log"
    conda:
        "../envs/qualimap.yaml"
    resources:
        mem='40G',
        javaopt="-Djava.awt.headless=true -Djava.io.tmpdir=results/tmp/"
    shell:
        "export JAVA_OPTS='{resources.javaopt}' && "
        "qualimap bamqc -bam {input.bam} "
        "-gtf {params.reference_gtf} "
        "-outdir {params.output_dir} "
        "-outfile {params.output_file} "
        "-outformat PDF:HTML --java-mem-size={resources.mem}> {log}"