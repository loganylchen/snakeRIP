#!/usr/bin/env bash

# 重定向标准输出和错误输出到日志文件
exec > "${snakemake_log[0]}" 2>&1

# 判断是否为双端测序
if [ "${snakemake_params[fq_n]}" -eq 2 ]; then
    fastp ${snakemake_params[extra]} \
        -i ${snakemake_input[fq1]} \
        -I ${snakemake_input[fq2]} \
        -o ${snakemake_output[fq1]} \
        -O ${snakemake_output[fq2]} \
        --thread ${snakemake_threads} \
        --html ${snakemake_output[qc_html]} \
        --json ${snakemake_output[qc_json]}
else
    fastp ${snakemake_params[extra]} \
        -i ${snakemake_input[fq1]} \
        -o ${snakemake_output[fq1]} \
        --thread ${snakemake_threads} \
        --html ${snakemake_output[qc_html]} \
        --json ${snakemake_output[qc_json]}
fi

