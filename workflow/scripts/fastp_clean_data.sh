#!/usr/bin/env bash

# set -x
# set -e

exec > "${snakemake_log[0]}" 2>&1

reads=(${snakemake_output[reads]})
echo ${reads[0]}


# 判断是否为双端测序
if [ "${snakemake_params[fq_n]}" -eq 2 ]; then
    fastp ${snakemake_params[extra]} \
        -i ${snakemake_input[fq1]} \
        -I ${snakemake_input[fq2]} \
        -o "${reads[0]}" \
        -O "${reads[1]}" \
        --thread ${snakemake[threads]} \
        --html ${snakemake_output[qc_html]} \
        --json ${snakemake_output[qc_json]}
else
    fastp ${snakemake_params[extra]} \
        -i ${snakemake_input[fq1]} \
        -o "${reads[0]}" \
        --thread ${snakemake[threads]} \
        --html ${snakemake_output[qc_html]} \
        --json ${snakemake_output[qc_json]}
fi

