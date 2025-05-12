#!/bin/bash


exec > "${snakemake_log[0]}" 2>&1

set -x
set -e

genome_reference="${snakemake_input[genome_reference]}"
merged_bam="${snakemake_input[merged_bam]}"
m5c_pileup_tmp="${snakemake_output[m5c_pileup_tmp]}"
noredundance_base="${snakemake_input[noredundance_base]}"
output_file_prefix="${snakemake_params[output_file_prefix]}"
threads="${snakemake[threads]}"

python /opt/conda/RNA-m5C/4_m5C_step-by-step_pileup/pileup_genome_multiprocessing_v1.4.py \
        -P ${threads} \
        -f ${genome_reference} \
        -i ${merged_bam} \
        -o ${m5c_pileup_tmp} 

sort -k 3,3 -k 1,2 -S 5G --parallel ${threads} -T ${m5c_pileup_tmp}.sort_tmp ${m5c_pileup_tmp} > ${m5c_pileup_tmp}.merged.sorted

python /opt/conda/RNA-m5C/4_m5C_step-by-step_pileup/m5C_pileup_formatter.py \
    --db ${noredundance_base} \
    -i ${m5c_pileup_tmp} \
    -o ${m5c_pileup_tmp}.formatted.txt \
    --CR ${m5c_pileup_tmp}.CR.txt

python /opt/conda/RNA-m5C/5_m5C_step-by_step-call_site/m5C_caller.py \
    -i {m5c_pileup_tmp}.formatted.txt \
    -o ${output_file_prefix} \
    --CR overall --method binomial 

#${output_file_prefix}.3.txt  ${output_file_prefix}.None.txt 
