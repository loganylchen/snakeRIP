import pandas as pd
import glob
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

samples = samples.loc[samples['tag']==config['tag'],:]

validate(samples, schema="../schemas/samples.schema.yaml")



rip_samples = samples.loc[samples['condition']=='RIP',:].index.to_list()
input_samples = samples.loc[samples['condition']=='input',:].index.to_list()



def check_raw_data(raw_data_string:str):
    if raw_data_string.endswith('.fq') or raw_data_string.endswith('.fq.gz') or raw_data_string.endswith('.fastq') or raw_data_string.endswith('.fastq.gz'):
        fq1, fq2 = raw_data_string.split(',')
        return 'fastq',fq1,fq2
    elif raw_data_string.startswith('SRR'):
        return 'sra',raw_data_string, None
    else:
        raise ValueError(f'{raw_data_string} is not a valide datatype')


def get_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['raw_data']
    data_type, *data = check_raw_data(raw_data)
    if data_type == 'fastq':
        return {
                'fq1':data[0],
                'fq2':data[1]
        }
    elif data_type == 'sra':
        return {
            'fq1': f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2': f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_2.fastq.gz"
        }



def get_clean_data(wildcards):
    return {
            'fq1': f"results/clean_fastq/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2': f"results/clean_fastq/{wildcards.sample}/{wildcards.sample}_2.fastq.gz",
            'index':"resources/star_genome",
        }

def get_sra(wildcards):
    return samples.loc[wildcards.sample].loc['raw_data']

def get_final_output():
    final_output = []
    final_output += expand("results/qc/{sample}/{sample}_rnaseq.pdf",sample=samples.index.to_list())
    final_output += expand("results/qc/{sample}/{sample}_bamqc.pdf",sample=samples.index.to_list())
    final_output.append(f"results/callpeak/MACS2_{config['tag']}_peaks.xls")
    final_output.append(f"results/callpeak/MACS2_{config['tag']}_peaks_motif.fasta")
    return final_output