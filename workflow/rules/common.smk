import pandas as pd
import glob
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")



def check_raw_data(raw_data_string:str):
    fq1, fq2 = raw_data_string.split(',')
    return fq1,fq2

def get_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['raw_data']
    data = check_raw_data(raw_data)
    return {
            'fq1':data[0],
            'fq2':data[1]
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
    final_output.append("callpeak/MACS2_peaks.xls")
    final_output.append("callpeak/MACS2_peaks_motif.fasta")
    return final_output