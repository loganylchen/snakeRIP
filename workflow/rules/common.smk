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



def is_pe(wildcards):
    if samples.loc[wildcards.sample].loc['seq_type']  == 'pe':
        return True
    else:
        return False
def check_raw_data(raw_data_string:str,seq_type:str):
    if raw_data_string.endswith('.fq') or raw_data_string.endswith('.fq.gz') or raw_data_string.endswith('.fastq') or raw_data_string.endswith('.fastq.gz'):
        if seq_type == 'se':
            fq = raw_data_string
            return 'fastq',seq_type,fq,None
        elif seq_type == 'pe':
            fq1, fq2 = raw_data_string.split(',')
            return 'fastq',seq_type,fq1,fq2
        else:
            raise ValueError(f'{seq_type} is not a valide seq_type')
    elif raw_data_string.startswith('SRR'):
        return 'sra',seq_type, raw_data_string, None
    else:
        raise ValueError(f'{raw_data_string} is not a valide datatype')
def get_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['raw_data']
    seq_type = samples.loc[wildcards.sample].loc['seq_type']
    data_type , seq_type,  *data = check_raw_data(raw_data,seq_type)
    if data_type == 'fastq':
        if seq_type =='se':
            return {
                'fq1': data[0]
            }
        elif seq_type =='pe':
            return {
                'fq1':data[0],
                'fq2':data[1]
            }
    else:
        if seq_type == 'se':
            return {
                'fq1':f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}.fastq.gz",
            }
        elif seq_type == 'pe':
            return {
                'fq1':f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
                'fq2':f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_2.fastq.gz"
            }


def get_clean_data(wildcards):
    if samples.loc[wildcards.sample].loc['seq_type'] == 'pe':
        return {
            'fq1': f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2': f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_2.fastq.gz",
            'index':"resources/star_genome",
        }
    elif samples.loc[wildcards.sample].loc['seq_type'] == 'se':
        return {
            'fq1': f"results/clean_fastq/{wildcards.sample}/{wildcards.sample}.fastq.gz",
            'index': "resources/star_genome",
        }
    else:
        raise ValueError(f'{wildcards.sample} is a wired name!')
def get_sra(wildcards):
    return samples.loc[wildcards.sample].loc['raw_data']

def get_final_output():
    contrasts = config['diffexp']['contrasts']
    subclasses = samples.loc[:,config['diffexp']['subclass']].unique()
    final_output = expand("results/star/{sample.sample_name}/ReadsPerGene.out.tab",sample=samples.itertuples())
    final_output += expand("results/hamr/{sample.sample_name}/hamr.mods.txt",sample=samples.itertuples())
    final_output += expand('results/modtect/{sample.sample_name}/modtect.combined.txt',sample=samples.itertuples())
    final_output.append("results/deseq2/count_matrix.rds")
    final_output.append("results/counts/count_matrix.tidy.featureCounts")
    final_output.append("results/wgcna/wgcna.rds")
    for key in contrasts:
        for subclass in subclasses:
            final_output.append(f"results/diffexp/{key}/{subclass}.diffexp.tsv")
            final_output.append(directory(f"results/enrichment/{key}_{subclass}"))
            final_output.append(f"results/visualization/Volcano.{key}_{subclass}.diffexp.pdf")
            final_output.append(f"results/visualization/Volcano.{key}_{subclass}.diffexp.png")
    final_output.append(f"results/visualization/PCA.pdf")
    final_output.append("callpeak/MACS2_peaks.xls")
    # final_output.append("results/counts/all.symbol.tsv")
    return final_output