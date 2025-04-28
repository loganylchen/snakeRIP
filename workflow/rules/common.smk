import pandas as pd
import glob
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

project = config['project']

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"Sample_name": str})
    .loc[lambda x: x["Project"] == project]
    .set_index("Sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


def check_raw_data(raw_data_string:str):
    if raw_data_string.endswith('.fq') or raw_data_string.endswith('.fq.gz') or raw_data_string.endswith('.fastq') or raw_data_string.endswith('.fastq.gz'):
        fqs = raw_data_string.split(',')
        if len(fqs) == 1:
            return fqs[0], ''
        elif len(fqs) == 2:
            return fqs[0],fqs[1]
        else:
            raise ValueError(f'{raw_data_string} has more than 2 files')
    else:
        raise ValueError(f'{raw_data_string} is not a valide datatype')


def get_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['Raw_data']
    fq1, fq2 = check_raw_data(raw_data)
    if fq2 == '':
        return {'fq1':fq1}
    else:
        return {'fq1':fq1, 'fq2':fq2}



def get_clean_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['Raw_data']
    fq1, fq2 = check_raw_data(raw_data)
    if fq2 == '':
        return {'fq1':f"{wildcards.project}/clean_data/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
           
        }
    else:
        return {'fq1':f"{wildcards.project}/clean_data/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2':f"{wildcards.project}/clean_data/{wildcards.sample}/{wildcards.sample}_2.fastq.gz",
        }

def get_fq_n(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['Raw_data']
    fq1, fq2 = check_raw_data(raw_data)
    if fq2 == '':
        return 1
    else:
        return 2


def get_final_output():
    final_output = ["resources/genome.fasta","resources/genelist.txt",]
    for sample in samples.index:
        if samples.loc[sample].loc['Sequence_type'] == 'm5C-BS-seq':
            final_output.append(f"{project}/clean_data/{sample}/{sample}_1.fastq.gz")
    return final_output