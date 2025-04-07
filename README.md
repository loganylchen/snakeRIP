# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `<description>`


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.

```mermaid
graph TD
    classDef startend fill:#F63B1184,stroke:#222333,stroke-width:2px;
    classDef process fill:#11F63B84,stroke:#222333,stroke-width:2px;
    classDef workflow fill:#E5F6FF,stroke:#73A6FF,stroke-width:2px
    classDef data fill:#FFF6CC,stroke:#FFBC52,stroke-width:2px
    classDef tempdata fill:#FFF111,stroke:#FFBC11,stroke-width:2px
    A(Start):::startend --> B{Sequencing type?}:::process;
    B --is canonical RNA-seq--> DATA_RNA_SEQ(RNA-seq):::data
    B --is canonical DNA-seq--> DATA_DNA_SEQ(WES/WGS):::data
    B --is Glori-seq data --> DATA_GLORI(GLORI-seq):::data
    
    %% B --  --> C[RNA-seq analysis];
    B -- is MERIP-seq --> D[MERIP-seq analysis];
    B -- is ChIP-seq --> E[ChIP-seq analysis];
    B -- is ATAC-seq --> F[ATAC-seq analysis];
    B -- is DNase-seq --> G[DNase-seq analysis];
    B -- is DNase-seq --> H[DNase-seq analysis];
    DATA_RNA_SEQ --> HAMR:::workflow
    DATA_RNA_SEQ --> ModTect:::workflow
    DATA_RNA_SEQ --> Inosine(Inosine Detection):::workflow
    DATA_DNA_SEQ --> Inosine
```