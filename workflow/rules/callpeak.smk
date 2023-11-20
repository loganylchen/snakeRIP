rule callpeak:
    input:
        treatment=expand("results/star/{rip}/{rip}.star.bam",rip=rip_samples),
        control=expand("results/star/{input}/{input}.star.bam",input=input_samples),
    output:
        multiext(f"callpeak/MACS2_{config['tag']}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    benchmark:
        f"benchmarks/macs2/callpeak_{config['tag']}.txt"
    log:
        f"logs/macs2/callpeak_{config['tag']}.log"
    params:
        config['params']['macs2']
    wrapper:
        "v2.13.0/bio/macs2/callpeak"


rule homer_annotatepeaks:
    input:
        peaks=f"callpeak/MACS2_{config['tag']}_peaks.narrowPeak",
        genome=config['reference']['genome'],
        # optional input files
        gtf=config['reference']['gtf'], # implicitly sets the -gtf flag
        # gene="", # implicitly sets the -gene flag for gene data file to add gene expression or other data types
        # motif_files="peaks_refs/motives.txt", # implicitly sets the -m flag
        # filter_motiv="", # implicitly sets the -fm flag
        # center="",  # implicitly sets the -center flag
        # nearest_peak="peaks_refs/b.peaks", # implicitly sets the -p flag
        # tag="",  # implicitly sets the -d flag for tagDirectories
        # vcf="", # implicitly sets the -vcf flag
        # bed_graph="", # implicitly sets the -bedGraph flag
        # wig="", # implicitly sets the -wig flag
        # map="", # implicitly sets the -map flag
        # cmp_genome="", # implicitly sets the -cmpGenome flag
        # cmp_Liftover="", # implicitly sets the -cmpLiftover flag
        # advanced_annotation=""  # optional, implicitly sets the -ann flag, see http://homer.ucsd.edu/homer/ngs/advancedAnnotation.html
    output:
        annotations=f"callpeak/MACS2_{config['tag']}_peaks.narrowPeak_annot.txt",
        # optional output, implicitly sets the -matrix flag, requires motif_files as input
        matrix=multiext(f"callpeak/MACS2_{config['tag']}_peaks.narrowPeak",
                        ".count.matrix.txt",
                        ".ratio.matrix.txt",
                        ".logPvalue.matrix.txt",
                        ".stats.txt"
                        ),
        # optional output, implicitly sets the -mfasta flag, requires motif_files as input
        mfasta=f"callpeak/MACS2_{config['tag']}_peaks_motif.fasta",
        # # optional output, implicitly sets the -mbed flag, requires motif_files as input
        mbed=f"callpeak/MACS2_{config['tag']}_peaks_motif.bed",
        # # optional output, implicitly sets the -mlogic flag, requires motif_files as input
        mlogic=f"callpeak/MACS2_{config['tag']}_peaks_motif.logic"
    threads: config['threads']['homer']
    benchmark:
        f"benchmarks/homer/callmotif_{config['tag']}.txt"
    params:
        mode="", # add tss, tts or rna mode and options here, i.e. "tss mm8"
        extra="-gid"  # optional params, see http://homer.ucsd.edu/homer/ngs/annotation.html
    log:
        f"logs/annotatePeaks/homer_annotatepeaks_{config['tag']}.log"
    wrapper:
        "v2.13.0/bio/homer/annotatePeaks"