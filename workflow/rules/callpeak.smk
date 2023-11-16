rule callpeak:
    input:
        treatment=expand("results/star/{rip}/{rip}.star.bam",rip=rip_samples),
        control=expand("results/star/{input}/{input}.star.bam",input=input_samples),
    output:
        multiext("callpeak/MACS2",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
    benchmark:
        "benchmarks/macs2/callpeak.txt"
    log:
        "logs/macs2/callpeak.log"
    params:
        config['params']['macs2']
    wrapper:
        "v2.13.0/bio/macs2/callpeak"

